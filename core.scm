(define-module (gemma core)
  #:use-module (ice-9 ports)
  #:use-module (ice-9 textual-ports)
  #:use-module (ice-9 rdelim)
  #:use-module (ice-9 match)
  #:use-module (srfi srfi-1)
  #:use-module (srfi srfi-26)
  #:use-module (srfi srfi-43)
  #:use-module ((gsl matrices) #:prefix mtx:)
  #:use-module ((gsl vectors) #:prefix vec:)
  #:use-module (gsl blas)
  #:use-module (system foreign)
  #:use-module (system foreign-library)
  #:use-module (rnrs bytevectors)
  #:use-module ((lmdb lmdb) #:prefix mdb:))

;; For speed.
(define %read-separated-lines-cache (make-hash-table))
(define (read-separated-lines file)
  "Read tab/space/comma-separated fields from FILE.
Return a list of lists of values."
  (if (hash-ref %read-separated-lines-cache file)
      (hash-ref %read-separated-lines-cache file)
      (hash-set!
       %read-separated-lines-cache
       file
       (call-with-port (open-input-file file)
         (lambda (port)
           (let read-lines ((line (first (%read-line port))))
             (if (eof-object? line)
                 '()
                 (cons (remove string-null?
                               (string-split
                                line (lambda (c) (memq c '(#\Tab #\Space #\,)))))
                       (read-lines (first (%read-line port)))))))))))

;; (first (read-separated-lines "/home/aartaka/git/GEMMA/example/mouse_hs1940.geno.txt"))

(define (read-bim file)
  (let ((lines (read-separated-lines file)))
    (map (lambda (split)
           (let (;; How are morgans/centimorgans marked?
                 (position (string->number (third split)))
                 (base-pair-coordinate (string->number (fourth split))))
             `(,(first split)
               ,(second split)
               ,position
               ,base-pair-coordinate
               ,@(cddddr split))))
         lines)))

;; (read-bim "/home/aartaka/git/GEMMA/example/mouse_hs1940.bim")

(define (string->num/nan str)
  "Convert the string to a valid float of NaN."
  (if (string= "NA" str)
      +nan.0
      (string->number str)))

(define (read-anno.txt file)
  (let ((lines (read-separated-lines file)))
    (map (lambda (split)
           (let ((marker (first split))
                 (pos (second split))
                 (chromosome (third split))
                 (unidentified? (fourth split)))
             (list marker
                   (string->num/nan pos)
                   (string->number chromosome)
                   (string->number unidentified?))))
         lines)))

;; (read-anno.txt "/home/aartaka/git/GEMMA/example/mouse_hs1940.anno.txt")

(define (geno.txt->lmdb geno.txt-file lmdb-dir)
  "Convert GENO.TXT-FILE to a n LMDB-DIR-located database.
Useful to speed up genotype matrix manipulation.
The keys of the resulting DB are marker names.
The values are `double' arrays with one value per individual."
  (mdb:call-with-env-and-txn
   lmdb-dir
   (lambda (env txn)
     (let* ((dbi (mdb:dbi-open txn #f 0))
            (lines (read-separated-lines geno.txt-file))
            (cursor (mdb:cursor-open txn dbi)))
       (with-exception-handler
           ;; Weird logic, but here we are: if there's a
           ;; "no-key-found" exception on CURSOR, fill the DB.
           (lambda (exc)
             (format #t "Keys not found, filling the database at ~s." lmdb-dir)
             (do ((lines lines (cdr lines)))
                 ((null? lines))
               (mdb:put txn dbi
                        (mdb:make-val (first (first lines)))
                        (let ((values (cdddr (first lines))))
                          (mdb:make-val (make-c-struct
                                         (make-list (length values) double)
                                         (map string->num/nan values))
                                        (* (length values)
                                           (sizeof double))))
                        mdb:+noodupdata+)))
         (lambda ()
           (mdb:cursor-first cursor)))
       (mdb:cursor-close cursor)
       (list (- (length (first lines)) 3)
             (map first lines))))
   #:mapsize (* 40 10485760)))

;; (geno.txt->lmdb "/home/aartaka/git/GEMMA/example/BXD_geno.txt" "/tmp/geno-mouse-lmdb/")

(define (vec-mean vec)
  (/ (vec:sum vec)
     (vec:length vec)))

(define (vec-replace-nan vec val)
  (vec:for-vec
   (lambda (index value)
     (when (nan? value)
       (vec:set! vec index val)))
   vec))

(define (string-na? str)
  (string=? "NA" str))

(define (useful-individuals pheno.txt)
  (let ((lines (read-separated-lines pheno.txt)))
    (do ((lines lines (cdr lines))
         (inds (list (not (string-na? (caar lines))))
               (cons (not (string-na? (caar lines)))
                     inds)))
        ((null? lines)
         (reverse inds)))))

(define* (useful-snps geno.txt pheno.txt #:key (miss-level 0.05) (maf-level 0.01))
  (let ((lines (read-separated-lines geno.txt))
        (useful-inds (useful-individuals pheno.txt))
        (useful-snp-table (make-hash-table)))
    (do ((snp 0 (1+ snp)))
        ((= snp (length lines)))
      (let* ((row (list-ref lines snp))
             (name (first row))
             (inds (cdddr row))
             (ind-count (length inds))
             (maf (do ((inds (cdddr row)
                             (cdr inds))
                       (useful-inds useful-inds
                                    (cdr useful-inds))
                       (maf (if (not (car useful-inds))
                                0
                                (string->number (car inds)))
                            (+ maf (if (not (car useful-inds))
                                       0
                                       (string->number (car inds))))))
                      ((null? inds) maf)))
             (miss-count (count string-na? inds))
             (maf (/ maf (* 2 (- ind-count miss-count)))))
        (when (and (< (/ miss-count ind-count) miss-level)
                   (< maf-level maf (- 1 maf-level)))
          (hash-set! useful-snp-table name #t))))
    useful-snp-table))

(define memcpy
  (foreign-library-function
   #f "memcpy"
   #:return-type '*
   #:arg-types (list '* '* size_t)))

(define (string-control? string)
  (string-any (lambda (c)
                (eq? 'Cc (char-general-category c)))
              string))

(define (cleanup-vector vec)
  (let ((mean (vec-mean vec)))
    (vec-replace-nan vec mean)
    (vec:add-constant! vec (- mean))))

(define (lmdb->genotypes-mtx lmdb-dir markers individuals)
  "Read the data from LMDB-DIR and convert it to GSL matrix.
The resulting matrix is #MARKERSxINDIVIDUALS sized."
  (let* ((mtx (mtx:alloc (length markers) individuals))
         (line-idx 0))
    (mdb:call-with-wrapped-cursor
     lmdb-dir
     (lambda (env txn dbi cursor)
       (mdb:for-cursor
        cursor
        (lambda (key data)
          ;; FIXME: It sometimes happens that LMDB table has one
          ;; or two corrupted rows. Ignoring them here
          (unless (string-control? (mdb:val-data-string key))
            (let* ((vec (vec:alloc individuals 0)))
              (memcpy (vec:ptr vec 0) (mdb:val-data data) (mdb:val-size data))
              (cleanup-vector vec)
              (mtx:vec->row! vec mtx line-idx)
              (set! line-idx (+ 1 line-idx)))))))
     #:mapsize (* 40 10485760))
    mtx))

;; (vector-ref (mtx->2d-vector (second (read-geno.txt "/home/aartaka/git/GEMMA/example/BXD_geno.txt"))) 0)
;; (vector-ref (mtx->2d-vector (second (read-geno.txt "/home/aartaka/git/GEMMA/example/mouse_hs1940.geno.txt"))) 0)

(define (kinship mtx n-useful-snps)
  "Calculate the kinship matrix for genotype MTX."
  (let ((result (mtx:alloc (mtx:columns mtx) (mtx:columns mtx) 0)))
    (dgemm! mtx mtx result #:beta 0 #:transpose-a +trans+)
    (mtx:scale! result (/ 1 n-useful-snps))
    result))
(define (kmain geno.txt pheno.txt lmdb-dir)
  (let* ((meta (geno.txt->lmdb geno.txt lmdb-dir))
         (useful-snps (useful-snps geno.txt pheno.txt))
         (mtx (lmdb->genotypes-mtx lmdb-dir (second meta) (first meta))))
    (kinship mtx (hash-count (cut or #t <> <>) useful-snps))))

;; (define kin (kmain "/home/aartaka/git/GEMMA/example/BXD_geno.txt" "/home/aartaka/git/GEMMA/example/BXD_pheno.txt" "/tmp/lmdb-bxd/"))
;; (define lmdb-dir "/tmp/lmdb-bxd/")
;; (define geno.txt "/home/aartaka/git/GEMMA/example/BXD_geno.txt")
;; (define pheno.txt "/home/aartaka/git/GEMMA/example/BXD_pheno.txt")
;; (define meta (geno.txt->lmdb geno.txt lmdb-dir))
;; (define mtx (lmdb->genotypes-mtx lmdb-dir (second meta) (first meta)))
;; (define kin (kinship mtx))
