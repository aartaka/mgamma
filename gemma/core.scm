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

;; TODO: Make a state machine parser instead? Something like guile-csv
;; TODO: Apply PEG to a whole file?
(define separators-char-set (list->char-set '(#\Tab #\Space #\,)))

(define (string-separate string)
  "Split the string on tab, space, or comma sequences.
Return a list of strings in between these separators."
  (let ((string-len (string-length string)))
    (let rec ((idx 0)
              (start-idx #f))
      (cond
       ((< idx string-len)
        (let* ((char (string-ref string idx))
               (separator? (char-set-contains? separators-char-set char)))
          (cond
           ((and start-idx separator?)
            (cons (substring string start-idx idx)
                  (rec (1+ idx) #f)))
           ((not (or start-idx separator?))
            (rec (1+ idx) idx))
           (else
            (rec (1+ idx) start-idx)))))
       ((and start-idx
             (< start-idx idx))
        (list (substring string start-idx idx)))
       (else
        '())))))

;; For speed.
(define %read-separated-lines-cache (make-hash-table))
(define (read-separated-lines file)
  "Read tab/space/comma-separated fields from FILE.
Return a list of lists of values."
  (if (hash-ref %read-separated-lines-cache file)
      (hash-ref %read-separated-lines-cache file)
      (begin
        (hash-set!
         %read-separated-lines-cache
         file
         (call-with-port (open-input-file file)
           (lambda (port)
             (let read-lines ((line (first (%read-line port))))
               (if (eof-object? line)
                   '()
                   (cons (string-separate line)
                         (read-lines (first (%read-line port)))))))))
        (hash-ref %read-separated-lines-cache file))))

;; (first (read-separated-lines "/home/aartaka/git/GEMMA/example/mouse_hs1940.geno.txt"))

(define (read-bim file)
  (let ((lines (read-separated-lines file)))
    (map (lambda (split)
           (let ( ;; How are morgans/centimorgans marked?
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
           (match split
             ((marker pos chromosome unidentified?)
              (list marker
                    (string->num/nan pos)
                    (string->number chromosome)
                    (string->number unidentified?)))))
         lines)))

;; (read-anno.txt "/home/aartaka/git/GEMMA/example/mouse_hs1940.anno.txt")

(define (geno.txt->lmdb geno.txt-file lmdb-dir)
  "Convert GENO.TXT-FILE to an LMDB-DIR-located database.
Useful to speed up genotype matrix manipulation.
The keys of the resulting DB are marker names.
The values are `double' arrays with one value per individual."
  (let ((lines (read-separated-lines geno.txt-file)))
    (unless (file-exists? (string-append lmdb-dir "data.mdb"))
      (mdb:call-with-wrapped-cursor
       lmdb-dir
       (lambda (env txn dbi cursor)
         (let* ((lines (read-separated-lines geno.txt-file))
                (double-size (sizeof double)))
           (format #t "Keys not found, filling the database at ~s.~%" lmdb-dir)
           (do ((lines lines (cdr lines)))
               ((null? lines))
             (mdb:put txn dbi
                      (mdb:make-val (first (first lines)))
                      (let* ((values (cdddr (first lines)))
                             (values-len (length values)))
                        (mdb:make-val (make-c-struct
                                       (make-list values-len double)
                                       (map string->num/nan values))
                                      (* values-len double-size)))
                      mdb:+noodupdata+))))
       #:mapsize (* 40 10485760)))
    (list (- (length (car lines)) 3)
          (map car lines))))

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

(define (cleanup-vector vec)
  (let ((mean (vec-mean vec)))
    (vec-replace-nan vec mean)
    (vec:add-constant! vec (- mean))))

(define (string-na? str)
  (string= "NA" str))

(define (useful-individuals pheno.txt)
  (let ((lines (read-separated-lines pheno.txt)))
    (do ((lines lines (cdr lines))
         (inds (list (not (string-na? (caar lines))))
               (cons (not (string-na? (caar lines)))
                     inds)))
        ((null? lines)
         (reverse inds)))))

(define (count-false lst)
  (do ((lst lst (cdr lst))
       (false-count
        (if (car lst)
            0
            1)
        (+ false-count
           (if (car lst)
               0
               1))))
      ((null? lst) false-count)))

(define* (useful-snps geno.txt pheno.txt #:key (miss-level 0.05) (maf-level 0.01))
  (let* ((lines (read-separated-lines geno.txt))
         (useful-inds (useful-individuals pheno.txt))
         (useful-snp-table (make-hash-table)))
    (do ((lines lines (cdr lines)))
        ((null? lines))
      (let* ((row (car lines))
             (name (first row))
             (inds (map string->number (cdddr row)))
             (ind-count (length inds))
             (maf (do ((inds inds (cdr inds))
                       (useful-inds useful-inds
                                    (cdr useful-inds))
                       (maf (if (car useful-inds)
                                (car inds)
                                0)
                            (+ maf (if (car useful-inds)
                                       (car inds)
                                       0))))
                      ((null? inds) maf)))
             (miss-count (count-false inds))
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
          (let* ((vec (vec:alloc individuals 0)))
            (memcpy (vec:ptr vec 0) (mdb:val-data data) (mdb:val-size data))
            (cleanup-vector vec)
            (mtx:vec->row! vec mtx line-idx)
            (set! line-idx (1+ line-idx))))))
     #:mapsize (* 40 10485760))
    mtx))

(define (kinship mtx n-useful-snps)
  "Calculate the kinship matrix for genotype MTX."
  (let ((result (mtx:alloc (mtx:columns mtx) (mtx:columns mtx) 0)))
    (dgemm! mtx mtx result #:beta 0 #:transpose-a +trans+)
    (mtx:scale! result (/ 1 n-useful-snps))
    result))

(define (kinship->cxx.txt kinship-mtx cxx.txt)
  (let ((last-row 0))
    (call-with-port (open-output-file cxx.txt)
      (lambda (port)
        (mtx:for-each
         (lambda (row column value)
           (when (not (= row last-row))
             (newline port))
           (format port "~f " value))
         kinship-mtx)))))

(define (cxx.txt->kinship cxx.txt)
  (let* ((lines (read-separated-lines cxx.txt))
         (lines-len (length lines))
         (line-len (length (first lines)))
         (result (mtx:alloc lines-len line-len)))
    (do ((row-idx 0 (1+ row-idx))
         (row lines (cdr lines)))
        ((= row-idx lines-len))
      (do ((column-idx 0 (1+ column-idx))
           (column (first row) (cdr column)))
          ((= column-idx line-len))
        (mtx:set! result row-idx column-idx (string->number (first column)))))
    result))

;; (cxx.txt->kinship "/home/aartaka/git/GEMMA/output/BXD_mouse.cXX.txt")

(define (kmain geno.txt pheno.txt lmdb-dir)
  (let* ((meta (geno.txt->lmdb geno.txt lmdb-dir))
         (useful-snps (useful-snps geno.txt pheno.txt))
         (mtx (lmdb->genotypes-mtx lmdb-dir (second meta) (first meta))))
    (kinship mtx (hash-count (cut or #t <> <>) useful-snps))))

;; (define kin (kmain "/home/aartaka/git/GEMMA/example/mouse_hs1940.geno.txt" "/home/aartaka/git/GEMMA/example/mouse_hs1940.pheno.txt" "/tmp/lmdb-hs/"))
;; (define lmdb-dir "/tmp/lmdb-bxd/")
;; (define geno.txt "/home/aartaka/git/GEMMA/example/BXD_geno.txt")
;; (define pheno.txt "/home/aartaka/git/GEMMA/example/BXD_pheno.txt")
;; (define meta (geno.txt->lmdb geno.txt lmdb-dir))
;; (define useful-snps-table (useful-snps geno.txt pheno.txt))
;; (define mtx (lmdb->genotypes-mtx lmdb-dir (second meta) (first meta)))
;; (define kin (kinship mtx (hash-count (cut or #t <> <>) useful-snps-table)))
