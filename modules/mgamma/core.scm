(define-module (mgamma core)
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
  #:use-module ((lmdb lmdb) #:prefix mdb:)
  #:use-module (json)
  #:export (mapsize
            geno.txt->lmdb
            geno.txt->genotypes-mtx
            lmdb->genotypes-mtx
            kinship-mtx
            kinship->cxx.txt
            kinship->lmdb
            lmdb->kinship
            cxx.txt->kinship
            useful-snps
            useful-individuals
            cleanup-mtx))

(define mapsize (make-parameter (* 100 10485760)))

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
  (let ((lines (read-separated-lines geno.txt-file))
        (float-size (sizeof float)))
    (unless (file-exists? (string-append lmdb-dir "data.mdb"))
      (mdb:call-with-wrapped-cursor
       lmdb-dir #f
       (lambda (env txn dbi cursor)
         (format #t "Keys not found, filling the database at ~s.~%" lmdb-dir)
         (mdb:put txn dbi
                  (mdb:make-val (string->pointer "meta" "UTF-8") 4)
                  (scm->json-string `(("type" . "geno")
                                      ("version" . "0.0.1")
                                      ("float" . #t))))
         (do ((lines lines (cdr lines)))
             ((null? lines))
           (mdb:put txn dbi
                    (mdb:make-val (first (first lines)))
                    (let* ((values (cdddr (first lines)))
                           (values-len (length values)))
                      (mdb:make-val (make-c-struct
                                     (make-list values-len float)
                                     (map string->num/nan values))
                                    (* values-len float-size)))
                    mdb:+noodupdata+)))
       #:mapsize (mapsize)))
    (list (- (length (car lines)) 3)
          (map car lines))))

;; (geno.txt->lmdb "/home/aartaka/git/GEMMA/example/BXD_geno.txt" "/tmp/geno-mouse-lmdb/")

(define (vec-mean vec)
  (let ((sum 0)
        (nans 0))
    (vec:for-vec
     (lambda (index value)
       (if (nan? value)
           (set! nans (1+ nans))
           (set! sum (+ value sum))))
     vec)
    (/ sum (- (vec:length vec) nans))))

(define (vec-replace-nan vec val)
  (vec:for-vec
   (lambda (index value)
     (when (nan? value)
       (vec:set! vec index val)))
   vec))

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

(define* (useful-snps genotypes-mtx markers pheno.txt #:key (miss-level 0.05) (maf-level 0.01))
  (let* ((useful-inds (useful-individuals pheno.txt))
         (ind-count (length useful-inds))
         (useful-snp-table (make-hash-table))
         (mtx-rows (mtx:rows genotypes-mtx))
         (mtx-cols (mtx:columns genotypes-mtx)))
    (do ((row 0 (1+ row))
         (markers markers (cdr markers)))
        ((= row mtx-rows))
      (let* ((name (car markers))
             (ind 0) ;; To allow using it in do for maf
             (maf (do ((ind ind (1+ ind))
                       (useful-inds useful-inds (cdr useful-inds))
                       (maf (if (car useful-inds)
                                (mtx:get genotypes-mtx row ind)
                                0)
                            (if (car useful-inds)
                                (+ maf (mtx:get genotypes-mtx row ind))
                                maf)))
                      ((= ind mtx-cols) maf)))
             (miss-count (do ((ind 0 (1+ ind))
                              (nans 0 (if (nan? (mtx:get genotypes-mtx row ind))
                                          (1+ nans)
                                          nans)))
                             ((= ind mtx-cols) nans)))
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

(define (pointer=? ptr1 ptr2 size)
  (bytevector=? (pointer->bytevector ptr1 size)
                (pointer->bytevector ptr2 size)))

(define* (lmdb->genotypes-mtx lmdb-dir)
  "Read the data from LMDB-DIR and convert it to GSL matrix."
  (let* ((mtx #f)
         (tmp-vec #f)
         (line-idx 0)
         (markers (list))
         (float-size (sizeof float))
         (meta #f))
    (mdb:call-with-wrapped-cursor
     lmdb-dir #f
     (lambda (env txn dbi cursor)
       (mdb:for-cursor
        cursor
        (lambda (key data)
          (let ((individuals (/ (mdb:val-size data) float-size))
                (meta? (pointer=? (string->pointer "meta" "UTF-8") (mdb:val-data key) 4)))
            (unless mtx
              ;; Markers x individuals
              (set! mtx (mtx:alloc (mdb:stat-entries (mdb:dbi-stat txn dbi))
                                   individuals)))
            (if meta?
                (set! meta (json-string->scm (mdb:val-data-string data)))
                (do ((i 0 (1+ i))
                     (inds (mdb:val-data-parse data (make-list individuals float))
                           (cdr inds)))
                    ((= i individuals))
                  (mtx:set! mtx line-idx i (car inds))))
            (unless meta?
              (set! markers (cons (mdb:val-data-string key) markers))
              (set! line-idx (1+ line-idx)))))))
     #:mapsize (mapsize))
    (when tmp-vec
      (vec:free tmp-vec))
    (list (or mtx (mtx:alloc 0 0 0))
          (reverse! markers)
          meta)))

(define (geno.txt->genotypes-mtx geno.txt)
  "Convert GENO.TXT-FILE to a proper genotype matrix
Return a (MATRIX MARKER-NAMES) list."
  (let* ((lines (read-separated-lines geno.txt))
         (mtx (mtx:alloc (length lines)
                         ;; The first 3 lines: marker, chr, and chr.
                         (- (length (first lines)) 3) +nan.0)))
    (do ((row 0 (1+ row))
         (lines lines (cdr lines)))
        ((null? lines))
      (do ((col 0 (1+ col))
           ;; cdddar: the elements starting from the third of the
           ;; current line.
           (inds (cdddar lines) (cdr inds)))
          ((null? inds))
        (let ((ind (string->number (first inds))))
          (when ind
            (mtx:set! mtx row col ind)))))
    (list mtx (map first lines))))

(define (cleanup-mtx mtx)
  (do ((mtx-rows (mtx:rows mtx))
       (mtx-cols (mtx:columns mtx))
       (row 0 (1+ row)))
      ((= row mtx-rows))
    (match (do ((col 0 (1+ col))
                (is-nan? (nan? (mtx:get mtx row 0))
                         (nan? (mtx:get mtx row col)))
                (sum 0 (if is-nan?
                           sum
                           (+ sum (mtx:get mtx row col))))
                (nans '() (if is-nan?
                              (cons col nans)
                              nans)))
               ((= col mtx-cols)
                (list sum nans)))
      ((sum nans)
       (let ((mean (/ sum (- mtx-cols (length nans)))))
         (do ((nans nans (cdr nans)))
             ((null? nans))
           (mtx:set! mtx row (car nans) mean))
         (do ((col 0 (1+ col)))
             ((= col mtx-cols))
           (mtx:set! mtx row col (- (mtx:get mtx row col) mean))))))))

(define (kinship-mtx mtx n-useful-snps)
  "Calculate the kinship matrix for genotype MTX."
  (let ((result (mtx:alloc (mtx:columns mtx) (mtx:columns mtx) 0)))
    (cleanup-mtx mtx)
    (dgemm! mtx mtx result #:beta 0 #:transpose-a +trans+)
    (mtx:scale! result (/ 1 n-useful-snps))
    result))

(define (kinship->lmdb kinship-mtx lmdb-dir)
  (let ((env (mdb:env-create #:mapsize (mapsize)))
        (float-size (sizeof float))
        (uint-size (sizeof unsigned-int))
        (rows (mtx:rows kinship-mtx)))
    (mdb:env-open env lmdb-dir)
    (do ((row-step 1 (1+ row-step))
         (initial-row 0 (+ initial-row row-step)))
        ((> initial-row rows))
      (let* ((txn (mdb:txn-begin env))
             (dbi (mdb:dbi-open txn)))
        (when (zero? (mdb:stat-entries (mdb:dbi-stat txn dbi)))
          (mdb:put txn dbi
                   (mdb:make-val (string->pointer "meta" "UTF-8") 4)
                   (scm->json-string `(("type" . "GRM")
                                       ("version" . "0.0.1")
                                       ("float" . #t)
                                       ("symmetric" . #t)))))
        (do ((row initial-row (1+ row))
             (row-size (* float-size (- rows initial-row))
                       (- row-size float-size)))
            ((or (= row rows)
                 (= row (+ initial-row row-step))))
          (mdb:put!
           txn dbi
           (mdb:make-val (make-c-struct (list unsigned-int) (list row))
                         uint-size)
           (mdb:make-val
            (make-c-struct
             (make-list (- rows row) float)
             (let rec ((col row))
               (if (= col rows)
                   '()
                   (cons (mtx:get kinship-mtx row col)
                         (rec (1+ col))))))
            row-size)))
        (mdb:txn-commit txn)))))

(define (lmdb->kinship lmdb-dir)
  (let ((mtx #f))
    (mdb:with-wrapped-cursor
     (lmdb-dir #f #:mapsize (mapsize))
     (env txn dbi cursor)
     (mdb:for-cursor
      cursor
      (lambda (key value)
        (unless mtx
          ;; 1- because there's meta record in addition to data.
          (set! mtx (mtx:alloc (1- (mdb:stat-entries (mdb:dbi-stat txn dbi)))
                               (1- (mdb:stat-entries (mdb:dbi-stat txn dbi)))
                               +nan.0)))
        (unless (pointer=? (string->pointer "meta" "UTF-8") (mdb:val-data key) 4)
          (let ((row (first (mdb:val-data-parse key (list unsigned-int)))))
            (do ((col row (1+ col))
                 (data (mdb:val-data-parse
                        value
                        (make-list (- (mtx:rows mtx) row) float))
                       (cdr data)))
                ((= col (mtx:columns mtx)))
              (mtx:set! mtx row col (car data))))))))
    (mtx:for-each
     (lambda (row column value)
       (when (nan? value)
         (mtx:set! mtx row column (mtx:get mtx column row))))
     mtx)
    mtx))

(define (kinship->cxx.txt kinship-mtx cxx.txt)
  (let ((last-row 0))
    (call-with-port (open-output-file cxx.txt)
      (lambda (port)
        (mtx:for-each
         (lambda (row column value)
           (when (not (= row last-row))
             (newline port)
             (set! last-row row))
           (format port "~f " value))
         kinship-mtx)))))

(define (cxx.txt->kinship cxx.txt)
  (let* ((lines (read-separated-lines cxx.txt))
         (lines-len (length lines))
         (line-len (length (first lines)))
         (result (mtx:alloc lines-len line-len)))
    (do ((row-idx 0 (1+ row-idx))
         (row lines (cdr row)))
        ((= row-idx lines-len))
      (do ((column-idx 0 (1+ column-idx))
           (column (first row) (cdr column)))
          ((= column-idx line-len))
        (mtx:set! result row-idx column-idx (string->number (first column)))))
    result))

;; (cxx.txt->kinship "/home/aartaka/git/GEMMA/output/BXD_mouse.cXX.txt")
