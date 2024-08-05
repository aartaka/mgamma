(define-module (mgamma io)
  #:use-module (mgamma config)
  #:use-module (rnrs bytevectors)
  #:use-module (srfi srfi-1)
  #:use-module (ice-9 ports)
  #:use-module (ice-9 textual-ports)
  #:use-module (ice-9 rdelim)
  #:use-module (ice-9 match)
  #:use-module (system foreign)
  #:use-module ((lmdb lmdb) #:prefix mdb:)
  #:use-module ((gsl vectors) #:prefix vec:)
  #:use-module ((gsl matrices) #:prefix mtx:)
  #:use-module (json)
  #:export (geno.txt->lmdb
            geno.txt->genotypes-mtx
            lmdb->genotypes-mtx
            kinship->cxx.txt
            kinship->lmdb
            lmdb->kinship
            cxx.txt->kinship
            pheno.txt->pheno-mtx
            covariates.txt->cvt-mtx
            eigenu.txt->eigenvectors
            snp-params->assoc.txt))

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
(define (read-separated-lines file)
  "Read tab/space/comma-separated fields from FILE.
Return a list of lists of values.
Necessary because the biological data often comes in BIMBAM format, a
frivolous mix of tab and comma separated values."
  (call-with-port (open-input-file file)
    (lambda (port)
      (remove
       null?
       (let read-lines ((line (first (%read-line port))))
         (if (eof-object? line)
             '()
             (cons (string-separate line)
                   (read-lines (first (%read-line port))))))))))

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

(define (geno.txt->lmdb geno.txt-file lmdb-dir)
  "Convert GENO.TXT-FILE to an LMDB-DIR-located database.
Useful to speed up genotype matrix manipulation.
The keys of the resulting DB are marker names.
The values are `float' arrays (converted to `double' matrices when
read back into Scheme) with one value per individual.
DB also contains a \"meta\" record specifying the \"type\",
\"version\", and whether the data if \"float\", as serialized JSON
object."
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

(define (pointer=? ptr1 ptr2 size)
  (bytevector=? (pointer->bytevector ptr1 size)
                (pointer->bytevector ptr2 size)))

(define* (lmdb->genotypes-mtx lmdb-dir)
  "Read the data from LMDB-DIR and convert it to GSL matrix.
Useful because LMDB-DIR-resident data often comes as floats, not
doubles."
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
            (when (and (not mtx)
                       (not meta?))
              ;; Markers x individuals
              (set! mtx (mtx:alloc
                         ;; Meta record is mandatory, I guess?
                         (1- (mdb:stat-entries (mdb:dbi-stat txn dbi)))
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

(define (kinship->lmdb kinship-mtx lmdb-dir)
  "Dump KINSHIP-MTX to an LMDB residing in LMDB-DIR."
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
  "Read kinship matrix from LMDB-DIR."
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
  "Dump KINSHIP-MTX to CXX.TXT file."
  (let ((last-row 0))
    (call-with-port (open-output-file cxx.txt)
      (lambda (port)
        (mtx:for-each
         (lambda (row column value)
           (when (not (= row last-row))
             (newline port)
             (set! last-row row))
           (format port "~f\t " value))
         kinship-mtx)))))

(define (txt->mtx file.txt)
  "Generic text-to-matrix reading.
Return a `mtx:mtx?' of doubles/f64."
  (let* ((lines (read-separated-lines file.txt))
         (lines-len (length lines))
         (line-len (length (first lines)))
         (result (mtx:alloc lines-len line-len)))
    (do ((row-idx 0 (1+ row-idx))
         (row lines (cdr row)))
        ((= row-idx lines-len))
      (do ((column-idx 0 (1+ column-idx))
           (column (first row) (cdr column)))
          ((= column-idx line-len))
        (mtx:set! result row-idx column-idx (string->num/nan (first column)))))
    result))

(define (cxx.txt->kinship cxx.txt)
  (txt->mtx cxx.txt))

(define (pheno.txt->pheno-mtx pheno.txt)
  (txt->mtx pheno.txt))

(define (covariates.txt->cvt-mtx covariates.txt)
  (txt->mtx covariates.txt))

(define (eigenu.txt->eigenvectors eigenu.txt)
  "Read eigenU (matrix of eigenvectors) from the EIGENU.TXT file."
  (txt->mtx eigenu.txt))

(define (snp-params->assoc.txt params-table assoc.txt)
  "Dump PARAMS-TABLE to the ASSOC.TXT file."
  (call-with-port (open-output-file assoc.txt)
    (lambda (p)
      (format p "rs\t maf\t beta\t se\t logl_H1\t l_remle\t p_wald~%")
      (hash-map->list
       (lambda (key value)
         (match value
           ((maf
             beta se
             lambda-remle
             p-wald
             logl-alt)
            (format p "~a\t ~s\t ~s\t ~s\t ~s\t ~s\t ~s~%"
                    key     maf  beta se   logl-alt lambda-remle p-wald))))
       params-table))))
