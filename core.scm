(define-module (gemma core)
  #:use-module (ice-9 ports)
  #:use-module (ice-9 textual-ports)
  #:use-module (ice-9 rdelim)
  #:use-module (srfi srfi-1)
  #:use-module (srfi srfi-43)
  #:use-module (gsl matrices)
  #:use-module (gsl vectors)
  #:use-module (gsl blas)
  #:use-module (system foreign)
  #:use-module (rnrs bytevectors)
  #:use-module ((lmdb lmdb) #:prefix mdb:))

(define (read-separated-lines file)
  (call-with-port (open-input-file file)
    (lambda (port)
      (let read-lines ((line (first (%read-line port))))
        (if (eof-object? line)
            '()
            (cons (remove (lambda (s)
                            (equal? "" s))
                          (string-split
                           line (lambda (c) (memq c '(#\Tab #\Space #\,)))))
                  (read-lines (first (%read-line port)))))))))

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
  (mdb:with-env-and-txn
   (lmdb-dir #:mapsize (* 40 10485760))
   (env txn)
   (let ((dbi (mdb:dbi-open txn #f 0))
         (lines (read-separated-lines geno.txt-file)))
     (let rec ((lines lines))
       (unless (null? lines)
         (mdb:put txn dbi
                  (mdb:make-val (first (first lines)))
                  (let ((values (cdddr (first lines))))
                    (mdb:make-val (make-c-struct
                                   (make-list (length values) float)
                                   (map string->num/nan values))
                                  (* (length values)
                                     (sizeof float))))
                  mdb:+noodupdata+)
         (rec (cdr lines))))
     (list (- (length (first lines)) 3)
           (map first lines)))))

;; (geno.txt->lmdb "/home/aartaka/git/GEMMA/example/BXD_geno.txt" "/tmp/geno-mouse-lmdb/")

(define (vec-mean vec)
  (let ((sum 0)
        (nans 0))
    (for-vec
     (lambda (index value)
       (if (nan? value)
           (set! nans (+ 1 nans))
           (set! sum (+ sum value))))
     vec)
    (if (= nans (vec-length vec))
        0
        (/ sum
           (- (vec-length vec) nans)))))

(define (vec-replace-nan vec val)
  (for-vec
   (lambda (index value)
     (when (nan? value)
       (vec-set! vec index val)))
   vec))

(define-public (vec-variance vec)
  (let ((mean (vec-mean vec))
        (expt-sum 0)
        (nans 0))
    (for-vec
     (lambda (idx val)
       (if (nan? val)
           (set! nans (+ 1 nans))
           (set! expt-sum (+ expt-sum (expt (- val mean) 2)))))
     vec)
    (/ expt-sum
       (- (vec-length vec) nans 1))))


;; (vector-ref (mtx->2d-vector (second (read-geno.txt "/home/aartaka/git/GEMMA/example/BXD_geno.txt"))) 0)

;; Get the line/vector
;; Find mean
;; Find variance
;; Plug mean in place of NA
;; Subtract the mean
;; Scale by 1/sqrt(geno_var)???
;; Multiply the matrix
(define (cleanup-vec vec)
  (let ((mean (vec-mean vec))
        (var (vec-variance vec)))
    ;; Replace NaNs with mean value.
    (vec-replace-nan vec mean)
    ;; Subtract mean from all the values, "center" them.
    (vec-add-constant! vec (- mean))
    ;; ??? -gk 2?
    ;; (vec-scale! vec (if (= var 0.0)
    ;;                     1
    ;;                     (/ 1 (sqrt var))))
    ))

(define (read-genotypes lmdb-dir markers individuals)
  (let* ((mtx (mtx-alloc (length markers) individuals))
         (line-idx 0))
    (mdb:call-with-env-and-txn
     lmdb-dir
     (lambda (env txn)
       (let ((dbi (mdb:dbi-open txn #f 0)))
         (mdb:with-cursor
          txn dbi (cursor)
          (mdb:for-cursor
           cursor
           (lambda (key data)
             (let* ((numbers (mdb:val-data-parse
                              data (make-list (/ (mdb:val-size data)
                                                 (sizeof float)) float)))
                    ;; FIXME: It sometimes happen that LMDB table has
                    ;; one or two corrupted rows. Ignoring them here
                    (vec (if (every (lambda (c)
                                      (not (eq? 'Cc (char-general-category c))))
                                    (string->list (mdb:val-data-string key)))
                             (vec-alloc (length numbers) numbers)
                             #f)))
               (when vec
                 (cleanup-vec vec)
                 (vec->mtx-row! vec mtx line-idx)
                 (set! line-idx (+ 1 line-idx))))))))))
    mtx))

;; (vector-ref (mtx->2d-vector (second (read-geno.txt "/home/aartaka/git/GEMMA/example/BXD_geno.txt"))) 0)
;; (vector-ref (mtx->2d-vector (second (read-geno.txt "/home/aartaka/git/GEMMA/example/mouse_hs1940.geno.txt"))) 0)

(define (kinship file lmdb-dir)
  (let* ((meta (geno.txt->lmdb file lmdb-dir))
         (mtx (read-genotypes lmdb-dir (second meta) (first meta)))
         (result (mtx-alloc (mtx-columns mtx) (mtx-columns mtx))))
    (dgemm! mtx mtx result #:beta 0 #:transpose-a +trans+)
    (mtx-scale! result (/ 1 (length (second meta))))
    result))

;; (let ((meta (geno.txt->lmdb "/home/aartaka/git/GEMMA/example/BXD_geno.txt" "/tmp/lmdb-bxd/")))
;;   (read-genotypes "/tmp/lmdb-bxd/" (second meta) (first meta)))
(mtx-get (kinship "/home/aartaka/git/GEMMA/example/BXD_geno.txt" "/tmp/lmdb-bxd/") 0 0)
