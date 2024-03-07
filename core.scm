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
       (unless (null-list? lines)
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
     (list (length (first lines))
           (map first lines)))))

;; (geno.txt->lmdb "/home/aartaka/git/GEMMA/example/mouse_hs1940.geno.txt" "/tmp/geno-mouse-lmdb/")

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
        (expt-sum 0))
    (for-vec
     (lambda (idx val)
       (set! expt-sum (+ expt-sum (expt (- val mean) 2))))
     vec)
    (/ expt-sum
       (- (vec-length vec) 1))))

(define (cleanup-vec vec)
  (let ((mean (vec-mean vec))
        (var (vec-variance vec)))
    ;; Replace NaNs with mean value.
    (vec-replace-nan vec mean)
    ;; Subtract mean from all the values, "center" them.
    (for-vec
     (lambda (idx val)
       (vec-set! vec idx (- val mean)))
     vec)
    ;; ???
    (for-vec
     (lambda (idx val)
       (vec-set! vec idx
                 (/ val
                    (if (= var 0.0)
                        1
                        (sqrt var)))))
     vec)))

(define (read-genotypes geno.txt lmdb-dir markers individuals)
  (let* ((lines (read-separated-lines geno.txt))
         (mtx (mtx-alloc (length markers) individuals))
         (line-idx 0))
    (mdb:call-with-env-and-txn
     lmdb-dir
     (lambda (env txn)
       (let ((dbi (mdb:dbi-open txn #f 0)))
         (mdb:call-with-cursor
          txn dbi
          (lambda (cursor)
            (mdb:for-cursor
             cursor
             (lambda (key data)
               (let* ((numbers (mdb:val-data-parse
                                data (make-list (/ (mdb:val-size data)
                                                   (sizeof float)) float)))
                      (vec (vec-alloc (length numbers) numbers)))
                 (cleanup-vec vec)
                 (vec->mtx-row! vec mtx line-idx)
                 (set! line-idx (+ 1 line-idx))))))))))
    mtx))

;; (vector-ref (mtx->2d-vector (second (read-geno.txt "/home/aartaka/git/GEMMA/example/BXD_geno.txt"))) 0)
;; (vector-ref (mtx->2d-vector (second (read-geno.txt "/home/aartaka/git/GEMMA/example/mouse_hs1940.geno.txt"))) 0)

(define (kinship file lmdb-dir)
  (let* ((meta (geno.txt->lmdb file lmdb-dir))
         (mtx (read-genotypes file lmdb-dir (length (second meta)) (first meta)))
         (result (mtx-alloc (mtx-columns mtx) (mtx-columns mtx) 0)))
    (dgemm! mtx mtx result #:beta 0 #:transpose-b +trans+)
    result))

;; (kinship "/home/aartaka/git/GEMMA/example/BXD_geno.txt" "/tmp/lmdb-bxd/")



