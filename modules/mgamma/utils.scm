(define-module (mgamma utils)
  #:use-module (srfi srfi-8)
  #:use-module (ice-9 match)
  #:use-module (system foreign)
  #:use-module ((gsl matrices) #:prefix mtx:)
  #:use-module ((gsl vectors) #:prefix vec:)
  #:use-module ((gsl blas) #:prefix blas:)
  #:use-module ((gsl eigensystems) #:prefix eigen:)
  #:use-module ((lapack lapack) #:prefix lapack:)
  #:export-syntax (inc! dec! dotimes dorange with-cleanup)
  #:export (2+
            vec-mean
            vec-replace-nan
            cleanup-mtx
            center-matrix!
            eigendecomposition
            eigendecomposition-zeroed
            submatrix
            submatrix->mtx!))

(define-syntax inc!
  (syntax-rules ()
    ((_ var)
     (set! var (1+ var)))
    ((_ var step)
     (set! var (+ step var)))))
(define-syntax dec!
  (syntax-rules ()
    ((_ var)
     (set! var (1- var)))
    ((_ var step)
     (set! var (- var step )))))

(define (2+ number)
  (+ number 2))

(define-syntax dotimes
  (syntax-rules ()
    ((_ (var upper-bound) body ...)
     (do ((var 0 (1+ var)))
         ((= var upper-bound))
       body ...))))
(define-syntax dorange
  (syntax-rules ()
    ((_ (var lower-bound upper-bound) body ...)
     (do ((var lower-bound (1+ var)))
         ((= var upper-bound))
       body ...))))

(define (vec-mean vec)
  "Mean value of the VEC, with all the NaNs ignored."
  (let ((sum 0)
        (nans 0))
    (vec:for-vec
     (lambda (index value)
       (if (nan? value)
           (set! nans (1+ nans))
           (set! sum (+ value sum))))
     vec)
    (if (or (zero? (vec:length vec))
            (= (vec:length vec) nans))
        0
        (/ sum (- (vec:length vec) nans)))))

(define (vec-replace-nan vec val)
  (vec:for-vec
   (lambda (index value)
     (when (nan? value)
       (vec:set! vec index val)))
   vec))

(define (cleanup-mtx mtx)
  "Plug mean value instead of NaNs, whenever possible.
Also subtract mean from all the values to 'center' them."
  (do ((mtx-rows (mtx:rows mtx))
       (mtx-cols (mtx:columns mtx))
       (row 0 (1+ row)))
      ((= row mtx-rows))
    (match (let rec ((col 0)
                     (nans '())
                     (sum 0))
             (if (= col mtx-cols)
                 (list sum nans)
                 (rec (1+ col)
                      (if (nan? (mtx:get mtx row col))
                          (cons col nans)
                          nans)
                      (if (nan? (mtx:get mtx row col))
                          sum
                          (+ sum (mtx:get mtx row col))))))
      ((sum nans)
       (let ((mean (if (= mtx-cols (length nans))
                       0
                       (/ sum (- mtx-cols (length nans))))))
         (do ((nans nans (cdr nans)))
             ((null? nans))
           (mtx:set! mtx row (car nans) mean))
         (do ((col 0 (1+ col)))
             ((= col mtx-cols))
           (mtx:set! mtx row col (- (mtx:get mtx row col) mean))))))))

(define (center-matrix! mtx)
  "Weird one directly copied from GEMMA's CenterMatrix."
  (let* ((size (mtx:rows mtx))
         (w (vec:alloc size 1.0))
         (gw (blas:gemv mtx w)))
    (blas:syr2! gw w mtx #:alpha (/ -1 size))
    (let ((d (blas:dot w gw)))
      (blas:syr! w mtx #:alpha (/ d (expt size 2))))
    ;; GEMMA says it's transpose, but I don't believe it --aartaka.
    (do ((i 0 (1+ i)))
        ((= i size))
      (do ((j 0 (1+ j)))
          ((= j i))
        (mtx:set! mtx i j (mtx:get mtx j i))))
    (vec:free w gw)
    mtx))

(define (eigendecomposition kinship)
  (let ((evalues-vec (vec:alloc (mtx:rows kinship)))
        (evectors-mtx (mtx:alloc (mtx:rows kinship) (mtx:rows kinship))))
    (lapack:dsyevr
     lapack:+row-major+ lapack:+eigenvalues-and-eigenvectors+ lapack:+range-all+ lapack:+lower+
     (mtx:rows kinship) (mtx:data kinship) (mtx:rows kinship)
     0.0 0.0 0 0 1.0e-7
     ;; Throwaway pointer M.
     (make-c-struct (list int) (list 0))
     (vec:data evalues-vec) (mtx:data evectors-mtx)
     (mtx:rows kinship)
     ;; NOTE: That's what GEMMA defines ISUPPZ like. No idea what this
     ;; means. --aartaka
     (make-c-struct
      (make-list (* 2 (mtx:rows kinship)) int)
      (make-list (* 2 (mtx:rows kinship)) 0)))
    (values evalues-vec evectors-mtx)))

(define (eigendecomposition-zeroed kinship)
  "Eigendecomposition, but zero the values below threshold.
Return two values:
- EVALUES-VEC
- EVECTORS-MTX"
  (receive (evalues-vec evectors-mtx)
      (eigendecomposition kinship)
    (do ((i 0 (1+ i)))
        ((= i (vec:length evalues-vec)))
      (when (< (abs (vec:get evalues-vec i)) 1e-10)
        ;; pylmm uses 1e-6 instead
        (vec:set! evalues-vec i 0)))
    (values evalues-vec evectors-mtx)))

(define-syntax-rule (with-cleanup form cleanup ...)
  (dynamic-wind
    (lambda ()
      #t)
    (lambda ()
      form)
    (lambda ()
      cleanup ...)))

(define (submatrix mtx start-row start-col rows cols)
  "Get a ROWSxCOLS submatrix of MTX, left upper corner at START-ROW, START-COL."
  (let ((new (mtx:alloc rows cols 0)))
    (do ((row 0 (1+ row))
         (src-row start-row (1+ src-row)))
        ((= row rows))
      (do ((col 0 (1+ col))
           (src-col start-col (1+ src-col)))
          ((= col cols))
        (mtx:set! new row col (mtx:get mtx src-row src-col))))
    new))

(define (submatrix->mtx! submatrix mtx mtx-start-row mtx-start-col rows cols)
  "Copy SUBMATRIX back into MXT starting on (mtx-start-row,mtx-start-col)
and copying a ROWSxCOLS chunk."
  (let ((new (mtx:alloc rows cols 0)))
    (do ((row 0 (1+ row))
         (mtx-row mtx-start-row (1+ mtx-row)))
        ((= row rows))
      (do ((col 0 (1+ col))
           (mtx-col mtx-start-col (1+ mtx-col)))
          ((= col cols))
        (mtx:set! mtx mtx-row mtx-col (mtx:get submatrix row col))))))

(define (subvector vec start len)
  (let ((new (vec:alloc len)))
    (do ((idx 0 (1+ idx))
         (vec-idx start (1+ vec-idx)))
        ((= idx len))
      (vec:set! new idx (vec:get vec vec-idx)))
    new))
