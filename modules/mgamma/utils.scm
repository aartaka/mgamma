(define-module (mgamma utils)
  #:use-module (ice-9 match)
  #:use-module (srfi srfi-1)
  #:use-module (srfi srfi-8)
  #:use-module (system foreign)
  #:use-module (system foreign-library)
  ;; #:use-module ((gsl core) #:prefix gsl:)
  #:use-module ((gsl matrices) #:prefix mtx:)
  #:use-module ((gsl vectors) #:prefix vec:)
  #:use-module ((gsl blas) #:prefix blas:)
  #:use-module ((gsl eigensystems) #:prefix eigen:)
  #:use-module ((lapack lapack) #:prefix lapack:)
  #:use-module (rnrs base)
  #:export-syntax (inc! dec!
                        approximately-zero?
                        dotimes dorange
                        with-cleanup with-gsl-free
                        define-parameterized)
  #:export (2+
            vec-mean
            vec-replace-nan
            cleanup-mtx
            center-matrix!
            eigendecomposition
            eigendecomposition-zeroed
            ;; cblas-dgemm
            submatrix
            submatrix->mtx!
            subvector
            subvector->vec!
            gsl-free))

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

(define (approximately-zero? x)
  (< (abs x) 1e-10))

(define-syntax dotimes
  (syntax-rules ()
    ((_ (var upper-bound . result-forms) body ...)
     (do ((var 0 (1+ var))
          (bound upper-bound))
         ((= var bound)
          . result-forms)
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
  (let* ((size (mtx:rows mtx)))
    (with-gsl-free
     ((w (vec:alloc size 1.0))
      (gw (blas:gemv mtx w)))
     (blas:syr2! gw w mtx #:alpha (/ -1 size))
     (blas:syr! w mtx #:alpha (/ (blas:dot w gw) (* size size)))
     ;; GEMMA says it's transpose, but I don't believe it --aartaka.
     (dotimes (i size)
       (dotimes (j i)
         (mtx:set! mtx i j (mtx:get mtx j i))))
     mtx)))


(define libopenblas (load-foreign-library "libopenblas.so"))
;; (define libopenblas (load-foreign-library "/gnu/store/vyglbwdr92nkglxzzbq9anagrby1g3g3-openblas-debug-0.3.20/lib/libopenblas.so" #:global? #t))
(define dsyevr-
  (foreign-library-function
   libopenblas "dsyevr_"
   #:arg-types '(* * * *
                   * * * *
                   * * * * *
                   * * * *
                   * * * *)
   #:return-type void))

(define (dsyevr jobz range uplo n a lda vl vu il iu abstol m w z ldz isuppz work lwork iwork liwork)
  (let  ((int->ptr
          (lambda (i)
            (assert (integer? i))
            (make-c-struct (list int) (list i))))
         (char->ptr
          (lambda (c)
            (assert (char? c))
            (make-c-struct (list uint8) (list (char->integer c)))))
         (inexact->ptr
          (lambda (f)
            (assert (real? f))
            (make-c-struct (list double) (list f))))
         (info (make-c-struct (list int) '(0))))
    (dsyevr- (char->ptr jobz) (char->ptr range) (char->ptr uplo)
             (int->ptr n) a (int->ptr lda)
             (inexact->ptr vl) (inexact->ptr vu)
             (int->ptr il) (int->ptr iu)
             (inexact->ptr abstol) m w z
             (int->ptr ldz) isuppz work lwork iwork liwork info)
    (unless (zero? (first (parse-c-struct info (list int))))
      (error 'dsyevr "Failed" (first (parse-c-struct info (list int)))))))

(load-extension "libmgamma.so" "init_eigendecomp")
;; (load-extension "/home/aartaka/git/mgamma/extension/libmgamma.so" "init_eigendecomp")

(define (eigendecomposition kinship)
  (receive (vals vecs status)
      (eigendecomp (char->integer #\V) (char->integer #\A) (char->integer #\L)
                   (mtx:rows kinship) (mtx:data kinship)
                   0.0 0.0 0 0 1.0e-7)
    (let ((evalues-vec (vec:alloc (mtx:rows kinship) vals))
          (evectors-mtx (mtx:alloc-square (mtx:rows kinship) vecs)))
      ;; FIXME: Doesn't give the right results yet—some signs are
      ;; inverted. Sigh. --aartaka
      (eigen:sort! evalues-vec evectors-mtx #:ascending)
      (mtx:transpose! evectors-mtx #f)
      (values evalues-vec evectors-mtx))))

(define (eigendecomposition-zeroed kinship)
  "Eigendecomposition, but zero the values below threshold.
Return two values:
- EVALUES-VEC
- EVECTORS-MTX"
  (receive (evalues-vec evectors-mtx)
      (eigendecomposition kinship)
    (vec:for-each
     (lambda (i val)
       (when (< val 1e-10)
         ;; pylmm uses 1e-6 instead
         (vec:set! evalues-vec i 0)))
     evalues-vec)
    (values evalues-vec evectors-mtx)))

;; (define %cblas-dgemm
;;   (foreign-library-function
;;    gsl:libgslcblas "cblas_dgemm"
;;    #:return-type void
;;    #:arg-types (list int int int
;;                      int int int
;;                      double '* int
;;                      '* int
;;                      double '* int)))

;; (define (cblas-dgemm transpose-a transpose-b alpha a b beta c)
;;   (assert (and (mtx:mtx? a) (mtx:mtx? b) (mtx:mtx? c)))
;;   (%cblas-dgemm
;;    blas:+row-major+
;;    (if transpose-a
;;        blas:+transpose+
;;        blas:+no-transpose+)
;;    (if transpose-b
;;        blas:+transpose+
;;        blas:+no-transpose+)
;;    (mtx:rows c) (mtx:columns c)
;;    (if transpose-a
;;        (mtx:rows a)
;;        (mtx:columns a))
;;    alpha (mtx:data a) (third (mtx:parts a))
;;    (mtx:data b) (third (mtx:parts b))
;;    beta (mtx:data c) (third (mtx:parts c))))

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
(define (subvector->vec! sub vec start len)
  (do ((idx 0 (1+ idx))
       (vec-idx start))
      ((= idx len))
    (vec:set! vec vec-idx (vec:get sub idx))))

(define (gsl-free . things)
  (for-each (lambda (thing)
              (if (mtx:mtx? thing)
                  (mtx:free thing)
                  (vec:free thing)))
            things))

(define-syntax-rule (with-gsl-free ((var init) ...) body ...)
  (let* ((var init) ...)
    (with-cleanup
     (begin body ...)
     (gsl-free var ...))))

;; Define a NAMEd procedure and PARAMETER-NAMEd parameter
;; variable. Proceeds with running BODY with ARGS when PARAMETER-NAMEd
;; variable is #false. When PARAMETER-NAME is `parameterize'd to a new
;; procedure, call this procedure on ARGS instead. Useful to override
;; a procedure (like replacing the results for testing or providing
;; shortcut data for long-running computation.)
(define-syntax define-parameterized
  (syntax-rules ()
    ((_ ((name parameter-name) . args) body ...)
     (begin
       (define parameter-name (make-parameter #f))
       (define (name . rest)
         (apply (or (parameter-name)
                    (lambda args body ...))
                rest))))
    ((_ (name . args) body ...)
     (define (name . args) body ...))
    ((_ name value)
     (define name (make-parameter value)))))
