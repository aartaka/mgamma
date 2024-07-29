(define-module (mgamma mvlmm)
  #:use-module (mgamma lmm)
  #:use-module (mgamma utils)
  #:use-module (mgamma config)
  #:use-module (srfi srfi-1)
  #:use-module (srfi srfi-8)
  #:use-module (srfi srfi-26)
  #:use-module (system foreign)
  #:use-module (system foreign-library)
  #:use-module ((gsl core) #:prefix gsl:)
  #:use-module ((gsl vectors) #:prefix vec:)
  #:use-module ((gsl matrices) #:prefix mtx:)
  #:use-module ((gsl blas) #:prefix blas:)
  #:use-module ((gsl linear-algebra) #:prefix linalg:)
  #:export (mvlmm-analyze))

(define (eigen-proc vg ve)
  "EigenProc"
  (let ((d-size (mtx:rows vg)))
    (mtx:with
     (ve-temp (mtx:rows ve) (mtx:columns ve) ve)
     (mtx:with
      (ve-h d-size d-size 0)
      (mtx:with
       (ve-hi d-size d-size 0)
       (let ((logdet-ve 0))
         (receive (ul dl)
             (eigendecomposition-zeroed ve-temp)
           (do ((i 0 (1+ i)))
               ((= i d-size))
             (let ((d (vec:get dl i)))
               (when (positive? d)
                 (set! logdet-ve (+ logdet-ve (log d)))
                 (mtx:with-column
                  (u-col ul i)
                  (set! d (sqrt d))
                  (blas:syr! u-col ve-h #:alpha d)
                  (blas:syr! u-col ve-hi #:alpha (/ 1 d))))))
           ;; Copied from `center-matrix!'
           (do ((i 0 (1+ i)))
               ((= i d-size))
             (do ((j 0 (1+ j)))
                 ((= j i))
               (mtx:set! ve-h i j (mtx:get ve-h j i))
               (mtx:set! ve-hi i j (mtx:get ve-hi j i))))
           (mtx:free ul)
           (vec:free dl))
         (let* ((vg-ve-hi (blas:gemm vg ve-hi))
                (lam (blas:gemm ve-hi vg-ve-hi)))
           (receive (ul dl)
               (eigendecomposition-zeroed lam)
             (vec:for-each
              (lambda (idx val)
                (when (negative? val)
                  (vec:set! dl idx 0)))
              dl)
             (let ((ultveh (blas:gemm ul ve-h #:transpose-a #t))
                   (ultvehi (blas:gemm ul ve-hi #:transpose-a #t)))
               (values dl ultveh ultvehi logdet-ve))))))))))

(define (calc-qi eval dl x)
  "CalcQi"
  (let* ((n-size (vec:length eval))
         (d-size (vec:length dl))
         (c-size (mtx:rows x))
         (dc-size (* d-size (1+ c-size)))
         (q (mtx:alloc dc-size dc-size))
         (qi (mtx:copy! q)))
    (do ((i 0 (1+ i)))
        ((= i c-size))
      (do ((j 0 (1+ j)))
          ((= j c-size))
        (do ((l 0 (1+ l)))
            ((= j d-size))
          (mtx:set!
           q (+ l (* i d-size)) (+ l (* j d-size))
           (if (< j i)
               (mtx:get q (+ l (* j d-size)) (+ l (* i d-size)))
               ;; Yes, it's horrible, but C++ version is even
               ;; worse --aartaka
               (let x-sum ((k 0)
                           (d 0.0))
                 (if (= k n-size)
                     d
                     (x-sum
                      (1+ k)
                      (+ d (/ (* (mtx:get x i k)
                                 (mtx:get x j k))
                              (1+ (* (vec:get dl l) (vec:get eval k)))))))))))))
    (receive (lu perms signum)
        (linalg:decompose q)
      (linalg:%invert lu perms qi)
      ;; FIXME: Leaks PERMS.
      (values qi (linalg:%determinant-log q)))))

(define (mph-initial eval x y)
  "MphInitial"
  (let* ((n-size (vec:length eval))
         (c-size (mtx:rows x))
         (d-size (mtx:rows y))
         (xt (mtx:transpose! x #t))
         (y-row-tmp (vec:alloc (mtx:columns y)))
         (vg (mtx:alloc d-size d-size))
         (ve (mtx:alloc d-size d-size))
         (b-sub (mtx:alloc d-size c-size)))
    ;; FIXME: GEMMA's MphInitial B is modified by submatrices, but
    ;; here it's modified in full.
    (do ((i 0 (1+ i)))
        ((= i d-size))
      (mtx:row->vec! y i y-row-tmp)
      (receive (lam logl-h0)
          (calc-lambda-null xt y-row-tmp eval)
        (receive (g e b)
            (calc-vg-ve-beta eval xt y-row-tmp lam)
          (mtx:set! vg i i g)
          (mtx:set! ve i i e))))
    ;; FIXME: Ignoring the d-size > 4 case for now
    (receive (dl ultveh ultvehi)
        (eigen-proc vg ve)
      (receive (qi logdet-q)
          (calc-qi eval dl x)
        (let ((utvehiy (blas:gemm ultvehi y))
              (xhiy (vec:calloc (* d-size c-size))))
          (do ((i 0 (1+ i)))
              ((= i d-size))
            (do ((j 0 (1+ j))
                 (d 0.0))
                ((= j c-size))
              (set! d 0.0)
              (do ((k 0 (1+ k)))
                  ((= j n-size))
                (set! d (+ d
                           (/ (* (mtx:get x j k)
                                 (mtx:get utvehiy i k))
                              (1+ (* (vec:get eval k)
                                     (vec:get dl i)))))))
              (vec:set! xhiy (+ i (* j d-size)) d)))
          (let* ((beta (blas:gemv qi xhiy))
                 (sub-beta (lambda (offset size)
                             (let ((new (vec:alloc size)))
                               (do ((i 0 (1+ i)))
                                   ((= i size))
                                 (vec:set! new i (vec:get beta (+ i offset))))))))
            (do ((i 0 (1+ i)))
                ((= i c-size))
              (let* ((sub (sub-beta (* i d-size) d-size))
                     (ultvehsub (blas:gemv ultveh sub)))
                (mtx:vec->column! ultvehsub b-sub i)
                (vec:free sub ultvehsub))))
          (vec:free dl xhiy)
          (mtx:free ultveh ultvehi qi))))
    ;; TODO: Ensure everything is freed properly.
    (values vg ve b-sub)))

(define (calc-x-hi-y eval dl x ultvehiy)
  "CalcXHiY"
  (let* ((n-size (vec:length eval))
         (c-size (mtx:rows x))
         (d-size (vec:length dl))
         (xhiy (vec:alloc (* d-size c-size))))
    (do ((i 0 (1+ i)))
        ((= i d-size))
      (let ((dl (vec:get dl i)))
        (do ((j 0 (1+ j)))
            ((= j c-size))
          (vec:set!
           xhiy (+ i (* j d-size))
           (let calc-d ((k 0)
                        (d 0))
             (if (= k n-size)
                 d
                 (calc-d (1+ k)
                         (+ d (/ (* (mtx:get x j k)
                                    (mtx:get ultvehiy i k))
                                 (1+ (* dl (vec:get eval k))))))))))))
    xhiy))

(define (mph-calc-logl eval xhiy dl ultvehiy qi)
  (let ((n-size (vec:length eval))
        (d-size (vec:length dl))
        (dc-size (mtx:rows qi))
        (logl 0))
    (do ((k 0 (1+ k)))
        ((= k n-size))
      (do ((i 0 (1+ i)))
          ((= i d-size))
        (let ((d (1+ (* (vec:get eval k)
                        (vec:get dl i)))))
          (set! logl (+ logl
                        (/ (expt (mtx:get ultvehiy i k) 2)
                           d)
                        (log d))))))
    (let* ((qiv (blas:gemv qi xhiy))
           (d (blas:dot xhiy qiv)))
      (vec:free qiv)
      (* -1/2 (- logl d)))))

(define (calc-omega eval dl)
  "CalcOmega"
  (let* ((n-size (vec:length eval))
         (d-size (vec:length dl))
         (omega-u (mtx:alloc d-size n-size))
         (omega-e (mtx:alloc d-size n-size)))
    (do ((k 0 (1+ k)))
        ((= k n-size))
      (let ((delta (vec:get eval k)))
        (do ((i 0 (1+ i)))
            ((= i d-size))
          (let* ((dl (vec:get dl i))
                 (du (/ dl (1+ (* delta dl))))
                 (de (* delta du)))
            (mtx:set! omega-u i k du)
            (mtx:set! omega-e i k de)))))
    (values omega-u omega-e)))

(define (update-rl-b xhiy qi ultvehib)
  "UpdateRL_B"
  (let* ((d-size (mtx:rows ultvehib))
         (c-size (mtx:columns ultvehib))
         (dc-size (mtx:rows qi))
         (b (blas:gemv qi xhiy)))
    (do ((i 0 (1+ i)))
        ((= i c-size))
      (mtx:with-column
       (ultvehib-col ultvehib i)
       ;; We don't have a luxury of vector views GEMMA abuses.
       ;; gsl_vector_view UltVehiB_col = gsl_matrix_column(UltVehiB, i);
       ;; gsl_vector_const_view b_subcol =
       ;;     gsl_vector_const_subvector(b, i * d_size, d_size);
       ;; gsl_vector_memcpy(&UltVehiB_col.vector, &b_subcol.vector);
       (vec:for-each
        (lambda (idx val)
          (vec:set! ultvehib-col idx (vec:get b (+ (* i d-size) idx))))
        ultvehib-col)
       (mtx:vec->column! ultvehib-col ultvehib i))
      (vec:free b)
      ultvehib)))

(define (update-u omega-e ultvehiy ultvehibx ultvehiu)
  "UpdateU"
  (mtx:copy! ultvehiy ultvehiu)
  (mtx:subtract! ultvehiu ultvehibx)
  (mtx:multiply-elements! ultvehiu omega-e))

(define (update-e ultvehiy ultvehibx ultvehiu ultvehie)
  "UpdateE"
  (mtx:copy! ultvehiy ultvehie)
  (mtx:subtract! ultvehie ultvehibx)
  (mtx:subtract! ultvehie ultvehiu))

(define (update-l-b x xxti ultvehiy ultvehiu ultvehibx ultvehib)
  "UpdateL_B"
  (let ((c-size (mtx:rows x))
        (d-size (mtx:rows ultvehiy)))
    (mtx:with
     (yux d-size c-size)
     (mtx:copy! ultvehiy ultvehibx)
     (mtx:subtract! ultvehibx ultvehiy)
     (blas:gemm! ultvehibx x yux #:beta 0 #:transpose-b #t)
     (blas:gemm! yux xxti ultvehib #:beta 0)
     ultvehib)))

(define (calc-sigma reml? eval dl x omega-u omega-e ultveh qi)
  "CalcSigma"
  (let* ((n-size (vec:length eval))
         (c-size (mtx:rows x))
         (d-size (vec:length dl))
         (dc-size (mtx:rows qi))
         (sigma-uu (mtx:calloc d-size d-size))
         (sigma-ee (mtx:calloc d-size d-size)))
    (dotimes (k n-size)
      (let ((fill (lambda (omega sigma)
                    (mtx:with-column
                     (omega-col omega k)
                     (vec:for-each
                      (lambda (idx val)
                        (mtx:set! sigma k k (+ (mtx:get sigma k k)
                                               val)))
                      omega-col)))))
        (fill omega-u sigma-uu)
        (fill omega-e sigma-ee)))
    (when reml?
      (mtx:with
       (mu dc-size d-size 0)
       (mtx:with
        (me dc-size d-size 0)
        (mtx:with
         (qim dc-size d-size)
         (dotimes (k n-size)
           (let ((delta (vec:get eval k)))
             (dotimes (i d-size)
               (let ((dl (vec:get dl i)))
                 (dotimes (j c-size)
                   (let ((d (/ (mtx:get x j k)
                               (1+ (* delta dl)))))
                     (mtx:set! me (+ i (* j d-size)) i d)
                     (mtx:set! mu (+ i (* j d-size)) i (* d dl))))))
             (blas:gemm! qi mu qim #:beta 0)
             (blas:gemm! mu qim sigma-uu #:transpose-a #t)
             (blas:gemm! qi me qim #:beta 0)
             (blas:gemm! me qim sigma-ee #:transpose-a #t)))))))
    (mtx:with
     (m d-size d-size)
     (blas:gemm! sigma-uu ultveh m #:beta 0)
     (blas:gemm! ultveh m sigma-uu #:transpose-a #t #:beta 0)
     (blas:gemm! sigma-ee ultveh m #:beta 0)
     (blas:gemm! ultveh m sigma-ee #:transpose-a #t #:beta 0))
    (values sigma-uu sigma-ee)))

(define (update-v eval u sigma-uu sigma-ee vg ve)
  "UpdateV"
  (mtx:fill! vg 0)
  (mtx:fill! ve 0)
  (let ((n-size (vec:length eval))
        (d-size (mtx:rows u)))
    (dotimes (k n-size)
      (let ((delta (vec:get eval k)))
        (unless (zero? delta)
          (mtx:with-column
           (u-col u k)
           (blas:dsyr! u-col vg #:alpha (/ 1 delta))))))
    (dotimes (i d-size)
      (dotimes (j i)
        (mtx:set! vg i j (mtx:get vg j i))
        (mtx:set! ve i j (mtx:get ve j i))))
    (mtx:add! vg sigma-uu)
    (mtx:add! ve sigma-ee)
    (mtx:scale! vg (/ 1 n-size))
    (mtx:scale! ve (/ 1 n-size))))

(define (mph-em reml? eval x y vg ve b)
  "MphEM"
  (let* ((n-size (vec:length eval))
         (c-size (mtx:rows x))
         (d-size (mtx:rows y))
         (dc-size (* d-size c-size))
         (xxt (blas:syrk x))
         (xxti (mtx:copy! xxt))
         (ultvehib (mtx:alloc d-size c-size))
         (ultvehibx (mtx:alloc d-size n-size))
         (ultvehiu (mtx:alloc d-size n-size))
         (ultvehie (mtx:alloc d-size n-size))
         (ultvehiy (mtx:alloc d-size n-size)))
    (do ((i 0 (1+ i)))
        ((= i c-size))
      (do ((j 0 (1+ j)))
          ((= j i))
        (mtx:set! xxt i j (mtx:get xxt j i))))
    (let ((logl-const (receive (lu pmt sig)
                          (linalg:decompose xxt)
                        (linalg:%invert lu pmt xxti)
                        (if reml?
                            (+ (* -1/2 (- n-size c-size) d-size (log (* 2 +pi+)))
                               (* 1/2 d-size (linalg:%determinant-log lu)))
                            (* -1/2 n-size d-size (log (* 2 +pi+))))))
          (logl-new #f))
      ;; TODO: Precision check
      (dotimes (t (em-iter))
        (receive (dl ultveh ultvehi logdet-ve)
            (eigen-proc vg ve)
          (receive (qi logdet-q)
              (calc-qi eval dl x)
            (blas:gemm! ultvehi y ultvehiy #:beta 0)
            (let* ((xhiy (calc-x-hi-y eval dl x ultvehiy)))
              (set! logl-new
                (+ logl-const
                   (mph-calc-logl eval xhiy dl ultvehiy qi)
                   ;; - 0.5 * (double)n_size * logdet_Ve;
                   (* -1/2 n-size logdet-ve)
                   (if reml?
                       (* -1/2 (- logdet-q (* c-size logdet-ve)))
                       0)))
              (receive (omega-u omega-e)
                  (calc-omega eval dl)
                (cond
                 (reml?
                  (update-rl-b xhiy qi ultvehib)
                  (blas:gemm! ultvehib x ultvehibx #:beta 0))
                 ((zero? t)
                  (blas:gemm! ultvehi b ultvehib #:beta 0)
                  (blas:gemm! ultvehib x ultvehibx #:beta 0)))
                (update-u omega-e ultvehiy ultvehibx ultvehiu)
                (unless reml?
                  (update-l-b x xxti ultvehiy ultvehiu ultvehibx ultvehib)
                  (blas:gemm! ultvehib x ultvehibx #:beta 0))
                (update-e ultvehiy ultvehibx ultvehiu ultvehie)
                (let* ((u-hat (blas:gemm ultveh ultvehiu #:transpose-a #t))
                       (e-hat (blas:gemm ultveh ultvehie #:transpose-a #t))
                       (b (blas:gemm ultveh ultvehib #:transpose-a #t)))
                  (receive (sigma-uu sigma-ee)
                      (calc-sigma reml? eval dl x omega-u omega-e ultveh qi)
                    (update-v eval u-hat sigma-uu sigma-ee vg ve)
                    (mtx:free ultveh ultvehi ultvehiy qi
                              omega-u omega-e u-hat e-hat b sigma-uu sigma-uu)
                    (vec:free dl xhiy))))))))
      (mtx:free xxt xxti ultvehib ultvehibx ultvehiu ultvehie ultvehiy)
      (values vg ve logl-new))))

(define (getindex i j d-size)
  "GetIndex"
  (let ((s (if (< j i) j i))
        (l (if (< j i) i j)))
    (+ (* s 1/2 (- (* 2 d-size) (1+ s)))
       l
       (- s))))

(define (update-vg-ve hessian gradient scale vg ve)
  "UpdateVgVe"
  (let ((v-size (/ (vec:length gradient) 2))
        (d-size (mtx:rows vg)))
    (vec:with
     (vec-v (* 2 v-size))
     (dotimes (i d-size)
       (dorange
        (j i d-size)
        (let ((v (getindex i j d-size)))
          (vec:set! vec-v v (mtx:get vg i j))
          (vec:set! vec-v (+ v v-size) (mtx:get ve i j)))))
     (blas:gemv! hessian gradient vec-v #:alpha (* -1 scale))
     (dotimes (i d-size)
       (dorange (j i d-size)
         (let* ((v (getindex i j d-size))
                (dg (vec:get vec-v v))
                (de (vec:get vec-v (+ v v-size))))
           (mtx:set! vg i j dg)
           (mtx:set! vg j i dg)
           (mtx:set! ve i j de)
           (mtx:set! ve j i de)))))))

(define (calc-hi-qi eval x vg ve)
  "CalcHiQi"
  (let* ((n-size (vec:length eval))
         (c-size (mtx:rows x))
         (d-size (mtx:rows vg))
         (dc-size (* d-size c-size))
         (hi-all (mtx:calloc 2 (* 2 n-size)))
         (qi (mtx:calloc dc-size dc-size)))
    (mtx:with
     (mat-dd d-size d-size)
     (receive (dl ultveh ultvehi logdet-ve)
         (eigen-proc vg ve)
       (let ((logdet-h (* n-size logdet-ve)))
         (dotimes (k n-size)
           (let ((delta (vec:get eval k)))
             (mtx:copy! ultvehi mat-dd)
             (dotimes (i d-size)
               (let* ((dl (vec:get dl i))
                      (d (1+ (* delta dl))))
                 (mtx:with-row
                  (mat-row mat-dd i)
                  (vec:scale! mat-row (/ 1 d))
                  (mtx:vec->row! mat-row mat-dd i))
                 (set! logdet-h (+ logdet-h (log d)))))
             (with-gsl-free
              ((hi-k (blas:gemm ultvehi mat-dd #:transpose-a #t)))
              ;; FIXME: is d-size always 2? hi-all only has 2 rows...
              (dotimes (row d-size)
                (do ((hi-k-column 0 (1+ hi-k-column))
                     (hi-all-column (* k d-size) (1+ hi-all-column)))
                    ((= hi-all-column (* (1+ k) d-size)))
                  (mtx:set! hi-all row hi-all-column
                            (mtx:get hi-k row hi-k-column)))))))
         (receive (qi logdet-q)
             (calc-qi eval dl x)
           (mtx:free qi)
           (let ((logdet-q (- logdet-q (* c-size logdet-h))))
             (dotimes (i c-size)
               (dotimes (j c-size)
                 (let ((qi-sub (submatrix
                                qi (* i d-size) (* j d-size)
                                d-size d-size)))
                   (if (< j i)
                       (let ((qi-sym (submatrix
                                      qi (* j d-size) (* i d-size)
                                      d-size d-size)))
                         (mtx:transpose! qi-sub qi-sym)
                         (submatrix->mtx!
                          qi-sym qi
                          (* j d-size) (* i d-size) d-size d-size))
                       (begin
                         (blas:gemm! qi-sub ultveh mat-dd #:beta 0)
                         (blas:gemm! ultveh mat-dd qi-sub #:beta 0)
                         (submatrix->mtx!
                          qi-sub qi
                          (* i d-size) (* j d-size) d-size d-size))))))
             (mtx:free ultveh ultvehi)
             (vec:free dl)
             (values hi-all qi logdet-h logdet-q))))))))

(define (calc-hiy-all y hi-all)
  "Calc_Hiy_all"
  (let* ((n-size (mtx:columns y))
         (d-size (mtx:rows y))
         (hiy-all (mtx:calloc d-size n-size)))
    (dotimes (k n-size)
      (with-gsl-free
       ((hi-k (submatrix hi-all 0 (* k d-size) d-size d-size)))
       (mtx:with-column
        (y-k y k)
        (mtx:with-column
         (hiy-k hiy-all k)
         (blas:gemv! hi-k y-k hiy-k #:beta 0)
         (mtx:vec->column! hiy-k hiy-all k)))))
    hiy-all))

(define (calc-xhi-all x hi-all)
  "Calc_xHi_all"
  (let* ((n-size (mtx:columns x))
         (c-size (mtx:rows x))
         (d-size (mtx:rows hi-all))
         (xhi-all (mtx:calloc (* 2 c-size) (* 2 n-size))))
    (dotimes (k n-size)
      (let ((hi-k (submatrix hi-all 0 (* k d-size) d-size d-size)))
        (dotimes (i c-size)
          (let ((xhi-sub
                 (submatrix
                  xhi-all
                  (* i d-size) (* k d-size) d-size d-size)))
            (mtx:copy! hi-k xhi-sub)
            (mtx:scale! xhi-sub (mtx:get x i k))
            (submatrix->mtx! xhi-sub xhi-all
                             (* i d-size) (* k d-size) d-size d-size)))))
    xhi-all))

(define (calc-xhiy y xhi-all)
  "Calc_xHiy"
  (let* ((n-size (mtx:columns y))
         (d-size (mtx:rows y))
         (dc-size (mtx:rows xhi-all))
         (xhiy (vec:calloc dc-size)))
    (dotimes (k n-size)
      (with-gsl-free
       ((xhi-k (submatrix xhi-all 0 (* k d-size) dc-size d-size))
        (y-k (mtx:column->vec! y k)))
       (blas:gemv! xhi-k y-k xhiy)))
    xhiy))

(define (calc-yhiy y hiy-all)
  "Calc_yHiy"
  (let ((n-size (mtx:columns y))
        (yhiy 0))
    (dotimes (k n-size)
      (mtx:with-column
       (y-k y k)
       (mtx:with-column
        (hiy-k hiy-all k)
        (set! yhiy (+ yhiy (blas:dot hiy-k y-k))))))
    yhiy))

(define (calc-xhidhiy eval xhi hiy i j)
  "Calc_xHiDHiy"
  (let* ((dc-size (mtx:rows xhi))
         (n-size (vec:length eval))
         (d-size (mtx:rows hiy))
         (xhidhiy-g (vec:calloc dc-size))
         (xhidhiy-e (vec:calloc dc-size)))
    (dotimes (k n-size)
      (mtx:with-column
       (xhi-col-i xhi (+ i (* k d-size)))
       (let* ((delta (vec:get eval k))
              (d (mtx:get hiy j k)))
         (blas:daxpy! (* d delta) xhi-col-i xhidhiy-g)
         (blas:daxpy! d xhi-col-i xhidhiy-e)
         (unless (= i j)
           (mtx:with-column
            (xhi-col-j xhi (+ j (* k d-size)))
            (let ((d (mtx:get hiy i k)))
              (blas:daxpy! (* d delta) xhi-col-j xhidhiy-g)
              (blas:daxpy! d xhi-col-j xhidhiy-e)))))))
    (values xhidhiy-g xhidhiy-e)))

(define (calc-xhidhiy-all eval xhi hiy)
  "Calc_xHiDHiy_all"
  (let* ((d-size (mtx:rows hiy))
         (dc-size (mtx:rows xhi))
         (v-size (* d-size (1+ d-size) 1/2))
         (xhidhiy-all-g (mtx:calloc dc-size v-size))
         (xhidhiy-all-e (mtx:calloc dc-size v-size)))
    (dotimes (i d-size)
      (dorange
       (j i d-size)
       (let ((v (getindex i j d-size)))
         (receive (xhidhiy-g xhidhiy-e)
             (calc-xhidhiy eval xhi hiy i j)
           (mtx:vec->column! xhidhiy-g xhidhiy-all-g v)
           (mtx:vec->column! xhidhiy-e xhidhiy-all-e v)
           (vec:free xhidhiy-g xhidhiy-e)))))
    (values xhidhiy-all-g xhidhiy-all-e)))

(define (calc-xhidhix eval xhi i j)
  "Calc_xHiDHix"
  (let* ((dc-size (mtx:rows xhi))
         (n-size (vec:length eval))
         (d-size (/ (mtx:columns xhi) n-size))
         (xhidhix-g (mtx:calloc dc-size dc-size))
         (xhidhix-e (mtx:calloc dc-size dc-size)))
    (with-gsl-free
     ((mat-dcdc (mtx:alloc-square dc-size))
      (mat-dcdc-t (mtx:alloc-square dc-size)))
     (dotimes (k n-size)
       (let ((delta (vec:get eval k)))
         (mtx:with-column
          (xhi-col-i xhi (+ i (* k d-size)))
          (mtx:with-column
           (xhi-col-j xhi (+ j (* k d-size)))
           (mtx:fill! mat-dcdc 0)
           (blas:ger! xhi-col-i xhi-col-i mat-dcdc)
           (mtx:transpose! mat-dcdc mat-dcdc-t)
           (mtx:add! xhidhix-e mat-dcdc)
           (mtx:scale! mat-dcdc delta)
           (mtx:add! xhidhix-e mat-dcdc)
           (unless (= i j)
             (mtx:add! xhidhix-e mat-dcdc-t)
             (mtx:scale! mat-dcdc-t delta)
             (mtx:add! xhidhix-g mat-dcdc-t)))))))
    (values xhidhix-g xhidhix-e)))

(define (calc-xhidhix-all eval xhi)
  "Calc_xHiDHix_all"
  (let* ((n-size (vec:length eval))
         (d-size (/ (mtx:columns xhi)
                    n-size))
         (v-size (* d-size (1+ d-size) 1/2))
         (dc-size (mtx:rows xhi))
         (xhidhix-all-g (mtx:calloc dc-size (* v-size d-size)))
         (xhidhix-all-e (mtx:calloc dc-size (* v-size d-size))))
    (dotimes (i d-size)
      (dorange
       (j i d-size)
       (let ((v (getindex i j d-size)))
         (receive (xhidhix-g xhidhix-e)
             (calc-xhidhix eval xhi i j)
           (submatrix->mtx!
            xhidhix-g xhidhix-all-g
            0 (* v dc-size) dc-size dc-size)
           (submatrix->mtx!
            xhidhix-e xhidhix-all-e
            0 (* v dc-size) dc-size dc-size)
           (mtx:free xhidhix-g xhidhix-e)))))
    (values xhidhix-all-g xhidhix-all-e)))

(define (calc-xhidhixqixhiy-all xhidhix-all-g xhidhix-all-e qixhiy)
  "Calc_xHiDHixQixHiy_all"
  (let* ((dc-size (mtx:rows xhidhix-all-g))
         (v-size (/ (mtx:columns xhidhix-all-g)
                    dc-size))
         (xhidhixqixhiy-all-g (mtx:alloc dc-size v-size))
         (xhidhixqixhiy-all-e (mtx:alloc dc-size v-size)))
    (dotimes (i v-size)
      (let ((xhidhix-g (submatrix xhidhix-all-g 0 (* i dc-size) dc-size dc-size))
            (xhidhix-e (submatrix xhidhix-all-e 0 (* i dc-size) dc-size dc-size)))
        (mtx:with-column
         (xhidhixqixhiy-g xhidhixqixhiy-all-g i)
         (mtx:with-column
          (xhidhixqixhiy-e xhidhixqixhiy-all-e i)
          (blas:gemv! xhidhix-g qixhiy xhidhixqixhiy-g)
          (blas:gemv! xhidhix-e qixhiy xhidhixqixhiy-e)
          (mtx:vec->column! xhidhixqixhiy-g xhidhixqixhiy-all-g i)
          (mtx:vec->column! xhidhixqixhiy-e xhidhixqixhiy-all-e i)))
        (mtx:free xhidhix-g xhidhix-e)))
    (values xhidhixqixhiy-all-g xhidhixqixhiy-all-e)))

(define (calc-xhidhidhiy eval hi xhi hiy i1 j1 i2 j2)
  "Calc_xHiDHiDHiy"
  (let* ((n-size (vec:length eval))
         (d-size (mtx:rows hiy))
         (dc-size (mtx:columns xhi))
         (v-size (* d-size (1+ d-size) 1/2))
         (xhidhidhiy-gg (mtx:alloc dc-size (expt v-size 2) 0))
         (xhidhidhiy-ee (mtx:copy! xhidhidhiy-gg))
         (xhidhidhiy-ge (mtx:copy! xhidhidhiy-gg)))
    (dotimes (k n-size)
      (mtx:with-column
       (xhi-col-i xhi (+ i1 (* k d-size)))
       (mtx:with-column
        (xhi-col-j xhi (+ j1 (* k d-size)))
        (let* ((delta (vec:get eval k))
               (d-hiy-i (mtx:get hiy i2 k))
               (d-hiy-j (mtx:get hiy j2 k))
               (d-hi-i1i2 (mtx:get hi i1 (+ i2 (* k d-size))))
               (d-hi-i1j2 (mtx:get hi i1 (+ j2 (* k d-size))))
               (d-hi-j1i2 (mtx:get hi j1 (+ i2 (* k d-size))))
               (d-hi-j1j2 (mtx:get hi j1 (+ j2 (* k d-size))))
               (equal? (= i2 j2)))
          (blas:axpy! (* delta delta d-hi-j1i2 d-hiy-j)
                      xhi-col-i xhidhidhiy-gg)
          (blas:axpy! (* d-hi-j1i2 d-hiy-j) xhi-col-i xhidhidhiy-ee)
          (blas:axpy! (* delta d-hi-j1i2 d-hiy-j) xhi-col-i xhidhidhiy-ge)
          (unless equal?
            (blas:axpy! (* delta delta d-hi-i1i2 d-hiy-j) xhi-col-j xhidhidhiy-gg)
            (blas:axpy! (* d-hi-i1i2 d-hiy-j) xhi-col-j xhidhidhiy-ee)
            (blas:axpy! (* delta d-hi-i1i2 d-hiy-j) xhi-col-j xhidhidhiy-ge))
          (unless (= i2 j2)
            (blas:axpy! (* delta delta d-hi-j1j2 d-hiy-i) xhi-col-i xhidhidhiy-gg)
            (blas:axpy! (* d-hi-j1j2 d-hiy-i) xhi-col-i xhidhidhiy-ee)
            (blas:axpy! (* delta d-hi-j1j2 d-hiy-j) xhi-col-i xhidhidhiy-ge)
            (unless equal?
              (blas:axpy! (* delta delta d-hi-i1j2 d-hiy-i) xhi-col-j xhidhidhiy-gg)
              (blas:axpy! (* d-hi-i1j2 d-hiy-i) xhi-col-j xhidhidhiy-ee)
              (blas:axpy! (* delta d-hi-i1j2 d-hiy-i) xhi-col-j xhidhidhiy-ge)))))))
    (values xhidhidhiy-gg xhidhidhiy-ee xhidhidhiy-ge)))

(define (calc-xhidhidhiy-all v-size eval hi xhi hiy)
  "Calc_xHiDHiDHiy_all"
  (let* ((dc-size (mtx:columns xhi))
         (d-size (mtx:rows hiy))
         (xhidhidhiy-all-gg (mtx:calloc dc-size (expt v-size 2)))
         (xhidhidhiy-all-ee (mtx:copy! xhidhidhiy-all-gg))
         (xhidhidhiy-all-ge (mtx:copy! xhidhidhiy-all-gg)))
    (dotimes (i1 d-size)
      (dorange
       (j1 i1 d-size)
       (let ((v1 (getindex i1 j1 d-size)))
         (dotimes (i2 d-size)
           (dorange
            (j2 i2 d-size)
            (let ((v2 (getindex i2 j2 d-size)))
              (receive (xhidhidhiy-gg xhidhidhiy-ee xhidhidhiy-ge)
                  (calc-xhidhidhiy eval hi xhi hiy i1 j1 i2 j2)
                (mtx:vec->column!
                 xhidhidhiy-gg xhidhidhiy-all-gg (+ v2 (* v1 v-size)))
                (mtx:vec->column!
                 xhidhidhiy-ee xhidhidhiy-all-ee (+ v2 (* v1 v-size)))
                (mtx:vec->column!
                 xhidhidhiy-ge xhidhidhiy-all-ge (+ v2 (* v1 v-size)))
                (vec:free xhidhidhiy-gg xhidhidhiy-ee xhidhidhiy-ge))))))))))

(define (calc-xhidhidhix eval hi xhi i1 j1 i2 j2)
  "Calc_xHiDHiDHix"
  (let* ((n-size (vec:length eval))
         (d-size (mtx:rows hi))
         (dc-size (mtx:rows xhi))
         (xhidhidhix-gg (mtx:alloc-square dc-size 0))
         (xhidhidhix-ee (mtx:alloc-square dc-size 0))
         (xhidhidhix-ge (mtx:alloc-square dc-size 0)))
    (mtx:with
     (mat-dcdc dc-size dc-size)
     (dotimes (k n-size)
       (mtx:with-column
        (xhi-col-i1 xhi (+ i1 (* k d-size)))
        (mtx:with-column
         (xhi-col-j1 xhi (+ j1 (* k d-size)))
         (mtx:with-column
          (xhi-col-i2 xhi (+ i2 (* k d-size)))
          (mtx:with-column
           (xhi-col-j2 xhi (+ j2 (* k d-size)))
           (let ((delta (vec:get eval k))
                 (d-hi-i1i2 (mtx:get hi i1 (+ i2 (* k d-size))))
                 (d-hi-i1j2 (mtx:get hi i1 (+ j2 (* k d-size))))
                 (d-hi-j1i2 (mtx:get hi j1 (+ i2 (* k d-size))))
                 (d-hi-j1j2 (mtx:get hi j1 (+ j2 (* k d-size))))
                 (equal? (= i1 j1)))
             (mtx:fill! mat-dcdc 0)
             (blas:ger! xhi-col-i1 xhi-col-j2 mat-dcdc #:alpha d-hi-j1i2)
             (mtx:add! xhidhidhix-ee mat-dcdc)
             (mtx:scale! mat-dcdc delta)
             (mtx:add! xhidhidhix-ge mat-dcdc)
             (mtx:scale! mat-dcdc delta)
             (mtx:add xhidhidhix-gg mat-dcdc)
             (unless equal?
               (mtx:fill! mat-dcdc 0)
               (blas:ger! xhi-col-j1 xhi-col-j2 mat-dcdc #:alpha d-hi-i1i2)
               (mtx:add! xhidhidhix-ee mat-dcdc)
               (mtx:scale! mat-dcdc delta)
               (mtx:add! xhidhidhix-ge mat-dcdc)
               (mtx:scale! mat-dcdc delta)
               (mtx:add! xhidhidhix-gg mat-dcdc))
             (unless (= i2 j2)
               (mtx:fill! mat-dcdc 0)
               (blas:ger! xhi-col-i1 xhi-col-i2 mat-dcdc #:alpha d-hi-j1j2)
               (mtx:add xhidhidhix-ee mat-dcdc)
               (mtx:scale! mat-dcdc delta)
               (mtx:add xhidhidhix-ge mat-dcdc)
               (mtx:scale! mat-dcdc delta)
               (mtx:add xhidhidhix-gg mat-dcdc)
               (unless equal?
                 (mtx:fill! mat-dcdc 0)
                 (blas:ger! xhi-col-j1 xhi-col-i2 mat-dcdc #:alpha d-hi-i1j2)
                 (mtx:add! xhidhidhix-ee mat-dcdc)
                 (mtx:scale! mat-dcdc delta)
                 (mtx:add! xhidhidhix-ge mat-dcdc)
                 (mtx:scale mat-dcdc delta)
                 (mtx:add! xhidhidhix-gg mat-dcdc))))))))))
    (values xhidhidhix-gg xhidhidhix-ee xhidhidhix-ge)))

(define (calc-xhidhidhix-all v-size eval hi xhi)
  "Calc_xHiDHiDHix_all"
  (let* ((d-size (mtx:columns xhi))
         (dc-size (mtx:rows xhi))
         (v-size (* d-size (1+ d-size) 1/2))
         (xhidhidhix-all-gg (mtx:alloc dc-size (* v-size v-size dc-size)))
         (xhidhidhix-all-ee (mtx:alloc dc-size (* v-size v-size dc-size)))
         (xhidhidhix-all-ge (mtx:alloc dc-size (* v-size v-size dc-size))))
    (dotimes (i1 d-size)
      (dorange
       (j1 i1 d-size)
       (let ((v1 (getindex i1 j1 d-size)))
         (dotimes (i2 d-size)
           (dorange
            (j2 i2 d-size)
            (let ((v2 (getindex i2 j2 d-size)))
              (unless (< v2 v1)
                (receive (xhidhidhix-gg1 xhidhidhix-ee1 xhidhidhix-ge1)
                    (calc-xhidhidhix eval hi xhi i1 j1 i2 j2)
                  (submatrix->mtx!
                   xhidhidhix-gg1 xhidhidhix-all-gg
                   0 (* d-size (+ v1 (* v2 v-size))) dc-size dc-size)
                  (submatrix->mtx!
                   xhidhidhix-ee1 xhidhidhix-all-ee
                   0 (* d-size (+ v1 (* v2 v-size))) dc-size dc-size)
                  (submatrix->mtx!
                   xhidhidhix-ge1 xhidhidhix-all-ge
                   0 (* d-size (+ v1 (* v2 v-size))) dc-size dc-size)
                  (mtx:free xhidhidhix-gg1 xhidhidhix-ee1 xhidhidhix-ge1)))))))))))

(define (calc-qivec-all qi vec-all-g vec-all-e)
  "Calc_QiVec_all"
  (let* ((qivec-all-g (vec:alloc (mtx:rows qi)))
         (qivec-all-e (vec:alloc (mtx:rows qi))))
    (dotimes (i (mtx:columns vec-all-g))
      (mtx:with-column
       (vec-g vec-all-g i)
       (mtx:with-column
        (vec-e vec-all-e i)
        (with-gsl-free
         ((quivec-g (blas:gemv qi vec-g))
          (quivec-e (blas:gemv qi vec-e)))
         (mtx:vec->column! quivec-g qivec-all-g i)
         (mtx:vec->column! quivec-e qivec-all-e i)))))))

(define (calc-qimat-all qi mat-all-g mat-all-e)
  "Calc_QiMat_all"
  (let* ((dc-size (mtx:rows qi))
         (v-size (/ (mtx:columns mat-all-g)
                    (mtx:rows mat-all-g)))
         (qimat-all-g (mtx:alloc dc-size (* v-size dc-size)))
         (qimat-all-e (mtx:copy! qimat-all-g)))
    (dotimes (i v-size)
      (with-gsl-free
       ((mat-g (submatrix mat-all-g 0 (* i dc-size) dc-size dc-size))
        (mat-e (submatrix mat-all-e 0 (* i dc-size) dc-size dc-size))
        (qimat-g (blas:gemm qi mat-g))
        (qimat-e (blas:gemm qi mat-e)))
       (submatrix->mtx! qimat-g qimat-all-g 0 (* i dc-size) dc-size dc-size)
       (submatrix->mtx! qimat-e qimat-all-e 0 (* i dc-size) dc-size dc-size)))
    (values qimat-all-g qimat-all-e)))

(define (calc-yhidhiy eval hiy i j)
  "Calc_yHiDHiy"
  (let* ((n-size (vec:length eval)))
    (let calc ((k 0)
               (yhidhiy-g 0)
               (yhidhiy-e 0))
      (if (= k n-size)
          ;; Untested, might be too smart to actually work. --aartaka.
          (values yhidhiy-g yhidhiy-e)
          (let ((delta (vec:get eval k))
                (d1 (mtx:get hiy i k))
                (d2 (mtx:get hiy j k)))
            (calc (1+ k)
                  (+ yhidhiy-g (* delta d1 d2 (if (= i j) 1 2)))
                  (+ yhidhiy-e (* d1 d2 (if (= i j) 1 2)))))))))

(define (calc-ypdpy eval hiy qixhiy xhidhiy-all-g xhidhiy-all-e
                    xhidhixqixhiy-all-g xhidhixqixhiy-all-e
                    i j)
  "Calc_yPDPy"
  (let* ((d-size (mtx:rows hiy))
         (v (getindex i j d-size)))
    (receive (ypdpy-g ypdpy-e)
        (calc-yhidhiy eval hiy i j)
      (values
       (+ ypdpy-g
          (- (* 2 (mtx:with-column (xhidhiy-g xhidhiy-all-g v)
                                   (blas:dot qixhiy xhidhiy-g))))
          (mtx:with-column (xhidhixqixhiy-g xhidhixqixhiy-all-g v)
                           (blas:dot qixhiy xhidhixqixhiy-g)))
       (+ ypdpy-e
          (- (* 2 (mtx:with-column (xhidhiy-e xhidhiy-all-e v)
                                   (blas:dot qixhiy xhidhiy-e))))
          (mtx:with-column (xhidhixqixhiy-e xhidhixqixhiy-all-e v)
                           (blas:dot qixhiy xhidhixqixhiy-e)))))))

(define (calc-tracehid eval hi i j)
  (let ((n-size (vec:length eval))
        (d-size (mtx:rows hi)))
    (let calc ((k 0)
               (thid-g 0)
               (thid-e 0))
      (if (= k n-size)
          (values thid-g thid-e)
          (let ((delta (vec:get eval k))
                (d (mtx:get hi j (+ i (* k d-size)))))
            (calc (1+ k)
                  (+ thid-g (* delta d (if (= i j) 1 2)))
                  (+ thid-e (* d (if (= i j) 1 2)))))))))

(define (calc-tracepd eval qi hi xhidhix-all-g xhidhix-all-e i j)
  (let* ((dc-size (mtx:rows qi))
         (d-size (mtx:rows hi))
         (v (getindex i j d-size)))
    (receive (tpd-g tpd-e)
        (calc-tracehid eval hi i j)
      (dotimes (k dc-size)
        (mtx:with-row
         (qi-row qi k)
         (mtx:with-column
          (xhidhix-g-col xhidhix-all-g (+ k (* v dc-size)))
          (mtx:with-column
           (xhidhix-e-col xhidhix-all-e (+ k (* v dc-size)))
           (set! tpd-g (- tpd-g (blas:ddot qi-row xhidhix-g-col)))
           (set! tpd-e (- tpd-e (blas:ddot qi-row xhidhix-e-col)))))))
      (values tpd-g tpd-e))))

(define (calc-yhidhidhiy eval hi hiy i1 j1 i2 j2)
  "Calc_yHiDHiDHiy"
  (let* ((n-size (vec:length eval))
         (d-size (mtx:rows hiy))
         (yhidhidhiy-gg 0)
         (yhidhidhiy-ee 0)
         (yhidhidhiy-ge 0))
    (dotimes (k n-size)
      (let ((delta (vec:get eval k))
            (d-hiy-i1 (mtx:get hiy i1 k))
            (d-hiy-j1 (mtx:get hiy j1 k))
            (d-hiy-i2 (mtx:get hiy i2 k))
            (d-hiy-j2 (mtx:get hiy j2 k))
            (d-hi-i1i2 (mtx:get hi i1 (+ i2 (* k d-size))))
            (d-hi-i1j2 (mtx:get hi i1 (+ j2 (* k d-size))))
            (d-hi-j1i2 (mtx:get hi j1 (+ i2 (* k d-size))))
            (d-hi-j1j2 (mtx:get hi j1 (+ j2 (* k d-size)))))
        (inc! yhidhidhiy-gg
              (* delta delta (+ (* d-hiy-i1 d-hi-j1i2 d-hiy-j2)
                                (if (= i1 j1)
                                    0
                                    (* d-hiy-j1 d-hi-i1i2 d-hiy-j2)))))
        (inc! yhidhidhiy-ee
              (+ (* d-hiy-i1 d-hi-j1j2 d-hiy-i2)
                 (if (= i1 j1)
                     0
                     (* d-hiy-j1 d-hi-i1i2 d-hiy-j2))))
        (inc! yhidhidhiy-ge
              (* delta (+ (* d-hiy-i1 d-hi-j1i2 d-hiy-j2)
                          (if (= i1 j1)
                              0
                              (* d-hiy-j1 d-hi-i1i2 d-hiy-j2)))))
        (unless (= i2 j2)
          (inc! yhidhidhiy-gg
                (* delta delta (+ (* d-hiy-i1 d-hi-j1j2 d-hiy-i2)
                                  (if (= i1 j1)
                                      0
                                      (* d-hiy-j1 d-hi-i1j2 d-hiy-i2)))))
          (inc! yhidhidhiy-ee
                (+ (* d-hiy-i1 d-hi-j1j2 d-hiy-i2)
                   (if (= i1 j1)
                       0
                       (* d-hiy-j1 d-hi-i1j2 d-hiy-i2))))
          (inc! yhidhidhiy-ee
                (* delta (+ (* d-hiy-i1 d-hi-j1j2 d-hiy-i2)
                            (if (= i1 j1)
                                0
                                (* d-hiy-j1 d-hi-i1j2 d-hiy-i2))))))))
    (values yhidhidhiy-gg yhidhidhiy-ee yhidhidhiy-ge)))

(define (calc-ypdpdpy eval hi xhi hiy qixhiy
                      xhidhiy-all-g xhidhiy-all-e
                      qixhidhiy-all-g qixhidhiy-all-e
                      xhidhixqixhiy-all-g xhidhixqixhiy-all-e
                      qixhidhixqixhiy-all-g qixhidhixqixhiy-all-e
                      xhidhidhiy-all-gg xhidhidhiy-all-ee xhidhidhiy-all-ge
                      xhidhidhix-all-gg xhidhidhix-all-ee xhidhidhix-all-ge
                      i1 j1 i2 j2)
  "Calc_yPDPDPy"
  (let* ((d-size (mtx:rows hi))
         (v-size (* d-size (1+ d-size) 1/2))
         (dc-size (mtx:rows xhi))
         (v1 (getindex i1 j1 d-size))
         (v2 (getindex i2 j2 d-size)))
    (vec:with
     (xhidhidhixqixhiy dc-size)
     (receive (ypdpdpy-gg ypdpdpy-ee ypdpdpy-ge)
         (calc-yhidhidhiy eval hi hiy i1 j1 i2 j2)
       (mtx:with-column
        (xhidhidhiy-gg1 xhidhidhiy-all-gg (+ v2 (* v1 v-size)))
        (mtx:with-column
         (xhidhidhiy-ee1 xhidhidhiy-all-ee (+ v2 (* v1 v-size)))
         (mtx:with-column
          (xhidhidhiy-ge1 xhidhidhiy-all-ge (+ v2 (* v1 v-size)))
          (mtx:with-column
           (xhidhidhiy-gg2 xhidhidhiy-all-gg (+ v2 (* v1 v-size)))
           (mtx:with-column
            (xhidhidhiy-ee2 xhidhidhiy-all-ee (+ v2 (* v1 v-size)))
            (mtx:with-column
             (xhidhidhiy-ge2 xhidhidhiy-all-ge (+ v2 (* v1 v-size)))
             (dec! ypdpdpy-gg (blas:dot qixhiy xhidhidhiy-gg1))
             (dec! ypdpdpy-ee (blas:dot qixhiy xhidhidhiy-ee1))
             (dec! ypdpdpy-ge (blas:dot qixhiy xhidhidhiy-ge1))
             (dec! ypdpdpy-gg (blas:dot qixhiy xhidhidhiy-gg2))
             (dec! ypdpdpy-ee (blas:dot qixhiy xhidhidhiy-ee2))
             (dec! ypdpdpy-ge (blas:dot qixhiy xhidhidhiy-ge2))))))))
       (mtx:with-column
        (xhidhiy-g1 xhidhiy-all-g v1)
        (mtx:with-column
         (xhidhiy-e1 xhidhiy-all-e v1)
         (mtx:with-column
          (qixhidhiy-g2 qixhidhiy-all-g v2)
          (mtx:with-column
           (qixhidhiy-e2 qixhidhiy-all-e v2)
           (dec! ypdpdpy-gg (blas:dot xhidhiy-g1 qixhidhiy-g2))
           (dec! ypdpdpy-ee (blas:dot xhidhiy-e1 qixhidhiy-e2))
           (dec! ypdpdpy-ge (blas:dot xhidhiy-g1 qixhidhiy-e2))
           (mtx:with-column
            (qixhidhiy-g1 qixhidhiy-all-g v1)
            (mtx:with-column
             (qixhidhiy-e1 qixhidhiy-all-e v1)
             (mtx:with-column
              (xhidhixqixhiy-g1 xhidhixqixhiy-all-g v1)
              (mtx:with-column
               (xhidhixqixhiy-e1 xhidhixqixhiy-all-e v1)
               (mtx:with-column
                (xhidhixqixhiy-g2 xhidhixqixhiy-all-g v2)
                (mtx:with-column
                 (xhidhixqixhiy-e2 xhidhixqixhiy-all-e v2)
                 (inc! ypdpdpy-gg (blas:dot xhidhixqixhiy-g1 qixhidhiy-g2))
                 (inc! ypdpdpy-gg (blas:dot xhidhixqixhiy-g2 qixhidhiy-g1))
                 (inc! ypdpdpy-ee (blas:dot xhidhixqixhiy-e1 qixhidhiy-e2))
                 (inc! ypdpdpy-ee (blas:dot xhidhixqixhiy-e2 qixhidhiy-e1))
                 (inc! ypdpdpy-ge (blas:dot xhidhixqixhiy-g1 qixhidhiy-e2))
                 (inc! ypdpdpy-ge (blas:dot xhidhixqixhiy-e2 qixhidhiy-g1))
                 ;; Asks for with-submatrix...
                 (let ((xhidhidhix-xx
                        (submatrix
                         xhidhidhix-all-gg
                         0 (* dc-size (+ v2 (* v1 v-size))) dc-size dc-size)))
                   (blas:dgemv! xhidhidhix-xx qixhiy xhidhidhixqixhiy #:beta 0)
                   (inc! ypdpdpy-gg (blas:dot xhidhidhixqixhiy qixhiy))
                   (mtx:free xhidhidhix-xx)
                   (set! xhidhidhix-xx
                     (submatrix
                      xhidhidhix-all-ee
                      0 (* dc-size (+ v2 (* v1 v-size))) dc-size dc-size))
                   (blas:dgemv! xhidhidhix-xx qixhiy xhidhidhixqixhiy #:beta 0)
                   (inc! ypdpdpy-ee (blas:dot xhidhidhixqixhiy qixhiy))
                   (mtx:free xhidhidhix-xx)
                   (set! xhidhidhix-xx
                     (submatrix
                      xhidhidhix-all-ge
                      0 (* dc-size (+ v2 (* v1 v-size))) dc-size dc-size))
                   (blas:dgemv! xhidhidhix-xx qixhiy xhidhidhixqixhiy #:beta 0)
                   (inc! ypdpdpy-ge (blas:dot xhidhidhixqixhiy qixhiy))
                   (mtx:free xhidhidhix-xx)
                   (mtx:with-column
                    (qixhidhixqixhiy-g1 qixhidhixqixhiy-all-g v1)
                    (mtx:with-column
                     (qixhidhixqixhiy-e1 qixhidhixqixhiy-all-e v1)
                     (dec! ypdpdpy-gg (blas:dot qixhidhixqixhiy-g1 xhidhixqixhiy-g2))
                     (dec! ypdpdpy-ee (blas:dot qixhidhixqixhiy-e1 xhidhixqixhiy-e2))
                     (dec! ypdpdpy-ge (blas:dot qixhidhixqixhiy-g1 xhidhixqixhiy-e2)))))))))))))))))))

(define (calc-tracehidhid eval hi i1 j1 i2 j2)
  "Calc_traceHiDHiD"
  (let ((n-size (vec:length eval))
        (d-size (mtx:rows hi))
        (thidhid-gg 0)
        (thidhid-ee 0)
        (thidhid-ge 0))
    (dotimes (k n-size)
      (let ((delta (vec:get eval k))
            (d-hi-i1i2 (mtx:get hi i1 (+ i2 (* k d-size))))
            (d-hi-i1j2 (mtx:get hi i1 (+ j2 (* k d-size))))
            (d-hi-j1i2 (mtx:get hi j1 (+ i2 (* k d-size))))
            (d-hi-j1j2 (mtx:get hi j1 (+ j2 (* k d-size))))
            (if/= (lambda (i j val)
                    (if (= i j)
                        0
                        val))))
        (inc! thidhid-gg
              (* delta delta (+ (* d-hi-i1j2 d-hi-j1i2)
                                (if/= i1 j1 (* d-hi-j1j2 d-hi-i1i2)))))
        (inc! thidhid-ee
              (+ (* d-hi-i1j2 d-hi-j1i2)
                 (if/= i1 j1 (* d-hi-j1j2 d-hi-i1i2))))
        (inc! thidhid-ge
              (* delta (+ (* d-hi-i1j2 d-hi-j1i2)
                          (if/= i1 j1 (* d-hi-j1j2 d-hi-i1i2)))))
        (unless (= i2 j2)
          (inc! thidhid-gg
                (* delta delta (+ (* d-hi-i1i2 d-hi-j1j2)
                                  (if/= i1 j1 (* d-hi-j1i2 d-hi-i1j2)))))
          (inc! thidhid-ee
                (+ (* d-hi-i1i2 d-hi-j1j2)
                   (if/= i1 j1 (* d-hi-j1i2 d-hi-i1j2))))
          (inc! thidhid-ge
                (* delta (+ (* d-hi-i1i2 d-hi-j1j2)
                            (if/= i1 j1 (* d-hi-j1i2 d-hi-i1j2))))))))
    (values thidhid-gg thidhid-ee thidhid-ge)))

(define (calc-tracepdpd
         eval qi hi xhi
         qixhidhix-all-g qixhidhix-all-e
         xhidhidhix-all-gg xhidhidhix-all-ee xhidhidhix-all-ge
         i1 j1 i2 j2)
  "Calc_tracePDPD"
  (let* ((dc-size (mtx:rows qi))
         (d-size (mtx:rows hi))
         (v-size (* d-size (1+ d-size) 1/2))
         (v1 (getindex i1 j1 d-size))
         (v2 (getindex i2 j2 d-size)))
    (receive (tpdpd-gg tpdpd-ee tpdpd-ge)
        (calc-tracehidhid eval hi i1 j1 i2 j2)
      (dotimes (i dc-size)
        (mtx:with-column
         (qi-row qi i)
         (mtx:with-column
          (xhidhidhix-gg-col xhidhidhix-all-gg (+ i (* dc-size (+ v2 (* v1 v-size)))))
          (mtx:with-column
           (xhidhidhix-ee-col xhidhidhix-all-ee (+ i (* dc-size (+ v2 (* v1 v-size)))))
           (mtx:with-column
            (xhidhidhix-ge-col xhidhidhix-all-ge (+ i (* dc-size (+ v2 (* v1 v-size)))))
            (dec! tpdpd-gg (blas:dot qi-row xhidhidhix-gg-col))
            (dec! tpdpd-ee (blas:dot qi-row xhidhidhix-ee-col))
            (dec! tpdpd-ge (blas:dot qi-row xhidhidhix-ge-col)))))))
      (dotimes (i dc-size)
        (mtx:with-row
         (qixhidhix-g-fullrow1 qixhidhix-all-g i)
         (mtx:with-row
          (qixhidhix-e-fullrow1 qixhidhix-all-e i)
          (let ((qixhidhix-g-row1 (subvector qixhidhix-g-fullrow1 (* v1 dc-size) dc-size))
                (qixhidhix-e-row1 (subvector qixhidhix-e-fullrow1 (* v1 dc-size) dc-size)))
            (mtx:with-column
             (qixhidhix-g-col2 qixhidhix-all-g (+ i (* v2 dc-size)))
             (mtx:with-column
              (qixhidhix-e-col2 qixhidhix-all-e (+ i (* v2 dc-size)))
              (inc! tpdpd-gg (blas:dot qixhidhix-g-row1 qixhidhix-g-col2))
              (inc! tpdpd-ee (blas:dot qixhidhix-e-row1 qixhidhix-e-col2))
              (inc! tpdpd-ge (blas:dot qixhidhix-g-row1 qixhidhix-e-col2))))))))
      (values tpdpd-gg tpdpd-ee tpdpd-ge))))

(define (calc-crt hessian-inv qi
                  qixhidhix-all-g qixhidhix-all-e
                  xhidhidhix-all-gg xhidhidhix-all-ee xhidhidhix-all-ge
                  d-size)
  "CalcCRT"
  (let* ((dc-size (mtx:rows qi))
         (c-size (/ dc-size d-size))
         (v-size (/ (mtx:rows hessian-inv) 2))
         (d-size (/ dc-size d-size))
         (b 0)
         (c 0)
         (d 0))
    (with-gsl-free
     ((m-dd (mtx:alloc-square d-size))
      (m-dcdc (mtx:alloc-square dc-size))
      (qi-si (mtx:alloc-square d-size))
      (qi-sub (mtx:alloc-square d-size))
      (qi-s (submatrix qi (* d-size (1- c-size)) (* d-size (1- c-size)) d-size d-size)))
     (mtx:copy! qi-s qi-sub)
     (receive (lu pmt signum)
         (linalg:decompose qi-sub)
       (linalg:%invert lu pmt qi-si)
       (dotimes (v1 v-size)
         (with-gsl-free
          ((qim-g1 (submatrix qixhidhix-all-g 0 (* v1 dc-size) dc-size dc-size))
           (qim-e1 (submatrix qixhidhix-all-e 0 (* v1 dc-size) dc-size dc-size))
           (qimqi-g1 (blas:gem qim-g1 qi))
           (qimqi-e1 (blas:gem qim-e1 qi))
           (qimqi-g1-s (submatrix qimqi-g1 (* d-size (1- c-size)) (* d-size (1- c-size)) d-size d-size))
           (qimqi-e1-s (submatrix qimqi-e1 (* d-size (1- c-size)) (* d-size (1- c-size)) d-size d-size))
           (qimqisqisi-g1 (blas:gem qimqi-g1-s qi-si))
           (qimqisqisi-e1 (blas:gem qimqi-e1-s qi-si)))
          (let ((trc-g1 0)
                (trc-e1 0)
                (trc-g2 0)
                (trc-e2 0)
                (trcc-gg 0)
                (trcc-ge 0)
                (trcc-ee 0)
                (trb-gg 0)
                (trb-ge 0)
                (trb-ee 0))
            (dotimes (k d-size)
              (dec! trc-g1 (mtx:get qimqisqisi-g1 k k)))
            (dotimes (k d-size)
              (dec! trc-e1 (mtx:get qimqisqisi-e1 k k)))
            (dorange
             (v2 v1 v-size)
             (with-gsl-free
              ((qim-g2 (submatrix qixhidhix-all-g 0 (* v2 dc-size) dc-size dc-size))
               (qim-e2 (submatrix qixhidhix-all-e 0 (* v2 dc-size) dc-size dc-size))
               (qimqi-g2 (blas:gem qim-g2 qi))
               (qimqi-e2 (blas:gem qim-e2 qi))
               (qimqi-g2-s (submatrix qimqi-g2 (* d-size (1- c-size)) (* d-size (1- c-size)) d-size d-size))
               (qimqi-e2-s (submatrix qimqi-e2 (* d-size (1- c-size)) (* d-size (1- c-size)) d-size d-size))
               (qimqisqisi-g2 (blas:gem qimqi-g2-s qi-si))
               (qimqisqisi-e2 (blas:gem qimqi-e2-s qi-si)))
              (dotimes (k d-size)
                (dec! trc-g2 (mtx:get qimqisqisi-g2 k k)))
              (dotimes (k d-size)
                (dec! trc-e2 (mtx:get qimqisqisi-e2 k k)))
              (blas:gem! qimqisqisi-g1 qimqisqisi-g2 m-dd #:beta 0)
              (dotimes (k d-size)
                (inc! trcc-gg (mtx:get m-dd k k)))
              (blas:gem! qimqisqisi-g1 qimqisqisi-e2 m-dd #:beta 0)
              (blas:gem! qimqisqisi-e1 qimqisqisi-g2 m-dd)
              (dotimes (k d-size)
                (inc! trcc-ge (mtx:get m-dd k k)))
              (blas:gem! qimqisqisi-e1 qimqisqisi-e2 m-dd #:beta 0)
              (dotimes (k d-size)
                (inc! trcc-ee (mtx:get m-dd k k)))
              (with-gsl-free
               ((qimqimqi-gg (blas:gem qim-g1 qimqi-g2))
                (qimqimqi-ge (blas:gem qim-g1 qimqi-e2))
                ;; Rebinding to add one more matrix to it.
                (qimqimqi-ge (begin
                               (blas:gem! qim-e1 qimqi-g2 qimqimqi-ge)
                               qimqimqi-ge))
                (qimqimqi-ee (blas:gem qim-e1 qimqi-e2))
                (qimqimqi-gg-s (submatrix qimqimqi-gg (* d-size (1- c-size)) (* d-size (1- c-size)) d-size d-size))
                (qimqimqi-ge-s (submatrix qimqimqi-ge (* d-size (1- c-size)) (* d-size (1- c-size)) d-size d-size))
                (qimqimqi-ee-s (submatrix qimqimqi-ee (* d-size (1- c-size)) (* d-size (1- c-size)) d-size d-size)))
               (blas:gem! qimqimqi-gg-s qi-si m-dd #:beta 0)
               (dotimes (k d-size)
                 (dec! trb-gg (mtx:get m-dd k k)))
               (blas:gem! qimqimqi-ge-s qi-si m-dd #:beta 0)
               (dotimes (k d-size)
                 (dec! trb-ge (mtx:get m-dd k k)))
               (blas:gem! qimqimqi-ee-s qi-si m-dd #:beta 0)
               (dotimes (k d-size)
                 (dec! trb-ee (mtx:get m-dd k k)))
               (with-gsl-free
                ((mm-gg (submatrix xhidhidhix-all-gg 0 (* dc-size (+ v2 (* v1 v-size))) dc-size dc-size))
                 (mm-ge (submatrix xhidhidhix-all-ge 0 (* dc-size (+ v2 (* v1 v-size))) dc-size dc-size))
                 (mm-ee (submatrix xhidhidhix-all-ge 0 (* dc-size (+ v2 (* v1 v-size))) dc-size dc-size))
                 (qimmqi-gg (begin
                              (blas:gem! qi mm-gg m-dcdc #:beta 0)
                              (blas:gem m-dcdc qi)))
                 (qimmqi-ge (begin
                              (blas:gem! qi mm-ge m-dcdc #:beta 0)
                              (blas:gem m-dcdc qi)))
                 (qimmqi-ee (begin
                              (blas:gem! qi mm-ee m-dcdc #:beta 0)
                              (blas:gem m-dcdc qi)))
                 (qimmqi-gg-s (submatrix qimmqi-gg (* d-size (1- c-size)) (* d-size (1- c-size)) d-size d-size))
                 (qimmqi-ge-s (submatrix qimmqi-ge (* d-size (1- c-size)) (* d-size (1- c-size)) d-size d-size))
                 (qimmqi-ee-s (submatrix qimmqi-ee (* d-size (1- c-size)) (* d-size (1- c-size)) d-size d-size)))
                (blas:gem! qimmqi-gg-s qi-si m-dd #:beta 0)
                (dotimes (k d-size)
                  (inc! trb-gg (mtx:get m-dd k k)))
                (blas:gem! qimmqi-ge-s qi-si m-dd #:beta 0)
                (dotimes (k d-size)
                  (inc! trb-ge (mtx:get m-dd k k)))
                (blas:gem! qimmqi-ee-s qi-si m-dd #:beta 0)
                (dotimes (k d-size)
                  (inc! trb-ee (mtx:get m-dd k k)))
                (let ((trd-gg (* 2 trb-gg))
                      (trd-ge (* 2 trb-ge))
                      (trd-ee (* 2 trb-ee))
                      (h-gg (* -1 (mtx:get hessian-inv v1 v2)))
                      (h-ge (* -1 (mtx:get hessian-inv v1 (+ v2 v-size))))
                      (h-ee (* -1 (mtx:get hessian-inv (+ v1 v-size) (+ v2 v-size)))))
                  (inc! b (+ (* h-gg trb-gg)
                             (* h-ge trb-ge)
                             (* h-ee trb-ee)))
                  (inc! c (+ (* h-gg (+ trcc-gg (* 1/2 trc-g1 trc-g2)))
                             (* h-ge (+ trcc-ge (* 1/2 trc-g1 trc-e2) (* 1/2 trc-e1 trc-g2)))
                             (* h-ee (+ trcc-ee (* 1/2 trc-e1 trc-e2)))))
                  (inc! d (+ (* h-gg (+ trcc-gg (* 1/2 trd-gg)))
                             (* h-ge (+ trcc-ge (* 1/2 trd-ge)))
                             (* h-ee (+ trcc-ee (* 1/2 trd-ee)))))))))))))))
    (let ((crt-a (- (* 2 d) c))
          (crt-b (* 2 b))
          (crt-c c))
      (values crt-a crt-b crt-c))))

(define (calc-dev reml? eval qi hi xhi hiy qixhiy hessian-inv gradient)
  "CalcDev"
  (let* ((dc-size (mtx:rows qi))
         (d-size (mtx:rows hi))
         (c-size (dc-size d-size))
         (v-size (* d-size (1+ d-size) 1/2))
         (hessian (mtx:alloc-square (* 2 v-size))))
    (receive (xhidhiy-all-g xhidhiy-all-e)
        (calc-xhidhiy-all eval xhi hiy)
      (receive (xhidhix-all-g xhidhix-all-e)
          (calc-xhidhix-all eval xhi)
        (receive (xhidhixqixhiy-all-g xhidhixqixhiy-all-e)
            (calc-xhidhixqixhiy-all xhidhix-all-g xhidhix-all-e qixhiy)
          (receive (xhidhidhiy-all-gg xhidhidhiy-all-ee xhidhidhiy-all-ge)
              (calc-xhidhidhiy-all v-size eval hi xhi hiy)
            (receive (xhidhidhix-all-gg xhidhidhix-all-ee xhidhidhix-all-ge)
                (calc-xhidhidhix-all v-size eval hi xhi)
              (receive (qixhidhiy-all-g qixhidhiy-all-e)
                  (calc-qivec-all qi xhidhiy-all-g xhidhiy-all-e)
                (receive (qixhidhixqixhiy-all-g qixhidhixqixhiy-all-e)
                    (calc-qivec-all qi xhidhixqixhiy-all-g xhidhixqixhiy-all-e)
                  (receive (qixhidhix-all-g qixhidhix-all-e)
                      (calc-qimat-all qi xhidhix-all-g xhidhix-all-e)
                    (dotimes (i1 d-size)
                      (dorange
                       (j1 i1 d-size)
                       (let ((v1 (getindex i1 j1 d-size)))
                         (receive (ypdpy-g ypdpy-e)
                             ;; Help me.
                             (calc-ypdpy eval hiy qixhiy xhidhiy-all-g xhidhiy-all-e
                                         xhidhixqixhiy-all-g xhidhixqixhiy-all-e
                                         i1 j1)
                           (receive (tx-g tx-e)
                               (if reml?
                                   (calc-tracepd eval qi hi xhidhix-all-g xhidhix-all-e i1 j1)
                                   (calc-tracehid eval hi i1 j1))
                             (vec:set! gradient v1 (+ (* -1/2 tx-g)
                                                      (* 1/2 ypdpy-g)))
                             (vec:set! gradient (+ v1 v-size) (+ (* -1/2 tx-e)
                                                                 (* 1/2 ypdpy-e)))))
                         (dotimes (i2 d-size)
                           (dorange
                            (j2 i2 d-size)
                            (let ((v2 (getindex i2 j2 d-size)))
                              (unless (< v2 v1)
                                (receive (ypdpdpy-gg ypdpdpy-ee ypdpdpy-ge)
                                    ;; Psychic damage taken
                                    (calc-ypdpdpy eval hi xhi hiy qixhiy
                                                  xhidhiy-all-g xhidhiy-all-e
                                                  qixhidhiy-all-g qixhidhiy-all-e
                                                  xhidhixqixhiy-all-g xhidhixqixhiy-all-e
                                                  qixhidhixqixhiy-all-g qixhidhixqixhiy-all-e
                                                  xhidhidhiy-all-gg xhidhidhiy-all-ee xhidhidhiy-all-ge
                                                  xhidhidhix-all-gg xhidhidhix-all-ee xhidhidhix-all-ge
                                                  i1 j1 i2 j2)
                                  (receive (txx-gg txx-ee txx-ge)
                                      (if reml?
                                          (calc-tracepdpd
                                           eval qi hi xhi
                                           qixhidhix-all-g qixhidhix-all-e
                                           xhidhidhix-all-gg xhidhidhix-all-ee xhidhidhix-all-ge
                                           i1 j1 i2 j2)
                                          (calc-tracehidhid eval hi i1 j1 i2 j2))
                                    (let ((dev2-gg (- (* 1/2 txx-gg) ypdpdpy-gg))
                                          (dev2-ee (- (* 1/2 txx-ee) ypdpdpy-ee))
                                          (dev2-ge (- (* 1/2 txx-ge) ypdpdpy-ge)))
                                      (mtx:set! hessian v1 v2 dev2-gg)
                                      (mtx:set! hessian (+ v1 v-size) (+ v2 v-size) dev2-ee)
                                      (mtx:set! hessian v1 (+ v2 v-size) dev2-ge)
                                      (mtx:set! hessian (+ v2 v-size) v1 dev2-ge)
                                      (unless (= v1 v2)
                                        (mtx:set! hessian v2 v1 dev2-gg)
                                        (mtx:set! hessian (+ v2 v-size) (+ v1 v-size) dev2-ee)
                                        (mtx:set! hessian v2 (+ v1 v-size) dev2-ge)
                                        (mtx:set! hessian (+ v1 v-size) v2 dev2-ge))))))))))))
                    (receive (lu pmt signum)
                        (linalg:decompose hessian)
                      (linalg:%invert hessian pmt hessian-inv)
                      (receive (crt-a crt-b crt-c)
                          (if (> c-size 1)
                              (calc-crt hessian-inv qi
                                        qixhidhix-all-g qixhidhix-all-g
                                        xhidhidhix-all-gg xhidhidhix-all-ee xhidhidhix-all-ge
                                        d-size)
                              (values 0 0 0))
                        (values hessian-inv crt-a crt-b crt-c)))))))))))))

(define (mph-nr reml? eval x y vg ve)
  (let* ((n-size (vec:length eval))
         (c-size (mtx:rows x))
         (d-size (mtx:rows y))
         (dc-size (* d-size c-size))
         (v-size (* d-size (1+ d-size) 1/2))
         (xxt (blas:syrk x))
         (hessian-inv (mtx:alloc-square (* 2 v-size))))
    (dotimes (i c-size)
      (dotimes (j i)
        (mtx:set! xxt i j (mtx:get xxt j i))))
    (vec:with
     (gradient (* 2 v-size))
     (receive (lu pmt signum)
         (linalg:decompose xxt)
       (let ((logl-const
              (if reml?
                  (+ (* -1/2 (- n-size c-size) d-size (log (* 2 +pi+)))
                     (* 1/2 d-size (linalg:%determinant-log lu)))
                  (* -1/2 n-size d-size (log (* 2 +pi+)))))
             (logl-new 0)
             (logl-old 0)
             (crt-a-save 0)
             (crt-b-save 0)
             (crt-c-save 0)
             (positive-definite? #t))
         (call-with-current-continuation
          (lambda (break)
            ;; TODO: `nr-precision' checks
            (dotimes (t (nr-iter)
                        (values hessian-inv crt-a-save crt-b-save crt-c-save logl-new))
              (unless (and (positive? t)
                           (or (< logl-new logl-old)
                               (not positive-definite?)
                               (< (- logl-new logl-old)
                                  (nr-precision))))
                (let ((vg-save (mtx:copy! vg))
                      (ve-save (mtx:copy! ve))
                      ;; KLUDGE: have to store them here for use in
                      ;; calc-dev down the line. Refactor it later
                      ;; --aartaka
                      (qi-save #f)
                      (hi-all-save #f)
                      (xhi-all-save #f)
                      (hiy-all-save #f)
                      (qixhiy-save #f))
                  (do ((scale 1.0 (/ scale 2))
                       (iter 0 (1+ iter))
                       ;; Necessary to make it run at least once.
                       (maybe-exit #f #t))
                      ((and maybe-exit
                            (and (or (not positive-definite?)
                                     (< logl-new logl-old)
                                     (> (- logl-new logl-old) 10))
                                 (< iter 10)
                                 (positive? t))))
                    (mtx:copy! vg-save vg)
                    (mtx:copy! ve-save ve)
                    (when (positive? t)
                      (update-vg-ve hessian-inv gradient scale vg ve))
                    (let ((check-definite
                           (lambda (mtx)
                             (mtx:with
                              (v-temp (mtx:rows ve) (mtx:columns ve) mtx)
                              (receive (u-temp d-temp)
                                  (eigendecomposition v-temp)
                                (dotimes (i d-size)
                                  (unless (positive? (vec:get d-temp i))
                                    (set! positive-definite? #f))))))))
                      (check-definite ve)
                      (check-definite vg))
                    (when positive-definite?
                      (receive (hi-all qi logdet-h logdet-q)
                          (calc-hi-qi eval x vg ve)
                        (set! qi-save qi)
                        (let* ((hiy-all (calc-hiy-all y hi-all))
                               (xhi-all (calc-xhi-all x hi-all))
                               (xhiy (calc-xhiy y xhi-all))
                               (qixhiy (blas:gemv qi xhiy))
                               (ypy (blas:dot qixhiy xhiy))
                               (ypy (- (calc-yhiy y hiy-all) ypy)))
                          (set! hi-all-save hi-all)
                          (set! xhi-all-save xhi-all)
                          (set! hiy-all-save hiy-all-save)
                          (set! qixhiy-save qixhiy)
                          (set! logl-new (- logl-const
                                            (* 1/2 logdet-h)
                                            (if reml?
                                                (* 1/2 logdet-q)
                                                0)
                                            (* 1/2 ypy)))))))
                  (cond
                   ((and (positive? t)
                         (or (< logl-new logl-old)
                             (not positive-definite?)))
                    (mtx:copy! vg-save vg)
                    (mtx:copy! vg-save ve)
                    (break hessian-inv crt-a-save crt-b-save crt-c-save logl-new))
                   ((< (- logl-new logl-old) (nr-precision))
                    (break hessian-inv crt-a-save crt-b-save crt-c-save logl-new)))
                  (set! logl-old logl-new)
                  (receive (hessian-inv crt-a crt-b crt-c)
                      (calc-dev
                       reml? eval qi-save hi-all-save
                       xhi-all-save hiy-all-save qixhiy-save hessian-inv gradient)
                    (mtx:scale! hessian-inv -1)
                    (set! crt-a-save crt-a)
                    (set! crt-b-save crt-b)
                    (set! crt-c-save crt-c))))))))))))

(define (mph-calc-beta eval w y vg ve)
  (let* ((n-size (vec:length eval))
         (c-size (mtx:rows w))
         (d-size (mtx:rows vg))
         (dc-size (* d-size c-size))
         (b (mtx:alloc d-size (1+ c-size)))
         (se-b (mtx:alloc d-size c-size)))
    (with-gsl-free
     ((beta (vec:alloc dc-size))
      (v-beta (mtx:alloc-square dc-size)))
     (receive (dl ultveh ultvehi)
         (eigen-proc vg ve)
       (with-gsl-free
        ((ultvehiy (blas:gem ultvehi y))
         (qi (calc-qi eval dl w))
         (whiy (vec:alloc dc-size 0))
         (qixhiy (vec:alloc dc-size))
         (qiwhiy (blas:gemv qi whiy)))
        (dotimes (i d-size)
          (dotimes (j c-size)
            (vec:set!
             whiy
             (+ i (* j d-size))
             (let calc-d ((k 0)
                          (d 0))
               (if (= k n-size)
                   d
                   (calc-d (1+ k)
                           (+ d
                              (* (mtx:get ultvehiy i k)
                                 (mtx:get w j k)
                                 (/ 1 (1+ (* (vec:get eval k)
                                             (vec:get dl i))))))))))))
        (blas:gemv! qi whiy qiwhiy #:beta 0)
        (dotimes (i c-size)
          (with-gsl-free
           ((qiwhiy-sub (subvector qiwhiy (* i d-size) d-size))
            (beta-sub (blas:gem ultveh qiwhiy-sub)))
           (dotimes (j c-size)
             (with-gsl-free
              ((qi-sub (submatrix qi (* i d-size) (* j d-size) d-size d-size))
               (v-beta-sub (if (< j i)
                               (with-gsl-free
                                ((vbeta-sym (submatrix v-beta (* j d-size) (* i d-size) d-size d-size)))
                                (mtx:transpose! vbeta-sym))
                               (with-gsl-free
                                ((qitemp-sub (blas:gemm qi-sub ultveh)))
                                (blas:gemm ultveh qitemp-sub)))))
              (submatrix->mtx! v-beta-sub v-beta (* j d-size) (* i d-size) d-size d-size)))
           (subvector->vec! beta-sub beta (* i d-size) d-size)))
        (dotimes (j (mtx:columns b))
          (dotimes (i (mtx:rows b))
            (mtx:set! b i j (vec:get beta (+ i (* j d-size))))
            (mtx:set! se-b i j
                      (sqrt (abs (mtx:get
                                  v-beta (+ i (* j d-size)) (+ i (* j d-size))))))))
        (values b se-b))))))

;; Duplicating from lmm.scm
(define gsl-cdf-chisq-q
  (foreign-library-function
   gsl:libgsl
   "gsl_cdf_chisq_Q"
   #:arg-types (list double double)
   #:return-type double))

(define (mph-calc-p eval x-row w y vg ve)
  (let* ((n-size (vec:length eval))
         (c-size (mtx:rows w))
         (d-size (mtx:rows vg))
         (dc-size (* d-size c-size)))
    (receive (dl ultveh ultvehi)
        (eigen-proc vg ve)
      (with-gsl-free
       ((ultveh ultveh)
        (ultvehi ultvehi)
        (qi (calc-qi eval dl w))
        (ultvehiy (blas:gemm ultvehi y))
        (xpx (mtx:alloc-square d-size))
        (xpy (vec:alloc d-size))
        (whix (mtx:alloc dc-size d-size))
        (whiy (vec:alloc dc-size))
        (qiwhix (mtx:alloc dc-size d-size)))
       (dotimes (i d-size)
         (receive (d1 d2)
             (let xp*calc ((k 0)
                           (d1 0)
                           (d2 0))
               (if (= k n-size)
                   (values d1 d2)
                   (let ((delta (vec:get eval k))
                         (dx (vec:get x-row k))
                         (dy (mtx:get ultvehiy i k)))
                     (xp*calc (1+ k)
                              (* dx dy (/ 1 (1+ (* delta dl))))
                              (* dx dx (/ 1 (1+ (* delta dl))))))))
           (vec:set! xpy i d1)
           (mtx:set! xpx i i d2)
           ;; FIXME: That's too stateful. But named let is not
           ;; exactly most readable eigher... --aartaka
           (dotimes (j c-size)
             (set! d1 0)
             (set! d2 0)
             (dotimes (k n-size)
               (let ((delta (vec:get eval k))
                     (dx (vec:get x-row k))
                     (dw (mtx:get w j k))
                     (dy (mtx:get ultvehiy i k)))
                 (inc! d1 (* dx dw (/ 1 (1+ (* delta dl)))))
                 (inc! d2 (* dy dw (/ 1 (1+ (* delta dl)))))))
             (mtx:set! whix (+ i (* j d-size)) i d1)
             (vec:set! whiy (+ i (* j d-size)) d2))))
       (blas:gemm! qi whix qiwhix #:beta 0)
       (blas:gemm! whix qiwhix xpx #:alpha -1 #:transpose-a #t)
       (blas:gemv! qiwhix whiy xpy #:alpha -1 #:transpose #t)
       (receive (xpx pmt signum)
           (linalg:decompose xpx)
         (linalg:solve xpx pmt xpy dl)
         (let ((vbeta (mtx:alloc-square d-size))
               (beta (blas:gemv ultveh dl #:transpose #t)))
           (linalg:%invert xpx pmt vbeta)
           ;; FIXME: Leaking permutations again.
           (blas:gemm! vbeta ultveh xpx #:beta 0)
           (blas:gemm! ultveh xpx vbeta #:beta 0)
           (let* ((d (blas:dot dl xpy))
                  (p-value (gsl-cdf-chisq-q d d-size)))
             (values p-value beta vbeta))))))))

(define gsl-cdf-chisq-qinv
  (foreign-library-function
   gsl:libgsl
   "gsl_cdf_chisq_Qinv"
   #:arg-types (list double double)
   #:return-type double))

(define (pcrt d-size p-value crt-a crt-b crt-c)
  (let* ((p-crt 0)
         (q d-size)
         (chisq (gsl-cdf-chisq-qinv p-value d-size))
         (a (/ crt-c (* 2 q (2+ q))))
         (b (1+ (/ (+ crt-a crt-b) (* 2 q))))
         (chisq-crt (/ (+ (* -1 b)
                          (sqrt (abs (+ (expt b 2)
                                        (* 4 a chisq)))))
                       (* 2 a))))
    (gsl-cdf-chisq-q chisq-crt d-size)))

(define (mvlmm-analyze markers useful-geno u eval utw uty)
  (let* ((n-size (mtx:rows uty))
         (d-size (mtx:columns uty))
         (c-size (mtx:columns utw))
         (dc-size (* d-size (1+ c-size)))
         (v-size (floor (* d-size (1+ c-size) 1/2)))
         (b-null (mtx:alloc d-size (1+ c-size)))
         (y (mtx:transpose! uty #t))
         (x (mtx:transpose! utw #t))
         (per-snp-params (make-hash-table (mtx:rows useful-geno))))
    (receive (vg ve b-sub)
        (mph-initial eval x y)
      (receive (vg ve logl-h0)
          (mph-em #:reml eval x y vg ve b-sub)
        (receive (logl-h0 hessian-inv crt-a crt-b crt-c)
            (mph-nr #:reml eval x y vg ve)
          (with-gsl-free
           ((x-sub (submatrix x 0 0 c-size n-size)))
           ;; KLUDGE: Maybe reuse the same matrix like GEMMA does?
           (mtx:free b-sub)
           (receive (b-sub se-b-null)
               (mph-calc-beta eval x-sub y vg ve)
             (receive (vg ve logl-h0)
                 (mph-em #:reml eval x y vg ve b-sub)
               (receive (logl-remle-h0 hessian-inv crt-a crt-b crt-c)
                   (mph-nr #:reml eval x y vg ve)
                 (mtx:free b-sub)
                 (receive (b-sub se-b-null)
                     (mph-calc-beta eval x-sub y vg ve)
                   (submatrix->mtx! b-sub b-null 0 0 d-size c-size)
                   ;; TODO: beta_remle_null etc.
                   (receive (vg-null ve-null logl-h0)
                       (mph-em #f eval x y vg ve b-sub)
                     (receive (logl-mle-h0 hessian-inv crt-a crt-b crt-c)
                         (mph-nr #:reml eval x y vg ve)
                       (mtx:free b-sub)
                       (receive (b-sub se-b-null)
                           (mph-calc-beta eval x-sub y vg ve)
                         (submatrix->mtx! b-sub b-null 0 0 d-size c-size)
                         (let ((utx (blas:gemm u useful-geno
                                               #:transpose-a #t #:transpose-b #t)))
                           (do ((i 0 (1+ i))
                                (markers markers (cdr markers)))
                               ((= i (mtx:rows useful-geno)))
                             (mtx:with-column
                              (x-row utx i)
                              (let* ((vg vg-null)
                                     (ve ve-null)
                                     (b b-null)
                                     (lol-h1 (mph-em #:reml eval x y vg ve b)))
                                (receive (p-wald beta vbeta)
                                    (mph-calc-p eval x-row x-sub y vg ve)
                                  (gsl-free beta vbeta)
                                  (when (< p-wald (p-nr))
                                    (receive (logl hessian-inv crt-a crt-b crt-c)
                                        (mph-nr #:reml eval x y vg ve)
                                      (receive (p-value beta vbeta)
                                          (mph-calc-p eval x-row x-sub y vg ve)
                                        (gsl-free vbeta)
                                        (set! p-wald p-value)
                                        (when (crt)
                                          (set! p-wald (pcrt d-size p-wald crt-a crt-b crt-c)))
                                        (let ((v-beta (unfold (cut = d-size <>)
                                                              (cut vec:get beta <>)
                                                              1+ 0))
                                              (v-vg '())
                                              (v-ve '())
                                              (v-vbeta '()))
                                          (dotimes (i d-size)
                                            (dotimes (j d-size)
                                              (set! v-vg (cons (mtx:get vg i j) v-vg))
                                              (set! v-ve (cons (mtx:get ve i j) v-ve))
                                              (set! v-vbeta (cons (mtx:get vbeta i j) v-vbeta))))
                                          (hash-set! per-snp-params
                                                     (car markers)
                                                     (list v-beta p-wald
                                                           0 ;;p-lrt
                                                           0 ;; p-score
                                                           (reverse! v-vg)
                                                           (reverse! v-ve)
                                                           (reverse! v-vbeta)))))))))))))))))))))))))
