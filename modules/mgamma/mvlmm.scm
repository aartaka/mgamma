(define-module (mgamma mvlmm)
  #:use-module (mgamma lmm)
  #:use-module (mgamma utils)
  #:use-module (mgamma config)
  #:use-module (srfi srfi-8)
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
  (let* ((n-size (vec:length eval))
         (c-size (mtx:rows x))
         (d-size (mtx:rows y))
         (xt (mtx:transpose! x #t))
         (y-row-tmp (vec:alloc (mtx:columns y)))
         (vg (mtx:alloc d-size d-size))
         (ve (mtx:alloc d-size d-size))
         (b (mtx:alloc d-size (1+ c-size))))
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
               (mtx:vec->column! ultvehsub b i)
               (vec:free sub ultvehsub))))
          (vec:free dl xhiy)
          (mtx:free ultveh ultvehi qi))))
    ;; TODO: Ensure everything is freed properly.
    (values vg ve b)))

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
         (ultvehie (mtx:alloc d-size n-size)))
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
            (let* ((ultvehiy (blas:gemm ultvehi y))
                   (xhiy (calc-x-hi-y eval dl x ultvehiy)))
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
      (mtx:free xxt xxti ultvehib ultvehibx ultvehiu ultvehie)
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
       (dotimes (j d-size)
         (unless (< j i)
           (let ((v (getindex i j d-size)))
             (vec:set! vec-v v (mtx:get vg i j))
             (vec:set! vec-v (+ v v-size) (mtx:get ve i j))))))
     (blas:gemv! hessian gradient vec-v #:alpha (* -1 scale))
     (dotimes (i d-size)
       (dotimes (j d-size)
         (unless (< j i)
           (let* ((v (getindex i j d-size))
                  (dg (vec:get vec-v v))
                  (de (vec:get vec-v (+ v v-size))))
             (mtx:set! vg i j dg)
             (mtx:set! vg j i dg)
             (mtx:set! ve i j de)
             (mtx:set! ve j i de))))))))

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
             (let ((hi-k (blas:gemm ultvehi mat-dd #:transpose-a #t)))
               ;; FIXME: is d-size always 2? hi-all only has 2 rows...
               (dotimes (row d-size)
                 (do ((hi-k-column 0 (1+ hi-k-column))
                      (hi-all-column (* k d-size) (1+ hi-all-column)))
                     ((= hi-all-column (* (1+ k) d-size)))
                   (mtx:set! hi-all row hi-all-column
                             (mtx:get hi-k row hi-k-column))))
               (mtx:free hi-k))))
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
      (let* ((hi-k (submatrix hi-all 0 (* k d-size) d-size d-size)))
        (mtx:with-column
         (y-k y k)
         (mtx:with-column
          (hiy-k hiy-all k)
          (blas:gemv! hi-k y-k hiy-k #:beta 0)
          (mtx:vec->column! hiy-k hiy-all k)))
        (mtx:free hi-k)))
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
      (let* ((xhi-k (submatrix xhi-all 0 (* k d-size) dc-size d-size))
             (y-k (mtx:column->vec! y k)))
        (blas:gemv! xhi-k y-k xhiy)
        (mtx:free xhi-k)
        (vec:free y-k)))
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
      (dotimes (j d-size)
        (unless (< j i)
          (let ((v (getindex i j d-size)))
            (receive (xhidhiy-g xhidhiy-e)
                (calc-xhidhiy eval xhi hiy i j)
              (mtx:vec->column! xhidhiy-g xhidhiy-all-g v)
              (mtx:vec->column! xhidhiy-e xhidhiy-all-e v)
              (vec:free xhidhiy-g xhidhiy-e))))))))

(define (calc-xhidhix eval xhi i j)
  "Calc_xHiDHix"
  (let* ((dc-size (mtx:rows xhi))
         (n-size (vec:length eval))
         (d-size (/ (mtx:columns xhi) n-size))
         (xhidhix-g (mtx:calloc dc-size dc-size))
         (xhidhix-e (mtx:calloc dc-size dc-size))
         (mat-dcdc (mtx:alloc-square dc-size))
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
            (mtx:add! xhidhix-g mat-dcdc-t))))))
    (mtx:free mat-dcdc mat-dcdc-t)
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
      (dotimes (j d-size)
        (unless (< j i)
          (let ((v (getindex i j d-size)))
            (receive (xhidhix-g xhidhix-e)
                (calc-xhidhix eval xhi i j)
              (submatrix->mtx!
               xhidhix-g xhidhix-all-g
               0 (* v dc-size) dc-size dc-size)
              (submatrix->mtx!
               xhidhix-e xhidhix-all-e
               0 (* v dc-size) dc-size dc-size)
              (mtx:free xhidhix-g xhidhix-e))))))
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

(define (calc-xhidhidhiy eval hi xhi xiy i1 j1 i2 j2)
  "Calc_xHiDHiDHiy")

(define (calc-xhidhidhiy-all v-size eval hi xhi hiy)
  "Calc_xHiDHiDHiy_all"
  (let* ((dc-size (mtx:columns xhi))
         (d-size (mtx:rows hiy))
         (xhidhidhiy-all-gg (mtx:calloc dc-size (expt v-size 2)))
         (xhidhidhiy-all-ee (mtx:copy! xhidhidhiy-all-gg))
         (xhidhidhiy-all-ge (mtx:copy! xhidhidhiy-all-gg)))
    (dotimes (i1 d-size)
      (dotimes (j1 d-size)
        (unless (< j1 i1)
          (let ((v1 (getindex i1 j1 d-size)))
            (dotimes (i2 d-size)
              (dotimes (j2 d-size)
                (unless (< j2 i2)
                  (let ((v2 (getindex i2 j2 d-size)))
                    (receive (xhidhidhiy-gg xhidhidhiy-ee xhidhidhiy-ge)
                        (calc-xhidhidhiy eval hi xhi hiy i1 j1 i2 j2)
                      (mtx:vec->column!
                       xhidhidhiy-gg xhidhidhiy-all-gg (+ v2 (* v1 v-size)))
                      (mtx:vec->column!
                       xhidhidhiy-ee xhidhidhiy-all-ee (+ v2 (* v1 v-size)))
                      (mtx:vec->column!
                       xhidhidhiy-ge xhidhidhiy-all-ge (+ v2 (* v1 v-size)))
                      (vec:free xhidhidhiy-gg xhidhidhiy-ee xhidhidhiy-ge))))))))))))

(define (calc-dev reml? eval qi hi xhi hiy qixhiy hessian-inv)
  "CalcDev"
  (let* ((dc-size (mtx:rows qi))
         (d-size (mtx:rows hi))
         (c-size (dc-size d-size))
         (v-size (* d-size (1+ d-size) 1/2)))
    (receive (xhidhiy-all-g xhidhiy-all-e)
        (calc-xhidhiy-all eval xhi hiy)
      (receive (xhidhix-all-g xhidhix-all-e)
          (calc-xhidhix-all eval xhi)
        (receive (xhidhixqixhiy-all-g xhidhixqixhiy-all-e)
            (calc-xhidhixqixhiy-all xhidhix-all-g xhidhix-all-e qixhiy)
          (receive (xhidhidhiy-all-gg xhidhidhiy-all-ee xhidhidhiy-all-ge)
              (calc-xhidhidhiy-all v-size eval hi xhi hiy)
            #t)
          )))))

(define (mph-nr reml? eval x y vg ve)
  (let* ((n-size (vec:length eval))
         (c-size (mtx:rows x))
         (d-size (mtx:rows y))
         (dc-size (* d-size c-size))
         (v-size (* d-size (1+ d-size) 1/2))
         (xxt (blas:syrk x))
         (hessian (mtx:alloc-square (* 2 v-size)))
         (gradient (vec:alloc (* 2 v-size))))
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
             (positive-definite? #t))
         (call-with-current-continuation
          (lambda (break)
            ;; TODO: `nr-precision' checks
            (dotimes (t (nr-iter))
              (unless (and (positive? t)
                           (or (< logl-new logl-old)
                               (not positive-definite?)
                               (< (- logl-new logl-old)
                                  (nr-precision))))
                (let ((vg-save (mtx:copy! vg))
                      (ve-save (mtx:copy! ve))
                      ;; KLUDGE: have to store them here for use in
                      ;; calc-deve down the line. Refactor it later
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
                    (unless (zero? t)
                      (update-vg-ve hessian gradient scale vg ve))
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
                    (break))
                   ((< (- logl-new logl-old) (nr-precision))
                    (break)))
                  (set! logl-old logl-new)
                  (calc-dev
                   reml? eval qi-save hi-all-save
                   xhi-all-save hiy-all-save qixhiy-save hessian)))))))))))

(define (mvlmm-analyze u eval utw uty)
  (let* ((n-size (mtx:rows uty))
         (d-size (mtx:columns uty))
         (c-size (mtx:columns utw))
         (dc-size (* d-size (1+ c-size)))
         (v-size (floor (* d-size (1+ c-size) 1/2)))
         (y (mtx:transpose! uty #t))
         (x (mtx:transpose! utw #t)))
    (receive (vg ve b)
        (mph-initial eval x y)
      (receive (vg ve logl-h0)
          (mph-em #t eval x y vg ve b)
        (receive (hessian crt-a crt-b crt-c logl-h0)
            (mph-nr #t eval x y vg ve)
          #t)))))
