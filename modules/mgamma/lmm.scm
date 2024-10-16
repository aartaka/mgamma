(define-module (mgamma lmm)
  #:use-module (mgamma utils)
  #:use-module (mgamma config)
  #:use-module (ice-9 match)
  #:use-module (srfi srfi-1)
  #:use-module (srfi srfi-8)
  #:use-module (srfi srfi-26)
  #:use-module (system foreign)
  #:use-module (system foreign-library)
  #:use-module ((gsl core) #:prefix gsl:)
  #:use-module ((gsl vectors) #:prefix vec:)
  #:use-module ((gsl matrices) #:prefix mtx:)
  #:use-module ((gsl blas) #:prefix blas:)
  #:use-module ((gsl root) #:prefix root:)
  #:use-module ((gsl linear-algebra) #:prefix linalg:)
  #:export (calc-lambda
            calc-lambda-null
            calc-vg-ve-beta
            lmm-analyze))

;; C++ function returns a size_t, while this one returns a signed
;; int. But then, something else is broken if we get a negative here.
(define (abindex a b n-covariates)
  (let* ((hi (max a b))
         (lo (min a b))
         (result
          ;; Infix version (doesn't make things easier):
          ;; (2 * cols - lo + 2) * (lo - 1) / 2 + hi - lo
          (+ (* (+ (* 2 (2+ n-covariates))
                   (- lo)
                   2)
                (1- lo)
                1/2)
             hi
             (- lo))))
    (when (negative? result)
      (warn (format #f "ABindex for ~a and ~a is ~a (negative), so something is broken" a b result)))
    result))

(define gsl-cdf-fdist-q
  (foreign-library-function
   gsl:libgsl
   "gsl_cdf_fdist_Q"
   #:arg-types (list double double double)
   #:return-type double))
(define gsl-cdf-chisq-q
  (foreign-library-function
   gsl:libgsl
   "gsl_cdf_chisq_Q"
   #:arg-types (list double double)
   #:return-type double))

(define (calc-pab! uab pab hi-eval n-covariates)
  (do ((p 0 (1+ p))) ;; rows phenotypes + covariates
      ((> p (1+ n-covariates)))
    (do ((a (1+ p) (1+ a))) ;; cols in p+1..rest
        ((> a (2+ n-covariates)))
      (do ((b a (1+ b))) ;; in a..rest
          ((> b (2+ n-covariates)))
        (let ((index-ab (abindex a b n-covariates)))
          (if (zero? p)
              (mtx:with-column
               (uab-col uab index-ab)
               (mtx:set! pab 0 index-ab (blas:ddot hi-eval uab-col)))
              (let* ((index-aw (abindex a p n-covariates))
                     (index-bw (abindex b p n-covariates))
                     (index-ww (abindex p p n-covariates))
                     (ps-ab (mtx:get pab (1- p) index-ab))
                     (ps-aw (mtx:get pab (1- p) index-aw))
                     (ps-bw (mtx:get pab (1- p) index-bw))
                     (ps-ww (mtx:get pab (1- p) index-ww))
                     (result (if (approximately-zero? ps-ww)
                                 ps-ab
                                 (- ps-ab (/ (* ps-aw ps-bw) ps-ww)))))
                (mtx:set! pab p index-ab result))))))))
(define (calc-ppab! uab pab ppab hi-hi-eval n-covariates)
  (do ((p 0 (1+ p))) ;; rows phenotypes + covariates
      ((> p (1+ n-covariates)))
    (do ((a (1+ p) (1+ a))) ;; cols in p+1..rest
        ((> a (2+ n-covariates)))
      (do ((b a (1+ b))) ;; in a..rest
          ((> b (2+ n-covariates)))
        (let ((index-ab (abindex a b n-covariates)))
          (if (zero? p)
              (mtx:with-column
               (uab-col uab index-ab)
               (mtx:set!
                ppab 0 index-ab
                (blas:dot hi-hi-eval uab-col)))
              (let* ((index-aw (abindex a p n-covariates))
                     (index-bw (abindex b p n-covariates))
                     (index-ww (abindex p p n-covariates))
                     (ps2-ab (mtx:get ppab (1- p) index-ab))
                     (ps-aw (mtx:get pab (1- p) index-aw))
                     (ps-bw (mtx:get pab (1- p) index-bw))
                     (ps-ww (mtx:get pab (1- p) index-ww))
                     (ps2-aw (mtx:get ppab (1- p) index-aw))
                     (ps2-bw (mtx:get ppab (1- p) index-bw))
                     (ps2-ww (mtx:get ppab (1- p) index-ww))
                     (result (if (approximately-zero? ps-ww)
                                 ps2-ab
                                 (-
                                  (+ ps2-ab
                                     (/ (* ps-aw ps-bw ps2-ww)
                                        (expt ps-ww 2)))
                                  (/ (+ (* ps-aw ps2-bw)
                                        (* ps-bw ps2-aw))
                                     ps-ww)))))
                (mtx:set! ppab p index-ab result))))))))
(define (calc-pppab! uab pab ppab pppab hi-hi-hi-eval n-covariates)
  (do ((p 0 (1+ p))) ;; rows phenotypes + covariates
      ((> p (1+ n-covariates)))
    (do ((a (1+ p) (1+ a))) ;; cols in p+1..rest
        ((> a (2+ n-covariates)))
      (do ((b a (1+ b))) ;; in a..rest
          ((> b (2+ n-covariates)))
        (let ((index-ab (abindex a b n-covariates)))
          (if (zero? p)
              (mtx:with-column
               (uab-col uab index-ab)
               (mtx:set!
                pppab 0 index-ab
                (blas:ddot hi-hi-hi-eval uab-col)))
              (let* ((index-aw (abindex a p n-covariates))
                     (index-bw (abindex b p n-covariates))
                     (index-ww (abindex p p n-covariates))
                     (ps3-ab (mtx:get pppab (1- p) index-ab))
                     (ps-aw (mtx:get pab (1- p) index-aw))
                     (ps-bw (mtx:get pab (1- p) index-bw))
                     (ps-ww (mtx:get pab (1- p) index-ww))
                     (ps2-aw (mtx:get ppab (1- p) index-aw))
                     (ps2-bw (mtx:get ppab (1- p) index-bw))
                     (ps2-ww (mtx:get ppab (1- p) index-ww))
                     (ps3-aw (mtx:get ppab (1- p) index-aw))
                     (ps3-bw (mtx:get ppab (1- p) index-bw))
                     (ps3-ww (mtx:get ppab (1- p) index-ww))
                     (result (if (approximately-zero? ps-ww)
                                 ps3-ab
                                 (+ (- ps3-ab
                                       (/ (* ps-aw ps-bw
                                             (expt ps2-ww 2))
                                          (expt ps-ww 3)))
                                    (- (/ (+ (* ps-aw ps3-bw)
                                             (* ps-bw ps3-aw)
                                             (* ps2-aw ps2-bw))
                                          ps-ww))
                                    (/
                                     (+ (* ps-aw ps2-bw ps2-ww)
                                        (* ps-bw ps2-aw ps2-ww)
                                        (* ps-aw ps-bw ps3-ww))
                                     (expt ps-ww 2))))))
                (mtx:set! pppab p index-ab result))))))))

(define (n-index n-covariates)
  (floor
   (* (+ n-covariates 3)
      (+ n-covariates 2)
      1/2)))

(define (calc-uab-null utw uty)
  "Calculate UAB for null/H0 model."
  (let* ((n-inds (mtx:rows utw))
         (n-covariates (mtx:columns utw))
         (uab (mtx:alloc n-inds (n-index n-covariates)))
         (tmp (vec:alloc n-inds)))
    (do ((a 1 (1+ a)))
        ((> a (2+ n-covariates)))
      (unless (= a (1+ n-covariates)) ;; The last column is phenotypes???
        (if (= a (2+ n-covariates))
            (vec:copy! uty tmp)
            (mtx:column->vec! utw (1- a) tmp))
        (do ((b a (1- b)))
            ((zero? b))
          (unless (= b (1+ n-covariates)) ;; Same as above
            (let* ((index-ab (abindex a b n-covariates))
                   (uab-col (mtx:column->vec! uab index-ab)))
              (if (= b (2+ n-covariates))
                  (vec:copy! uty uab-col)
                  (let ((column (mtx:column->vec! utw (1- b))))
                    (vec:copy! column uab-col)
                    (vec:free column)))
              (vec:multiply! uab-col tmp)
              (mtx:vec->column! uab-col uab index-ab))))))
    (vec:free tmp)
    uab))

;; uab is reused/modified from the result of calc-uab-null
(define (calc-uab-alt! utw uty-col utx-col uab)
  "Calculate Uab for alternative model."
  (let* ((n-inds (mtx:rows utw))
         (n-covariates (mtx:columns utw))
         (tmp (vec:alloc (mtx:rows uab))))
    (do ((b 1 (1+ b)))
        ((= b (+ n-covariates 3)))
      (let ((index-ab (abindex (1+ n-covariates) b n-covariates)))
        (cond
         ((= b (2+ n-covariates))
          (vec:copy! uty-col tmp))
         ((= b (1+ n-covariates))
          (vec:copy! utx-col tmp))
         (else
          (mtx:column->vec! utw (1- b) tmp)))
        (vec:multiply tmp utx-col)
        (mtx:vec->column! tmp uab index-ab)))
    (vec:free tmp)
    uab))

;; This function exists for the sole reason of closing over UAB et
;; al. while passing the generated functions to GSL root solvers.
(define (make-log-functions reml? calc-null? n-inds n-covariates uab eval)
  "Used to create functions closed over UAB, EVAL etc.
GEMMA uses FUNC_PARAMS struct and GSL root setters for that, but we
have closures for that in Scheme."
  (list
   ;; LogRL_dev1/LogL_dev1
   (lambda (l)
     (let ((nc-total (if calc-null?
                         n-covariates
                         (1+ n-covariates)))
           (n-index (n-index n-covariates))
           ;; REML-only
           (df (- n-inds n-covariates
                  (if calc-null?
                      0
                      1))))
       (vec:with
        (hi-eval (vec:length eval) 1)
        (vec:with
         (v-temp (vec:length eval) eval)
         (vec:scale! v-temp l)
         (vec:add-constant! v-temp 1)
         (vec:divide! hi-eval v-temp)
         (vec:with
          (hi-hi-eval (vec:length eval) hi-eval)
          (vec:multiply! hi-hi-eval hi-eval)
          (vec:fill! v-temp 1)
          (mtx:with
           (pab (2+ n-covariates) n-index)
           (calc-pab! uab pab hi-eval n-covariates)
           (mtx:with
            (ppab (2+ n-covariates) n-index)
            (calc-ppab! uab pab ppab hi-hi-eval n-covariates)
            (let* ((trace-hi (blas:dot hi-eval v-temp))
                   (trace-hik (/ (- n-inds trace-hi) l))
                   ;; REML-only
                   (trace-p
                    (do ((i 0 (1+ i))
                         (trace-p
                          trace-hi
                          (let ((index-ww (abindex (1+ i) (1+ i) n-covariates)))
                            (- trace-p
                               (/ (mtx:get ppab i index-ww)
                                  (mtx:get pab i index-ww))))))
                        ((= i nc-total)
                         trace-p)))
                   (trace-pk (/ (- df trace-p) l))
                   (index-ww (abindex (2+ n-covariates) (2+ n-covariates)
                                      n-covariates))
                   (p-yy (mtx:get pab nc-total index-ww))
                   (pp-yy (mtx:get ppab nc-total index-ww))
                   (y-pkp-y (/ (- p-yy pp-yy) l))
                   (result (+ (* -1/2 (if reml? trace-pk trace-hik))
                              (/ (* 1/2 (if reml? df n-inds) y-pkp-y)
                                 p-yy))))
              (real-part result)))))))))
   ;; LogRL_dev2/LogL_dev2
   (lambda (l)
     (let ((nc-total (if calc-null?
                         n-covariates
                         (1+ n-covariates)))
           (n-index (n-index n-covariates))
           ;; REML-only
           (df (- n-inds n-covariates
                  (if calc-null?
                      0
                      1)))
           (eval-len (vec:length eval)))
       (vec:with
        (hi-eval eval-len 1)
        (vec:with
         (v-temp eval-len eval)
         (vec:scale! v-temp l)
         (vec:add-constant! v-temp 1)
         (vec:divide! hi-eval v-temp)
         (vec:with
          (hi-hi-eval eval-len hi-eval)
          (vec:multiply! hi-hi-eval hi-eval)
          (vec:with
           (hi-hi-hi-eval eval-len hi-hi-eval)
           (vec:multiply! hi-hi-hi-eval hi-eval)
           (vec:fill! v-temp 1)
           (mtx:with
            (pab (2+ n-covariates) n-index)
            (calc-pab! uab pab hi-eval n-covariates)
            (mtx:with
             (ppab (2+ n-covariates) n-index)
             (calc-ppab! uab pab ppab hi-hi-eval n-covariates)
             (mtx:with
              (pppab (2+ n-covariates) n-index)
              (calc-pppab! uab pab ppab pppab hi-hi-hi-eval n-covariates)
              (let* ((trace-hi (blas:dot hi-eval v-temp))
                     (trace-hi-hi (blas:dot hi-hi-eval v-temp))
                     (trace-hik-hik (/
                                     (+ n-inds
                                        trace-hi-hi
                                        (- (* 2 trace-hi)))
                                     (expt l 2)))
                     ;; REML-only
                     (trace-p trace-hi)
                     (trace-pp trace-hi-hi)
                     (_ (when reml?
                          (do ((i 0 (1+ i)))
                              ((= i nc-total))
                            (let* ((index-ww (abindex (1+ i) (1+ i) n-covariates))
                                   (ps-ww (mtx:get pab i index-ww))
                                   (ps2-ww (mtx:get ppab i index-ww))
                                   (ps3-ww (mtx:get pppab i index-ww)))
                              (dec! trace-p (/ ps2-ww ps-ww))
                              (inc! trace-pp (- (/ (expt ps2-ww 2)
                                                   (expt ps-ww 2))
                                                (/ (* 2 ps3-ww)
                                                   ps-ww)))))))
                     (trace-pkpk (/ (+ df
                                       trace-pp
                                       (- (* 2 trace-p)))
                                    (expt l 2)))
                     (index-ww (abindex (2+ n-covariates) (2+ n-covariates) n-covariates))
                     (p-yy (mtx:get pab nc-total index-ww))
                     (pp-yy (mtx:get ppab nc-total index-ww))
                     (ppp-yy (mtx:get pppab nc-total index-ww))
                     (ypkpy (/ (- p-yy pp-yy) l))
                     (ypkpkpy (/ (+ p-yy
                                    ppp-yy
                                    (- (* 2 pp-yy)))
                                 (expt l 2)))
                     (result (- (* 1/2 (if reml?
                                           trace-pkpk
                                           trace-hik-hik))
                                (* 1/2
                                   (if reml? df n-inds)
                                   (- (* 2 ypkpkpy p-yy)
                                      (expt ypkpy 2))
                                   (/ 1 (expt p-yy 2))))))
                (real-part result)))))))))))
   ;; LogRL_dev12/LogL_dev12
   (lambda (l)
     (let ((nc-total (if calc-null?
                         n-covariates
                         (1+ n-covariates)))
           (n-index (n-index n-covariates))
           ;; REML-only
           (df (- n-inds n-covariates
                  (if calc-null?
                      0
                      1)))
           (eval-len (vec:length eval)))
       (vec:with
        (hi-eval eval-len 1)
        (vec:with
         (v-temp eval-len eval)
         (vec:scale! v-temp l)
         (vec:add-constant! v-temp 1)
         (vec:divide! hi-eval v-temp)
         (vec:with
          (hi-hi-eval eval-len hi-eval)
          (vec:multiply! hi-hi-eval hi-eval)
          (vec:with
           (hi-hi-hi-eval eval-len hi-hi-eval)
           (vec:multiply! hi-hi-hi-eval hi-eval)
           (vec:fill! v-temp 1)
           (mtx:with
            (pab (2+ n-covariates) n-index)
            (calc-pab! uab pab hi-eval n-covariates)
            (mtx:with
             (ppab (2+ n-covariates) n-index)
             (calc-ppab! uab pab ppab hi-hi-eval n-covariates)
             (mtx:with
              (pppab (2+ n-covariates) n-index)
              (calc-pppab! uab pab ppab pppab hi-hi-hi-eval n-covariates)
              (let* ((trace-hi (blas:dot hi-eval v-temp))
                     (trace-hi-hi (blas:dot hi-hi-eval v-temp))
                     (trace-hik (/ (- n-inds trace-hi) l))
                     (trace-hik-hik (/ (+ n-inds
                                          trace-hi-hi
                                          (- (* 2 trace-hi)))
                                       (expt l 2)))
                     ;; REML-only
                     (trace-p trace-hi)
                     (trace-pp trace-hi-hi)
                     (_ (when reml?
                          (do ((i 0 (1+ i)))
                              ((= i nc-total))
                            (let* ((index-ww (abindex (1+ i) (1+ i) n-covariates))
                                   (ps-ww (mtx:get pab i index-ww))
                                   (ps2-ww (mtx:get ppab i index-ww))
                                   (ps3-ww (mtx:get pppab i index-ww)))
                              (set! trace-p (- trace-p (/ ps2-ww ps-ww)))
                              (set! trace-pp
                                (+ trace-pp
                                   (- (/ (expt ps2-ww 2)
                                         (expt ps-ww 2))
                                      (/ (* 2 ps3-ww)
                                         ps-ww))))))))
                     ;; REML-only
                     (trace-pk (/ (- df trace-p) l))
                     (trace-pkpk (/
                                  (+ df
                                     trace-pp
                                     (- (* 2 trace-p)))
                                  (expt l 2)))
                     (index-ww (abindex (2+ n-covariates) (2+ n-covariates) n-covariates))
                     (p-yy (mtx:get pab nc-total index-ww))
                     (pp-yy (mtx:get ppab nc-total index-ww))
                     (ppp-yy (mtx:get pppab nc-total index-ww))
                     (ypkpy (/ (- p-yy pp-yy) l))
                     (ypkpkpy (/ (+ p-yy
                                    ppp-yy
                                    (- (* 2 pp-yy)))
                                 (expt l 2)))
                     (dev1 (+ (* -1/2 (if reml?
                                          trace-pk
                                          trace-hik))
                              (/ (* 1/2 (if reml? df n-inds) ypkpy)
                                 p-yy)))
                     (dev2 (- (* 1/2 (if reml?
                                         trace-pkpk
                                         trace-hik-hik))
                              (/ (* 1/2
                                    (if reml? df n-inds)
                                    (- (* 2 ypkpkpy p-yy)
                                       (expt ypkpy 2)))
                                 (expt p-yy 2)))))
                (values (real-part dev1) (real-part dev2))))))))))))
   ;; LogRL_f/LogL_f
   (lambda (l)
     (let ((nc-total (if calc-null?
                         n-covariates
                         (1+ n-covariates)))
           (n-index (n-index n-covariates))
           ;; REML-only
           (df (- n-inds n-covariates
                  (if calc-null?
                      0
                      1)))
           (eval-len (vec:length eval)))
       (vec:with
        (hi-eval eval-len 1)
        (vec:with
         (v-temp eval-len eval)
         (vec:scale! v-temp l)
         (vec:add-constant! v-temp 1)
         (vec:divide! hi-eval v-temp)
         (mtx:with
          (pab (2+ n-covariates) n-index)
          (calc-pab! uab pab hi-eval n-covariates)
          (mtx:with
           (iab (2+ n-covariates) n-index)
           (let* ((logdet-h
                   (do ((i 0 (1+ i))
                        (logdet-h 0
                                  (+ logdet-h
                                     (log (abs (vec:get v-temp i))))))
                       ((= i eval-len)
                        logdet-h)))
                  ;; This is looking weird because logdet-h should be
                  ;; calculate with old v-temp before it's modified
                  ;; here. It's also because I don't want to introduce
                  ;; yet another layer of nesting to this already
                  ;; horrible function.
                  (_ (begin
                       (vec:fill! v-temp 1.0)
                       (calc-pab! uab iab v-temp n-covariates)))
                  (logdet-hiw (do ((i 0 (1+ i))
                                   (logdet-hiw
                                    0
                                    (let ((index-ww (abindex (1+ i) (1+ i) n-covariates)))
                                      (+ logdet-hiw
                                         (log (mtx:get pab i index-ww))
                                         (- (log (mtx:get iab i index-ww)))))))
                                  ((= i nc-total)
                                   logdet-hiw)))
                  (index-yy (abindex (2+ n-covariates) (2+ n-covariates) n-covariates))
                  (p-yy (mtx:get pab nc-total index-yy))
                  (c (if reml?
                         (* 1/2
                            df
                            (- (log df)
                               (log (* 2 +pi+))
                               1))
                         (* 1/2
                            n-inds
                            (- (log n-inds)
                               (log (* 2 +pi+))
                               1))))
                  (result (- c
                             (* 1/2 logdet-h)
                             (if reml? (/ logdet-hiw 2) 0)
                             (* 1/2 (if reml? df n-inds) (log p-yy)))))
             (real-part result))))))))))

(define (calc-lambda reml? calc-null? n-inds n-covariates uab eval)
  "Calculate lambda for null (when CALC-NULL?) or alternative model.
(UAB from `calc-uab-alt!' or `calc-uab-null'.)
Return (LAMBDA LOGF) values."
  (match (make-log-functions reml? calc-null? n-inds n-covariates uab eval)
    ((log-dev1 log-dev2 log-dev12 log-f)
     (let* ((lambda-interval (/ (log (/ (l-max) (l-min))) (n-regions)))
            (sign-changes
             (remove (cut not <>)
                     (unfold (cut = (n-regions) <>)
                             (lambda (i)
                               (let* ((lambda-l (* (l-min) (exp (* lambda-interval i))))
                                      (lambda-h (* (l-min) (exp (* lambda-interval (1+ i)))))
                                      (dev1-l (log-dev1 lambda-l))
                                      (dev1-h (log-dev1 lambda-h)))
                                 ;; If sign flips b/w dev1-l & dev1-h
                                 (if (<= (* dev1-l dev1-h) 0)
                                     (list lambda-l lambda-h)
                                     #f)))
                             1+ 0)))
            ;; Yes, that's necessary because GSL throws a lot of
            ;; errors about infinite, misshaped, etc. functions.
            (handler (gsl:set-error-handler-off!)))
       (if (null? sign-changes)
           (let ((logf-l (log-f (l-min)))
                 (logf-h (log-f (l-max))))
             (if (>= logf-l logf-h)
                 (values (l-min) logf-l)
                 (values (l-max) logf-h)))
           (let rec ((sign-changes sign-changes)
                     (lam +nan.0)
                     (logf +nan.0))
             (cond
              ((null? sign-changes)
               (gsl:set-error-handler! handler)
               (let ((logf-l (log-f (l-min)))
                     (logf-h (log-f (l-max)))
                     (lam lam)
                     (logf logf))
                 (when (> logf-l logf)
                   (set! lam (l-min))
                   (set! logf logf-l))
                 (when (> logf-h logf)
                   (set! lam (l-max))
                   (set! logf logf-h))
                 (values lam logf)))
              (else
               (match (car sign-changes)
                 ((lambda-l lambda-h)
                  ;; Allocating a new solver every time. Wasteful?
                  ;; GEMMA uses root:set! instead. --aartaka
                  (root:call-with
                   root:+brent-solver+
                   (lambda (solver)
                     (do ((i 0 (1+ i))
                          (approximation (root:iterate! solver)
                                         (root:iterate! solver)))
                         ((or (= i 100)
                              (eq? approximation #f)
                              (root:test-interval solver 0 1e-1))))
                     (if (root:test-interval solver 0 1e-1)
                         (let* ((old-root (root:root solver)))
                           (root:with
                            (polisher root:+newton-polisher+
                                      #:function log-dev1
                                      #:derivative log-dev2
                                      #:function+derivative log-dev12
                                      #:approximate-root old-root)
                            (receive (approximation prev-approximation)
                                (do ((i 0 (1+ i))
                                     (approximation (root:iterate! polisher)
                                                    (root:iterate! polisher))
                                     (prev-approximation #f approximation))
                                    ((or (= i 100)
                                         (eq? approximation #f)
                                         (> approximation (l-max))
                                         (< approximation (l-min))
                                         (and prev-approximation
                                              (root:test-delta
                                               approximation prev-approximation
                                               0 1e-5)))
                                     (values approximation prev-approximation)))
                              (let ((root prev-approximation))
                                (if root
                                    (let* ((l (min (l-max)
                                                   (max root
                                                        (l-min))))
                                           (logf-l (log-f l)))
                                      (cond
                                       ((and (nan? lam)
                                             (nan? logf))
                                        (rec (cdr sign-changes) l logf-l))
                                       ((< logf logf-l)
                                        (rec (cdr sign-changes) l logf-l))
                                       (else
                                        (rec (cdr sign-changes) lam logf))))
                                    (rec (cdr sign-changes) lam logf))))))
                         (rec (cdr sign-changes) lam logf)))
                   #:function log-dev1
                   #:upper lambda-h
                   #:lower lambda-l)))))))))))

(define (calc-lambda-null reml? utw uty-col eval)
  "Calculate lambda/logf for null model (without Uab)."
  (let* ((n-covariates (mtx:columns utw))
         (n-inds (mtx:rows utw)))
    (with-gsl-free
     ((uab (calc-uab-null utw uty-col)))
     (calc-lambda reml? #t n-inds n-covariates uab eval))))

(define (calc-rlwald l n-inds n-covariates uab eval)
  "Calculate (BETA SE P-WALD) for a given UAB and Lambda."
  (let ((n-index (n-index n-covariates))
        (df (- n-inds n-covariates 1))
        (eval-len (vec:length eval)))
    (vec:with
     (hi-eval eval-len 1)
     (vec:with
      (v-temp eval-len eval)
      (vec:scale! v-temp l)
      (vec:add-constant! v-temp 1)
      (vec:divide! hi-eval v-temp)
      (mtx:with
       (pab (2+ n-covariates) n-index)
       (calc-pab! uab pab hi-eval n-covariates)
       (let* ((index-yy (abindex (2+ n-covariates) (2+ n-covariates) n-covariates))
              (index-xx (abindex (1+ n-covariates) (1+ n-covariates) n-covariates))
              (index-xy (abindex (2+ n-covariates) (1+ n-covariates) n-covariates))
              (p-yy (mtx:get pab n-covariates index-yy))
              (p-xx (mtx:get pab n-covariates index-xx))
              (p-xy (mtx:get pab n-covariates index-xy))
              (px-yy (mtx:get pab (1+ n-covariates) index-yy))
              (beta (/ p-xy p-xx))
              (tau (/ df px-yy))
              (se (abs (/ 1 (* tau p-xx))))
              (p-wald (gsl-cdf-fdist-q (* tau (- p-yy px-yy)) 1 df)))
         (values beta se p-wald)))))))

(define (calc-vg-ve-beta l eval utw uty-col)
  (let* ((n-covariates (mtx:columns utw))
         (n-inds (mtx:rows utw))
         (n-index (n-index n-covariates))
         (uab (calc-uab-null utw uty-col))
         (hiw (mtx:copy! utw)))
    (vec:with
     (v-temp (vec:length eval) eval)
     (vec:scale! v-temp l)
     (vec:add-constant! v-temp 1)
     (vec:with
      (hi-eval (vec:length eval) 1)
      (vec:divide! hi-eval v-temp)
      (vec:with
       (hiw-col-tmp (mtx:rows hiw))
       (do ((i 0 (1+ i)))
           ((= i (mtx:columns hiw)))
         (mtx:column->vec! hiw i hiw-col-tmp)
         (vec:multiply! hiw-col-tmp hi-eval)
         (mtx:vec->column! hiw-col-tmp hiw i)))
      (let ((whiw (blas:gemm hiw utw #:transpose-a #t))
            (whiy (blas:gemv hiw uty-col #:transpose #t))
            (pab (mtx:alloc (2+ n-covariates) n-index)))
        (calc-pab! uab pab hi-eval n-covariates)
        (receive (lu perms signum)
            (linalg:decompose whiw)
          (let* ((beta (linalg:solve lu perms whiy))
                 (index-yy (abindex (2+ n-covariates) (2+ n-covariates) n-covariates))
                 (p-yy (mtx:get pab n-covariates index-yy))
                 (ve (/ p-yy (- n-inds n-covariates)))
                 (vg (* ve l)))
            (mtx:free uab hiw whiw pab)
            (vec:free whiy)
            ;; TODO: se_beta, Vbeta etc.
            (values vg ve beta))))))))

(define (calc-rl-score l n-inds n-covariates eval uab)
  (let ((n-index (n-index n-covariates))
        (df (- n-inds n-covariates 1)))
    (vec:with
     (v-temp (vec:length eval) eval)
     (vec:scale! v-temp l)
     (vec:add-constant! v-temp 1)
     (vec:with
      (hi-eval (vec:length eval) 1)
      (vec:divide! hi-eval v-temp)
      (mtx:with
       (pab (2+ n-covariates) n-index)
       (calc-pab! uab pab hi-eval n-covariates)
       (let* ((index-yy (abindex (2+ n-covariates) (2+ n-covariates) n-covariates))
              (index-xx (abindex (1+ n-covariates) (1+ n-covariates) n-covariates))
              (index-xy (abindex (2+ n-covariates) (1+ n-covariates) n-covariates))
              (p-yy (mtx:get pab n-covariates index-yy))
              (p-xx (mtx:get pab n-covariates index-xx))
              (p-xy (mtx:get pab n-covariates index-xy))
              (px-yy (mtx:get pab (1+ n-covariates) index-yy))
              (beta (/ p-xy p-xx))
              (se (sqrt (abs (/ 1 (* (/ df px-yy) p-xx)))))
              (p-score (gsl-cdf-fdist-q
                        (/ (* n-inds (expt p-xy 2))
                           (* p-yy p-xx)) 1 df)))
         (values beta se p-score)))))))

(define (calc-pve eval utw uty-col l trace-g)
  (let* ((n-covariates (mtx:columns utw))
         (n-inds (mtx:rows utw))
         (n-index (n-index n-covariates)))
    (with-gsl-free
     ((uab (calc-uab-null utw uty-col)))
     (match (make-log-functions #:reml #:null n-inds n-covariates uab eval)
       ((_ logrl-dev2 . _)
        (let* ((se (sqrt (abs (/ -1 (logrl-dev2 l)))))
               (pve (/ (* trace-g l)
                       (1+ (* trace-g l))))
               (pve-se (* trace-g
                          (/ 1 (expt (1+ (* trace-g l)) 2))
                          se)))
          (values pve pve-se)))))))

(define (lmm-analyze markers useful-geno useful-inds useful-snps
                     u eval utw uty
                     n-covariates trace-g)
  (let* ((n-markers (length markers))
         (n-phenotypes (mtx:columns utw))
         (n-useful-inds (mtx:columns useful-geno))
         (uty-col (mtx:column->vec! uty 0))
         (uab (calc-uab-null utw uty-col))
         (utx (blas:gemm u useful-geno #:transpose-a #t #:transpose-b #t))
         (per-snp-params (make-hash-table n-markers)))
    (mtx:with
     (b n-phenotypes n-covariates)
     (mtx:with
      (se-b 1 n-covariates)
      (mtx:with-row
       (se-beta se-b 0)
       (receive (l-mle logl-mle)
           (calc-lambda #f #:null n-useful-inds n-covariates uab eval)
         (l-mle-null l-mle)
         (log-mle-null logl-mle)
         (receive (vg-mle-null ve-mle-null beta)
             (calc-vg-ve-beta (l-mle-null) eval utw uty-col)
           ;; TODO: Skipping (se_)beta_mle_null for now.
           (receive (l-remle log-remle)
               (calc-lambda #:reml #:null n-useful-inds n-covariates uab eval)
             (l-remle-null l-remle)
             (log-remle-null log-remle)
             ;; TODO: Skipping (se_)beta_remle_null for now.
             (receive (p p-se)
                 (calc-pve eval utw uty-col (l-remle-null) trace-g)
               ;; TODO: Print it?
               (pve p)
               (pve-se p-se))))))))
    (vec:with
     (tmp (mtx:rows utx))
     (do ((i 0 (1+ i))
          (markers markers (cdr markers)))
         ((= i n-markers))
       (mtx:column->vec! utx i tmp)
       (calc-uab-alt! utw uty-col tmp uab)
       (receive (beta se p-score)
           (calc-rl-score (l-mle-null) n-useful-inds n-covariates eval uab)
         (receive (lam-alt logl-alt)
             (calc-lambda #f #f n-useful-inds n-covariates uab eval)
           (let ((p-lrt (gsl-cdf-chisq-q (* 2 (- logl-alt (log-mle-null))) 1)))
             (hash-set! per-snp-params (car markers)
                        (list
                         ;; Maf
                         (hash-ref useful-snps (car markers))
                         beta se
                         ;; Omitting lambda-remle
                         lam-alt
                         ;; Omitting p-wald
                         p-lrt p-score
                         logl-alt)))))))
    per-snp-params))
