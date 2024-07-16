(define-module (mgamma core)
  #:use-module (mgamma config)
  #:use-module (mgamma io)
  #:use-module (mgamma useful)
  #:use-module (mgamma utils)
  #:use-module (rnrs base)
  #:use-module (srfi srfi-1)
  #:use-module (srfi srfi-8)
  #:use-module (srfi srfi-26)
  #:use-module (ice-9 match)
  #:use-module (ice-9 format)
  #:use-module (system foreign)
  #:use-module (system foreign-library)
  #:use-module ((gsl core) #:prefix gsl:)
  #:use-module ((gsl matrices) #:prefix mtx:)
  #:use-module ((gsl vectors) #:prefix vec:)
  #:use-module ((gsl blas) #:prefix blas:)
  #:use-module ((gsl root) #:prefix root:)
  #:use-module ((lmdb lmdb) #:prefix mdb:)
  #:use-module (json)
  #:export (kinship-mtx
            analyze))

(gsl:set-error-handler
 (lambda* (#:optional (reason "unknown reason") (file "unknown-file") (line -1) (errno -1) #:rest rest)
   (let ((error-text
          (format #f "Error ~d (~a at ~a:~d): ~a"
                  errno (gsl:strerror errno)
                  (if (pointer? file)
                      (pointer->string file)
                      file)
                  line
                  (if (pointer? reason)
                      (pointer->string reason)
                      reason))))
     (display error-text)
     (newline)
     (error 'gsl error-text))))

(define (kinship-mtx geno-mtx markers useful-snps)
  "Calculate the kinship matrix for GENO-MTX.
Only calculate if for USEFUL-SNPS out of MARKERS."
  (let* ((n-useful-snps (hash-count (lambda (k v) #t) useful-snps))
         ;; Because we need to sort the useful SNPs into their own matrix.
         (intermediate-mtx (mtx:alloc n-useful-snps (mtx:columns geno-mtx)))
         (tmp-vec (vec:alloc (mtx:columns geno-mtx) 0)))
    (do ((i 0 (1+ i))
         (result-i 0 (if (hash-ref useful-snps (car markers))
                         (1+ result-i)
                         result-i))
         (markers markers (cdr markers)))
        ((null? markers))
      (when (hash-ref useful-snps (car markers))
        (mtx:row->vec! geno-mtx i tmp-vec)
        (mtx:vec->row! tmp-vec intermediate-mtx result-i)))
    (cleanup-mtx intermediate-mtx)
    (let ((result (blas:gemm intermediate-mtx intermediate-mtx #:transpose-a blas:+trans+)))
      (mtx:scale! result (/ 1 n-useful-snps))
      (mtx:free intermediate-mtx)
      (vec:free tmp-vec)
      result)))

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
              (let* ((uab-col (mtx:column->vec! uab index-ab))
                     (result (blas:ddot hi-eval uab-col)))
                (vec:free uab-col)
                (mtx:set! pab 0 index-ab result))
              (let* ((index-aw (abindex a p n-covariates))
                     (index-bw (abindex b p n-covariates))
                     (index-ww (abindex p p n-covariates))
                     (ps-ab (mtx:get pab (1- p) index-ab))
                     (ps-aw (mtx:get pab (1- p) index-aw))
                     (ps-bw (mtx:get pab (1- p) index-bw))
                     (ps-ww (mtx:get pab (1- p) index-ww))
                     (result (if (zero? ps-ww)
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
              (let* ((uab-col (mtx:column->vec! uab index-ab))
                     (result (blas:dot hi-hi-eval uab-col)))
                (vec:free uab-col)
                (mtx:set! ppab 0 index-ab result))
              (let* ((index-aw (abindex a p n-covariates))
                     (index-bw (abindex b p n-covariates))
                     (index-ww (abindex p p n-covariates))
                     (ps2-ab (mtx:get pab (1- p) index-ab))
                     (ps-aw (mtx:get pab (1- p) index-aw))
                     (ps-bw (mtx:get pab (1- p) index-bw))
                     (ps-ww (mtx:get pab (1- p) index-ww))
                     (ps2-aw (mtx:get ppab (1- p) index-aw))
                     (ps2-bw (mtx:get ppab (1- p) index-bw))
                     (ps2-ww (mtx:get ppab (1- p) index-ww))
                     (result (if (zero? ps-ww)
                                 ps2-ab
                                 (-
                                  (+ ps2-ab
                                     (/ (* ps-aw ps-bw ps2-ww)
                                        ps-ww
                                        ps-ww))
                                  (/ (+ (* ps-aw ps2-bw)
                                        (* ps-bw ps2-aw))
                                     ps-ww)))))
                (mtx:set! pab p index-ab result))))))))

(define (calc-pppab! uab pab ppab pppab hi-hi-hi-eval n-covariates)
  (do ((p 0 (1+ p))) ;; rows phenotypes + covariates
      ((> p (1+ n-covariates)))
    (do ((a (1+ p) (1+ a))) ;; cols in p+1..rest
        ((> a (2+ n-covariates)))
      (do ((b a (1+ b))) ;; in a..rest
          ((> b (2+ n-covariates)))
        (let ((index-ab (abindex a b n-covariates)))
          (if (zero? p)
              (let* ((uab-col (mtx:column->vec! uab index-ab))
                     (result (blas:ddot hi-hi-hi-eval uab-col)))
                (vec:free uab-col)
                (mtx:set! pab 0 index-ab result))
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
                     (result (if (zero? ps-ww)
                                 ps3-ab
                                 (+ (- ps3-ab
                                       (/ (* ps-aw ps-bw ps2-ww)
                                          (expt ps-ww 3)))
                                    (- (/ (+ (* ps-aw ps3-bw)
                                             (* ps-bw ps3-aw)
                                             (* ps2-aw ps2-bw))
                                          ps-ww))
                                    (/
                                     (+ (* ps-aw ps2-bw ps2-ww)
                                        (* ps-bw ps2-aw ps2-ww)
                                        (* ps-aw ps-bw ps3-ww))
                                     (* ps-ww ps-ww))))))
                (mtx:set! pab p index-ab result)))))))
  ;; TODO
  #f)

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
(define (calc-uab-alt! utw uty-col utx-col uab n-covariates)
  "Calculate Uab for alternative model."
  (let* ((n-inds (mtx:rows utw))
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

(define (calc-covariate-pheno y w useful-pheno-mtx cvt-mtx useful-individuals)
  "Put USEFUL-PHENO-MTX data into Y and CVT-MTX into W.
Only include the data for USEFUL-INDIVIDUALS."
  (mtx:copy! useful-pheno-mtx y)
  (do ((n-covariates (if cvt-mtx
                         (mtx:columns cvt-mtx)
                         1))
       (ci-test 0 (if (car inds)
                      (1+ ci-test)
                      ci-test))
       (i 0 (1+ i))
       (inds useful-individuals (cdr inds)))
      ((null? inds))
    (when (car inds)
      (do ((cvt 0 (1+ cvt)))
          ((= cvt n-covariates))
        (mtx:set! w ci-test cvt (mtx:get cvt-mtx i cvt))))))

;; This function exists for the sole reason of closing over UAB et
;; al. while passing the generated functions to GSL root solvers.
(define (make-log-functions calc-null? n-inds n-covariates uab eval)
  "Used to create functions closed over UAB, EVAL etc.
GEMMA uses FUNC_PARAMS struct and GSL root setters for that, but we
have closures for that in Scheme."
  (list
   ;; LogRL_dev1
   (lambda (l)
     (let ((nc-total (if calc-null?
                         n-covariates
                         (1+ n-covariates)))
           (n-index (n-index n-covariates))
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
                   (y-pkp-y (/ (- p-yy pp-yy) l)))
              (real-part
               (+ (* -1/2 trace-pk)
                  (/ (* 1/2 df y-pkp-y)
                     p-yy)))))))))))
   ;; LogRL_dev2
   (lambda (l)
     (let ((nc-total (if calc-null?
                         n-covariates
                         (1+ n-covariates)))
           (n-index (n-index n-covariates))
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
                     (trace-p trace-hi)
                     (trace-pp trace-hi-hi)
                     (_ (do ((i 0 (1+ i)))
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
                                           ps-ww)))))))
                     (trace-pkpk (/ (+ df
                                       trace-pp
                                       (- (* 2 trace-p)))
                                    (* l l)))
                     (index-ww (abindex (2+ n-covariates) (2+ n-covariates) n-covariates))
                     (p-yy (mtx:get pab nc-total index-ww))
                     (pp-yy (mtx:get ppab nc-total index-ww))
                     (ppp-yy (mtx:get pppab nc-total index-ww))
                     (ypkpy (/ (- p-yy pp-yy) l))
                     (ypkpkpy (/ (+ pp-yy
                                    ppp-yy
                                    (- (* 2 pp-yy)))
                                 (* l l))))
                (real-part
                 (- (* 1/2 trace-pkpk)
                    (* 1/2
                       df
                       (- (* 2 ypkpkpy p-yy)
                          (expt ypkpy 2))
                       (/ 1 (expt p-yy 2)))))))))))))))
   ;; LogRL_dev12
   (lambda (l)
     (let ((nc-total (if calc-null?
                         n-covariates
                         (1+ n-covariates)))
           (n-index (n-index n-covariates))
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
                     (trace-p trace-hi)
                     (trace-pp trace-hi-hi)
                     (_ (do ((i 0 (1+ i)))
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
                                           ps-ww)))))))
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
                                 (* l l)))
                     (dev1 (+ (* -1/2 trace-pk)
                              (/ (* 1/2 df ypkpy)
                                 p-yy)))
                     (dev2 (- (* 1/2 trace-pkpk)
                              (* 1/2
                                 df
                                 (- (* 2 ypkpkpy p-yy)
                                    (expt ypkpy 2))
                                 (/ 1 (expt p-yy 2))))))
                (values (real-part dev1) (real-part dev2))))))))))))
   ;; LogRL_f
   (lambda (l)
     (let ((nc-total (if calc-null?
                         n-covariates
                         (1+ n-covariates)))
           (n-index (n-index n-covariates))
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
           (vec:fill! v-temp 1.0)
           (calc-pab! uab iab v-temp n-covariates)
           (let* ((logdet-h
                   (do ((i 0 (1+ i))
                        (logdet-h 0
                                  (+ logdet-h
                                     (log (abs (vec:get v-temp i))))))
                       ((= i eval-len)
                        logdet-h)))
                  (logdet-hiw (do ((i 0 (1+ i))
                                   (logdet-hiw
                                    0
                                    (let ((index-ww (abindex (1+ i) (1+ i) n-covariates)))
                                      (+ logdet-hiw
                                         (log (mtx:get pab i index-ww))
                                         (- (log (mtx:get iab i index-ww)))))))
                                  ((= i nc-total)
                                   logdet-hiw)))
                  (index-ww (abindex (2+ n-covariates) (2+ n-covariates) n-covariates))
                  (p-yy (mtx:get pab nc-total index-ww))
                  (c (* 1/2
                        df
                        (- (log df)
                           (log (* 2 +pi+))
                           1)))
                  (index-yy (abindex (2+ n-covariates) (2+ n-covariates)
                                     n-covariates))
                  (result (- c
                             (/ logdet-h 2)
                             (/ logdet-hiw 2)
                             (* 1/2 df (log p-yy)))))
             (real-part result))))))))))

(define (calc-lambda calc-null? n-inds n-covariates uab eval)
  "Calculate lambda for null (when CALC-NULL?) or alternative model (UAB from `calc-uab-alt!')
Return (LAMBDA LOGF) values."
  (match (make-log-functions calc-null? n-inds n-covariates uab eval)
    ((log-rl-dev1 log-rl-dev2 log-rl-dev12 log-rl-f)
     (let* ((lambda-interval (/ (log (/ (l-max) (l-min))) (n-regions)))
            (sign-changes
             (remove (cut not <>)
                     (unfold (cut = (n-regions) <>)
                             (lambda (i)
                               (let* ((lambda-l (* (l-min) (exp (* lambda-interval i))))
                                      (lambda-h (* (l-min) (exp (* lambda-interval (1+ i)))))
                                      (dev1-l (log-rl-dev1 lambda-l))
                                      (dev1-h (log-rl-dev1 lambda-h)))
                                 ;; If sign flips b/w dev1-l & dev1-h
                                 (if (<= (* dev1-l dev1-h) 0)
                                     (list lambda-l lambda-h)
                                     #f)))
                             1+ 0))))
       (if (null? sign-changes)
           (let ((logf-l (log-rl-f (l-min)))
                 (logf-h (log-rl-f (l-max))))
             (if (>= logf-l logf-h)
                 (values (l-min) logf-l)
                 (values (l-max) logf-h)))
           (let rec ((sign-changes sign-changes)
                     (lam +nan.0)
                     (logf +nan.0))
             (cond
              ((null? sign-changes)
               (values lam logf))
              (else
               (match (car sign-changes)
                 ((lambda-l lambda-h)
                  ;; Yes, that's necessary because GSL throws a lot of
                  ;; errors about infinite, misshaped, etc. functions.
                  (let ((handler (gsl:set-error-handler-off!)))
                    ;; Allocating a new solver every time. Wasteful?
                    ;; GEMMA uses root:set! instead. --aartaka
                    (root:call-with
                     root:+brent-solver+
                     (lambda (solver)
                       (do ((i 0 (1+ i))
                            (approximation (root:iterate! solver)
                                           (root:iterate! solver)))
                           ;; GEMMA checks for interval too, but let's
                           ;; leave it for later.
                           ((or (= i 100)
                                (eq? approximation #f))))
                       (if (root:test-interval solver 0 1e-1)
                           (let* ((handler (gsl:set-error-handler-off!))
                                  (old-root (root:root solver))
                                  (root (or (root:optimize
                                             root:+newton-polisher+ 100 1e-5
                                             #:function log-rl-dev1
                                             #:derivative log-rl-dev2
                                             #:function+derivative log-rl-dev12
                                             #:approximate-root old-root)
                                            old-root))
                                  (_ (gsl:set-error-handler! handler)))
                             (if root
                                 (let* ((l (min (max root
                                                     (l-min))
                                                (l-max)))
                                        (logf-l (log-rl-f l)))
                                   (cond
                                    ((and (nan? lam)
                                          (nan? logf))
                                     (rec (cdr sign-changes) l logf-l))
                                    ((< logf logf-l)
                                     (rec (cdr sign-changes) l logf-l))
                                    (else
                                     (rec (cdr sign-changes) lam logf))))
                                 (rec (cdr sign-changes) lam logf)))
                           (rec (cdr sign-changes) lam logf)))
                     #:function log-rl-dev1
                     #:upper lambda-h
                     #:lower lambda-l))))))))))))

(define (calc-lambda-null utw uty-col eval)
  "Calculate lambda/logf for null model (without Uab)."
  (let* ((n-covariates (mtx:columns utw))
         (n-inds (mtx:rows utw))
         (n-index (n-index n-covariates))
         (uab (calc-uab-null utw uty-col)))
    (dynamic-wind
      (lambda ()
        #t)
      (lambda ()
        (calc-lambda #t n-inds n-covariates uab eval))
      (lambda ()
        (mtx:free uab)))))

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

(define (analyze geno markers kinship eigenvectors pheno cvt)
  "Return the per-snp params for MARKERS in GENO.
Use KINSHIP, EIGENVECTORS , PHENO, and CVT (all matrices) for
computations, but mostly clean them up into new ones and use those.
KINSHIP is computed from GENO when #f.
EIGENVECTORS are computed from KINSHIP when #f."
  (let* ((useful-individuals (useful-individuals pheno cvt))
         (useful-geno (useful-geno-mtx geno useful-individuals))
         (useful-pheno (useful-pheno-mtx pheno useful-individuals))
         (useful-snps (useful-snps useful-geno markers useful-pheno cvt))
         (cvt (or cvt
                  ;; This is not useful-kinship so that
                  ;; calc-covariate-pheno gets the right size of
                  ;; cvt-mtx.
                  (mtx:alloc (mtx:columns geno) 1 1)))
         (useful-geno (useful-geno-mtx-per-snps useful-geno markers useful-snps))
         (useful-kinship (if kinship
                             (useful-kinship-mtx kinship useful-individuals)
                             (kinship-mtx useful-geno markers useful-snps)))
         (n-covariates (mtx:columns cvt))
         (n-phenotypes (mtx:columns useful-pheno))
         (n-useful-inds (mtx:rows useful-kinship))
         (n-markers (mtx:rows useful-geno))
         (y (mtx:alloc n-useful-inds n-phenotypes 0))
         (w (mtx:alloc n-useful-inds n-covariates 0))
         (b (mtx:alloc n-phenotypes n-covariates))
         (n-index (n-index n-covariates))
         (per-snp-params (make-hash-table n-markers)))
    (cleanup-mtx useful-geno)
    (calc-covariate-pheno y w useful-pheno cvt useful-individuals)
    (center-matrix! useful-kinship)
    (receive (eval u)
        (eigendecomposition-zeroed useful-kinship)
      (let* ((u
              ;; Necessary because GSL/LAPACKE don't always provide
              ;; the right signs for eigenvectors.
              (or eigenvectors u))
             (utw (blas:gemm u w #:transpose-a blas:+transpose+))
             (uty (blas:gemm u y #:transpose-a blas:+transpose+))
             (y-col (mtx:column->vec! y 0))
             (uty-col (mtx:column->vec! uty 0))
             (uab (calc-uab-null utw uty-col))
             (utx (blas:gemm u useful-geno #:transpose-a blas:+transpose+ #:transpose-b blas:+transpose+)))
        (when (= 1 n-phenotypes)
          (receive (lam logl-h0)
              (calc-lambda-null utw uty-col eval)
            (l-mle-null lam)
            (log-mle-null logl-h0))
          ;; TODO
          #f)
        (vec:with
         (tmp (mtx:rows utx))
         (do ((i 0 (1+ i))
              (markers markers (cdr markers)))
             ((= i n-markers))
           (mtx:column->vec! utx i tmp)
           (calc-uab-alt! utw uty-col tmp uab n-covariates)
           (receive (lam logl-alt)
               (calc-lambda #f n-useful-inds n-covariates uab eval)
             (receive (beta se p-wald)
                 (calc-rlwald lam n-useful-inds n-covariates uab eval)
               (hash-set! per-snp-params (car markers)
                          (list
                           ;; Maf
                           (hash-ref useful-snps (car markers))
                           beta se
                           lam
                           p-wald
                           logl-alt))))))))
    per-snp-params))

;; (define geno (geno.txt->genotypes-mtx "/home/aartaka/git/GEMMA/example/BXD_geno.txt"))
;; (define geno-mtx (first geno))
;; (define geno-markers (second geno))
;; (define pheno-mtx (pheno.txt->pheno-mtx "/home/aartaka/git/GEMMA/example/BXD_pheno.txt"))
;; (define cvt-mtx (covariates.txt->cvt-mtx "/home/aartaka/git/GEMMA/example/mouse_hs1940_snps_anno.txt"))

;; (define kinship (kinship-mtx geno-mtx geno-markers (useful-snps geno-mtx geno-markers pheno-mtx #f)))
;; (define eigen (eigenu.txt->eigenvectors "/home/aartaka/git/GEMMA/output/BXD.eigenU.txt"))
;; (define useful-inds (useful-individuals pheno-mtx #f))
;; (define params (analyze geno-mtx geno-markers kinship eigen pheno-mtx #f))
;; (begin (hash-map->list (lambda (key value)
;;                          (format #t "~a: ~s~%" key value))
;;                        params)
;;        #t)
;; (define uab (analyze geno-mtx geno-markers kinship
;;                      pheno-mtx #f))
;; (define uty (analyze geno-mtx geno-markers kinship
;;                      pheno-mtx #f))
;; (define uty-col (analyze geno-mtx geno-markers kinship
;;                          pheno-mtx #f))
;; (vec:for-each (lambda (i e)
;;                 (newline)
;;                 (display e))
;;               $299)
;; (let ((prev-row 0))
;;   (mtx:for-each (lambda (i j e)
;;                   (when (not (= i prev-row))
;;                     (newline)
;;                     (set! prev-row i))
;;                   (format #t "~s " e))
;;                 $459))
;; (define rs3707673 (hash-ref params "rs3707673"))
;; (hash-ref params "rs6336442")
;; (hash-map->list (lambda (key value)
;;                   (format #t "~a: ~s~%"
;;                           key value))
;;                 params)
