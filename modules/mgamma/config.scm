(define-module (mgamma config)
  #:export (mapsize
            n-regions l-min l-max
            l-mle-null log-mle-null l-remle-null log-remle-null
            pve pve-se
            p-nr em-iter em-precision nr-iter nr-precision crt
            +pi+))

(define mapsize (make-parameter (* 100 10485760)))

;; LMM params
(define n-regions (make-parameter 10))
(define l-min (make-parameter 1e-5))
(define l-max (make-parameter 1e+5))
(define l-mle-null (make-parameter 0))
(define log-mle-null (make-parameter 0))
(define l-remle-null (make-parameter 0))
(define log-remle-null (make-parameter 0))
(define pve (make-parameter 0))
(define pve-se (make-parameter 0))

;; mvLMM params
(define p-nr (make-parameter 0.001))
(define em-iter (make-parameter 10000))
(define em-precision (make-parameter 0.0001))
(define nr-iter (make-parameter 100))
(define nr-precision (make-parameter 0.0001))
(define crt (make-parameter #f))

;; Defining π here (as the largest precision fraction from Wikipedia),
;; because I haven't found it in Guile standard lib (at least in
;; the manual).
(define +pi+ 245850922/78256779)
