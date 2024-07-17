(define-module (mgamma config)
  #:export (mapsize
            n-regions l-min l-max
            l-mle-null log-mle-null
            +pi+))

(define mapsize (make-parameter (* 100 10485760)))

;; LMM params
(define n-regions (make-parameter 10))
(define l-min (make-parameter 1e-5))
(define l-max (make-parameter 1e+5))
(define l-mle-null (make-parameter 0))
(define log-mle-null (make-parameter 0))

;; Defining π here (as the largest precision fraction from Wikipedia),
;; because I haven't found it in Guile standard lib (at least in
;; the manual).
(define +pi+ 245850922/78256779)
