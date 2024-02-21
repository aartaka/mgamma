(define-module (gemma stat)
  #:use-module
  #:use-module (srfi srfi-1))

(define (variance values)
  (let ((vector-mean values))
    (fold + 0 (map ))))
