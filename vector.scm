(define-module (gemma vector)
  #:use-module (srfi srfi-1)
  #:use-module (srfi srfi-43))

(define-public (vector-filter pred? v)
  (let* ((new-size (vector-count (lambda (i x) (not (pred? x))) v))
         (new-vec (make-vector new-size)))
    (let filter-loop ((old-idx 0)
                      (new-idx 0))
      (unless (= old-idx (vector-length v))
        (unless (pred? (vector-ref v old-idx))
          (vector-set! new-vec new-idx (vector-ref v old-idx)))
        (filter-loop
         (+ 1 old-idx)
         (if (pred? (vector-ref v old-idx))
             new-idx
             (+ 1 new-idx)))))
    new-vec))

;; (vector-filter nan? #(1 2 3 4 5 +nan.0))

(define-public (vector-mean vals)
  (/ (vector-fold (lambda (i acc elem)
                    (+ elem acc))
                  0
                  (vector-filter nan? vals))
     ;; Length
     (vector-count (lambda (i x) (not (nan? x))) vals)))

(define-public (vector-variance vals)
  (let ((mean (vector-mean vals)))
    (/ (fold + 0 (map (lambda (x)
                        (expt (- x mean) 2))
                      (vector->list vals)))
       (- (vector-length vals) 1))))

;; (vector-variance #(8 22 61))

(define-public (vector-replace-nan vec replacement)
  (vector-map (lambda (idx elem)
                (if (nan? elem)
                    replacement
                    elem))
              vec))

;; (vector-replace-nan #(1 2 3 4 5 +nan.0) 888)
