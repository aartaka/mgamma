(define-module (gemma matrix)
  #:use-module (srfi srfi-1))

(define (tree->matrix tree)
  (list->vector (map list->vector tree)))

(define-public (matrix-transpose matrix)
  (let ((column-length (vector-length (vector-ref matrix 0))))
    (tree->matrix (let conversion ((idx 0))
                    (if (= idx column-length)
                        '()
                        (cons (map (lambda (row)
                                     (vector-ref row idx))
                                   (vector->list matrix))
                              (conversion (+ 1 idx))))))))

;; (matrix-transpose #(#(6 4 24)
;;                     #(1 −9 8)))

(define (matrix-add m1 m2)
  (unless (and (= (vector-length m1)
                  (vector-length m2))
               (= (vector-length (vector-ref m1 0))
                  (vector-length (vector-ref m2 0))))
    (error "Matrices are not of the same dimensions!"))
  (let ((add-rows (lambda (r1 r2)
                    (map + (vector->list r1) (vector->list r2)))))
    (tree->matrix (map (lambda (r1 r2)
                         (add-rows r1 r2))
                       (vector->list m1)
                       (vector->list m2)))))

;; (matrix-add
;;  #(#(3 8)
;;    #(4 6))
;;  #(#(4 0)
;;    #(1 -9)))

(define (matrix-scalar-multiply scalar m1)
  (tree->matrix
   (map (lambda (row)
          (map (lambda (e) (* e scalar))
               (vector->list row)))
        (vector->list m1))))

;; (matrix-scalar-multiply 2 #(#(4 0) #(1 -9)))

(define (matrix-negate m1)
  (matrix-scalar-multiply -1 m1))

;; (matrix-negate
;;  #(#(2 -4)
;;    #(7 10)))

(define (matrix-subtract m1 m2)
  (matrix-add m1 (matrix-negate m2)))

;; (matrix-subtract
;;  #(#(3 8)
;;    #(4 6))
;;  #(#(4 0)
;;    #(1 -9)))

(define (matrix-row m row)
  (list-ref (vector->list m) row))

(define (matrix-column m column)
  (list->vector
   (map (lambda (v)
          (vector-ref v column))
        (vector->list m))))

(define (matrix-x m)
  (vector-length (matrix-row m 0)))

(define (matrix-y m)
  (vector-length (matrix-column m 0)))

(define (matrix-ref m r c)
  (vector-ref
   (vector-ref m c)
   r))

(define (dot-product v1 v2)
  (fold
   + 0
   (map (lambda (e1 e2)
          (* e1 e2))
        (vector->list v1)
        (vector->list v2))))

(define-public (matrix-multiply m1 m2)
  (let* ((result-x (matrix-x m2))
         (result-y (matrix-y m1)))
    (tree->matrix
     (let rec-y ((y-idx 0))
       (if (= y-idx result-y)
           '()
           (cons (let rec-x ((x-idx 0))
                   (if (= x-idx result-x)
                       '()
                       (cons (dot-product (matrix-row m1 x-idx)
                                          (matrix-column m2 y-idx))
                             (rec-x (+ 1 x-idx)))))
                 (rec-y (+ 1 y-idx))))))))

;; (matrix-multiply
;;  #(#(1 2 3)
;;    #(4 5 6))
;;  #(#(7 8)
;;    #(9 10)
;;    #(11 12)))

(define (matrix2-determinant m)
  (- (* (matrix-ref m 0 0)
        (matrix-ref m 1 1))
     (* (matrix-ref m 0 1)
        (matrix-ref m 1 0))))

;; (matrix2-determinant
;;  #(#(3 8)
;;    #(4 6)))

(define (matrix-except m row column)
  (tree->matrix
   (let rec-y ((y-idx 0))
     (cond
      ((= y-idx (matrix-y m))
       '())
      ((= y-idx row)
       (rec-y (+ 1 y-idx)))
      (else
       (cons (let rec-x ((x-idx 0))
               (cond
                ((= x-idx (matrix-x m))
                 '())
                ((= x-idx column)
                 (rec-x (+ 1 x-idx)))
                (else
                 (cons (matrix-ref m x-idx y-idx)
                       (rec-x (+ 1 x-idx))))))
             (rec-y (+ 1 y-idx))))))))

;; (matrix-except
;;  #(#(6 1 1)
;;    #(4 -2 5)
;;    #(2 8 7))
;;  0 1)

(define (matrix-determinant m)
  (if (= 2 (matrix-x m) (matrix-y m))
      (matrix2-determinant m)
      (let rec ((x-idx 0))
        (if (= x-idx (matrix-x m))
            0
            (+
             ((if (even? x-idx) + -)
              (* (matrix-ref m 0 x-idx)
                 (matrix-determinant (matrix-except m x-idx 0))))
             (rec (+ 1 x-idx)))))))

;; (matrix-determinant
;;  #(#(6 1 1)
;;    #(4 -2 5)
;;    #(2 8 7)))

;; (matrix-determinant
;;  #(#(6 1 1 1)
;;    #(4 -2 5 1)
;;    #(2 8 7 1)
;;    #(1 8 4 5 )))

(define (identity-matrix size)
  (tree->matrix
   (let rec-y ((y 0))
     (if (= y size)
         '()
         (cons (let rec-x ((x 0))
                 (if (= x size)
                     '()
                     (cons (if (= y x)
                               1
                               0)
                           (rec-x (+ 1 x)))))
               (rec-y (+ 1 y)))))))

;; (identity-matrix 3)

(define (matrix-minors m)
  (tree->matrix
   (let rec-y ((y 0))
     (if (= y (matrix-y m))
         '()
         (cons (let rec-x ((x 0))
                 (if (= x (matrix-x m))
                     '()
                     (cons
                      (matrix-determinant (matrix-except m y x))
                      (rec-x (+ 1 x)))))
               (rec-y (+ 1 y)))))))

;; (matrix-minors #(#(3 0 2)
;;                  #(2 0 -2)
;;                  #(0 1 1)))

(define (matrix-cofactors m)
  (tree->matrix
   (let rec-y ((y 0))
     (if (= y (matrix-y m))
         '()
         (cons (let rec-x ((x 0))
                 (if (= x (matrix-x m))
                     '()
                     (cons
                      (if (or (and (even? x) (odd? y))
                              (and (even? y) (odd? x)))
                          (- (matrix-ref m x y))
                          (matrix-ref m x y))
                      (rec-x (+ 1 x)))))
               (rec-y (+ 1 y)))))))

(define (matrix-inverse m)
  (matrix-scalar-multiply
   (/ 1 (matrix-determinant m))
   (matrix-transpose
    (matrix-cofactors
     (matrix-minors m)))))

;; (matrix-inverse #(#(3 0 2)
;;                   #(2 0 -2)
;;                   #(0 1 1)))

(define (center-matrix mat)
  ;; w <- rep(1, nrow(mat))
  ;; #  gsl_blas_dgemv(CblasNoTrans, 1.0, G, w, 0.0, Gw);
  ;; mat %*% w -> matw
  ;; #  gsl_blas_dsyr2(CblasUpper, -1.0 / (double)G->size1, Gw, w, G);
  ;; foo <- mat - (matw %*% t(w) + w %*% t(matw)) / nrow(mat)
  ;; #  gsl_blas_ddot(w, Gw, &d);
  ;; d <- t(w) %*% matw
  ;; A = α * x * t(x) + A
  ;; #  gsl_blas_dsyr(CblasUpper, d / ((double)G->size1 * (double)G->size1), w, #G);
  ;; out <- foo + (w %*% t(w)) * as.numeric(d / (nrow(mat)^2))
  ;; return(out)
  (let* ((w (vector (make-vector (matrix-y mat) 1)))
         (matw (matrix-multiply mat w))
         (foo
          (matrix-subtract
           mat
           (matrix-scalar-multiply
            (matrix-add
             (matrix-multiply matw (matrix-transpose w))
             (matrix-multiply w (matrix-transpose matw)))
            (/ 1 (matrix-y mat)))))
         (d (dot-product (matrix-transpose w) matw))
         ;; out <- foo + (w %*% t(w)) * as.numeric(d / (nrow(mat)^2))
         (out (matrix-add
               foo
               (matrix-scalar-multiply
                (matrix-multiply w (matrix-transpose w))
                (/ d (expt (matrix-y mat) 2))))))
    out))
