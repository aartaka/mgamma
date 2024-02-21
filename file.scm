(define-module (gemma core)
  #:use-module (ice-9 ports)
  #:use-module (ice-9 textual-ports)
  #:use-module (srfi srfi-1)
  #:use-module (srfi srfi-43)
  #:use-module (gemma matrix)
  #:use-module (gemma vector))

(define (read-separated-lines file)
  (call-with-port (open-input-file file)
    (lambda (port)
      (let read-lines ((line (get-line port)))
        (if (eof-object? line)
            '()
            (cons (remove (lambda (s)
                            (zero? (string-length s)))
                          (string-split
                           line (lambda (c) (memv c '(#\Tab #\Space #\,)))))
                  (read-lines (get-line port))))))))

;; (first (read-separated-lines "/home/aartaka/git/GEMMA/example/mouse_hs1940.geno.txt"))

(define (read-bim file)
  (let ((lines (read-separated-lines file)))
    (map (lambda (split)
           (let (;; How are morgans/centimorgans marked?
                 (position (string->number (third split)))
                 (base-pair-coordinate (string->number (fourth split))))
             `(,(first split)
               ,(second split)
               ,position
               ,base-pair-coordinate
               ,@(cddddr split))))
         lines)))

;; (read-bim "/home/aartaka/git/GEMMA/example/mouse_hs1940.bim")

(define (string->num/nan str)
  (if (string= "NA" str)
      +nan.0
      (string->number str)))

(define (read-anno.txt file)
  (let ((lines (read-separated-lines file)))
    (map (lambda (split)
           (let ((marker (first split))
                 (pos (second split))
                 (chromosome (third split))
                 (unidentified? (fourth split)))
             (list marker
                   (string->num/nan pos)
                   (string->number chromosome)
                   (string->number unidentified?))))
         lines)))

;; (read-anno.txt "/home/aartaka/git/GEMMA/example/mouse_hs1940.anno.txt")

(define (read-geno.txt file)
  (let ((lines (read-separated-lines file)))
    (map (lambda (line)
           (let* ((numbers (list->vector (map string->num/nan (cdddr line))))
                  (mean (vector-mean numbers))
                  (un-nan (vector-replace-nan numbers mean)))
             (list (first line)
                   (second line)
                   (third line)
                   un-nan)))
         lines)))

;; (first (read-geno.txt "/home/aartaka/git/GEMMA/example/BXD_geno.txt"))

(define (kinship genotypes)
  (let ((cleaned-up (vector-map
                     (lambda (i gs)
                       (let* ((mean (vector-mean gs))
                              (un-nan (vector-replace-nan gs mean))
                              (centered (vector-map (lambda (i elem)
                                                      (- elem mean))
                                                    un-nan))
                              (var (vector-variance centered)))
                         (vector-map (lambda (i elem)
                                       (/ elem
                                          (if (= var 0)
                                              1
                                              (sqrt var))))
                                     centered)))
                     genotypes)))
    (matrix-multiply (matrix-transpose cleaned-up)
                     cleaned-up)))

;; (vector-ref (kinship (list->vector (map fourth (read-geno.txt "/home/aartaka/git/GEMMA/example/mouse_hs1940.geno.txt")))) 0)



