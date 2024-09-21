(define-module (mgamma core)
  #:use-module (mgamma config)
  #:use-module (mgamma io)
  #:use-module (mgamma useful)
  #:use-module (mgamma utils)
  #:use-module (mgamma lmm)
  #:use-module (mgamma mvlmm)
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
          (format #f "~a:~d ~a (errno ~d): ~a"
                  (if (pointer? file)
                      (pointer->string file)
                      file)
                  line
                  (gsl:strerror errno)
                  errno
                  (if (pointer? reason)
                      (pointer->string reason)
                      reason))))
     (display error-text)
     (newline)
     (apply error 'gsl error-text rest))))

(define (kinship-mtx geno-mtx markers useful-snps)
  "Calculate the kinship matrix for GENO-MTX.
Only calculate if for USEFUL-SNPS out of MARKERS."
  (let* ((n-useful-snps (hash-count (lambda (k v) #t) useful-snps))
         (tmp-vec (vec:alloc (mtx:columns geno-mtx) 0)))
    (with-gsl-free
     ;; Because we need to sort the useful SNPs into their own matrix.
     ((intermediate-mtx (mtx:alloc n-useful-snps (mtx:columns geno-mtx)))
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
       result))))

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

(define (analyze geno markers kinship eigenvectors pheno pheno-nums cvt)
  "Return the per-snp params for MARKERS in GENO.
Use KINSHIP, EIGENVECTORS , PHENO, and CVT (all matrices) for
computations, but mostly clean them up into new ones and use those.
KINSHIP is computed from GENO when #f.
EIGENVECTORS are computed from KINSHIP when #f.
In case PHENO-NUMS is a non-empty non-false list of zero-based
numbers, run multivariate LMM on the data instead of univariate."
  (let* ((useful-individuals (useful-individuals pheno cvt))
         (useful-geno (useful-geno-mtx geno useful-individuals))
         (useful-pheno (useful-pheno-mtx pheno useful-individuals pheno-nums))
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
         (n-phenotypes (if pheno-nums
                           (length pheno-nums)
                           ;; (mtx:columns useful-pheno)
                           1))
         (n-useful-inds (mtx:rows useful-kinship))
         (n-markers (mtx:rows useful-geno))
         (y (mtx:alloc n-useful-inds n-phenotypes 0))
         (w (mtx:alloc n-useful-inds n-covariates 0)))
    (cleanup-mtx useful-geno)
    (calc-covariate-pheno y w useful-pheno cvt useful-individuals)
    (center-matrix! useful-kinship)
    (receive (eval u)
        (eigendecomposition-zeroed useful-kinship)
      (let* ((u
              ;; Necessary because GSL/LAPACKE don't always provide
              ;; the right signs for eigenvectors.
              (or eigenvectors u))
             (utw (blas:gemm u w #:transpose-a #t))
             (uty (blas:gemm u y #:transpose-a #t))
             (y-col (mtx:column->vec! y 0))
             (uty-col (mtx:column->vec! uty 0)))
        (if (= 1 n-phenotypes)
            (lmm-analyze markers useful-geno useful-individuals useful-snps u eval utw uty-col n-covariates)
            (mvlmm-analyze markers useful-geno useful-snps u eval utw uty))))))

;; (define geno (geno.txt->genotypes-mtx "/home/aartaka/Downloads/iron/iron_geno.txt"))
;; (define geno-mtx (first geno))
;; (define geno-markers (second geno))
;; (define pheno-mtx (pheno.txt->pheno-mtx "/home/aartaka/Downloads/iron/iron_pheno.txt"))
;; (define cvt-mtx (covariates.txt->cvt-mtx "/home/aartaka/git/GEMMA/example/mouse_hs1940_snps_anno.txt"))

;; (define kinship (kinship-mtx geno-mtx geno-markers (useful-snps geno-mtx geno-markers pheno-mtx #f)))
;; (define eigen (eigenu.txt->eigenvectors "/home/aartaka/git/GEMMA/output/iron.eigenU.txt"))
;; (define useful-inds (useful-individuals pheno-mtx #f))
;; (define params (analyze geno-mtx geno-markers kinship #f pheno-mtx '(0 1) #f))
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
