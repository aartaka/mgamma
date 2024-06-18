;; Runnable from the repository root directory with
;; guile --debug -L . test/test.scm
(define-module (mgamma tests)
  #:use-module (srfi srfi-1)
  #:use-module (srfi srfi-64)
  #:use-module ((mgamma core) #:prefix core:)
  #:use-module ((gsl vectors) #:prefix vec:)
  #:use-module ((gsl matrices) #:prefix mtx:)
  #:use-module ((lmdb lmdb) #:prefix mdb:))

(test-begin "io")
(define lines ((@@ (mgamma core) read-separated-lines)  "test/separated-lines.txt"))
(test-assert (every (lambda (l)
                      (equal? '("a" "b" "c") l))
                    lines))
(test-end "io")

(test-begin "vec-helpers")
(test-equal 2.5 ((@@ (mgamma core) vec-mean) (vec:alloc 4 #(1 2 3 4))))
(test-equal 0 ((@@ (mgamma core) vec-mean) (vec:alloc 0 #())))
;; Should it be 0, actually?—hiding NaNs is not nice. --aartaka
(test-equal 0 ((@@ (mgamma core) vec-mean) (vec:alloc 3 #(+nan.0 +nan.0 +nan.0))))
(test-equal #(3.0 3.0 3.0)
  (vec:->vector (let ((vec (vec:alloc 3 #(+nan.0 +nan.0 +nan.0))))
                  ((@@ (mgamma core) vec-replace-nan) vec 3)
                  vec)))
(test-end "vec-helpers")

(test-begin "mtx-helpers")
(define (cleanup-and-return rows cols fill)
  (let ((mtx (mtx:alloc rows cols fill)))
    ((@@ (mgamma core) cleanup-mtx) mtx)
    (mtx:->2d-vector mtx)))
(test-equal #(#(-0.5 0.5)
              #(-0.5 0.5))
  (cleanup-and-return
   2 2 #((1 2)
         (3 4))))
(test-equal #(#(0.0 0.0)
              #(-0.5 0.5))
  (cleanup-and-return
   2 2 #((1 +nan.0)
         (3 4))))
(test-equal #(#(0.0 0.0)
              #(0.0 0.0))
  (cleanup-and-return
   2 2 #((+nan.0 +nan.0)
         (+nan.0 +nan.0))))
(test-equal #(#(0.5 -0.5)
              #(-0.5 -0.5))
  (mtx:->2d-vector
   ((@@ (mgamma core) center-matrix!)
    (mtx:alloc 2 2 #((1 2)
                     (3 4))))))
(test-end "mtx-helpers")

(test-begin "import-export")
(define geno-data
  ((@@ (mgamma core) geno.txt->genotypes-mtx)
   "test/geno.txt"))
(define geno-mtx (first geno-data))
(test-equal #(#(1.0 1.0 0.0 0.0 0.0 1.0 1.0 0.0 1.0)
              #(1.0 +nan.0 0.0 0.0 0.0 1.0 1.0 0.0 1.0)
              #(1.0 1.0 0.0 0.0 +nan.0 1.0 1.0 0.0 1.0))
  (mtx:->2d-vector geno-mtx))
(test-equal '("rs31443144" "rs6269442" "rs32285189")
  (second geno-data))
;; Don't need to test the implementation of
;; pheno.txt->pheno-mtx/covariates.txt->cvt-mtx because they are the
;; same function internally.
(define kin ((@@ (mgamma core) cxx.txt->kinship)
             "test/cxx.txt"))
(test-equal 0.2429357387 (mtx:get kin 0 0))
(test-equal 0.2589379741 (mtx:get kin 3 3))
(test-equal 0.2538312374 (mtx:get kin 6 6))

((@@ (mgamma core) kinship->cxx.txt)
 kin
 "test/another.cxx.txt")
(test-assert (mtx:equal? kin
                         ((@@ (mgamma core) cxx.txt->kinship)
                          "test/another.cxx.txt")))
(unless (file-exists? "/tmp/mgamma-kin/")
  (mkdir "/tmp/mgamma-kin/"))
((@@ (mgamma core) kinship->lmdb)
 kin
 "/tmp/mgamma-kin")
;; Can't compare with mtx:equal?, because LMDB storage format looses
;; precision due to float storage.
(define-syntax-rule (test-approx a b)
  (test-assert (< (abs (- a b)) 1e-5)))

(define lmdb-kin ((@@ (mgamma core) lmdb->kinship)
                  "/tmp/mgamma-kin/"))
(test-approximate (mtx:get kin 0 0) (mtx:get lmdb-kin 0 0) 1e-5)
(test-approximate (mtx:get kin 3 3) (mtx:get lmdb-kin 3 3) 1e-5)
(test-approximate (mtx:get kin 6 6) (mtx:get lmdb-kin 6 6) 1e-5)

(unless (file-exists? "/tmp/mgamma-geno/")
  (mkdir "/tmp/mgamma-geno/"))
((@@ (mgamma core) geno.txt->lmdb)
 "test/geno.txt"
 "/tmp/mgamma-geno")
(define lmdb-geno
  ((@@ (mgamma core) lmdb->genotypes-mtx)
   "/tmp/mgamma-geno"))
;; Can't compare with mtx:equal? again.
(define lmdb-geno-mtx (first lmdb-geno))
(test-approximate (mtx:get geno-mtx 0 0) (mtx:get lmdb-geno-mtx 0 0) 1e-5)
(test-approximate (mtx:get geno-mtx 3 3) (mtx:get lmdb-geno-mtx 3 3) 1e-5)
(test-approximate (mtx:get geno-mtx 6 6) (mtx:get lmdb-geno-mtx 6 6) 1e-5)
(test-equal 3 (length (second lmdb-geno)))
(test-end "import-export")

(test-begin "usefulness")
(define pheno-mtx (mtx:alloc 7 2 0))
(mtx:set! pheno-mtx 5 0 +nan.0)
(mtx:set! pheno-mtx 2 1 +nan.0)
(mtx:set! pheno-mtx 3 0 +nan.0)
(define useful-inds (core:useful-individuals pheno-mtx #f))
(test-equal '(#t #t #t #f #t #f #t)
  useful-inds)
(define useful-pheno ((@@ (mgamma core) useful-pheno-mtx)
                      pheno-mtx
                      useful-inds))
(test-equal 5 (mtx:rows useful-pheno))
;; Because it's not on the first column, this NaN is retained. GEMMA
;; does that too? --aartaka
(test-equal +nan.0 (mtx:get useful-pheno 2 1))
(test-equal 0.0 (mtx:get useful-pheno 0 0))
(test-equal 0.0 (mtx:get useful-pheno 3 1))
(test-equal 0.0 (mtx:get useful-pheno 1 0))
(define useful-kin ((@@ (mgamma core) useful-kinship-mtx)
                    kin
                    useful-inds))
(test-equal '(5 5) (mtx:dimensions useful-kin))
(test-equal -0.00219780305 (mtx:get useful-kin 3 0))
(test-equal -0.01736104752 (mtx:get useful-kin 0 4))
;; Not the cleanest test there might be, because geno-mtx is 3x10,
;; while inds are 7. But it works in a more or less sensible
;; way. --aartaka
(define useful-geno ((@@ (mgamma core) useful-geno-mtx)
                     geno-mtx
                     useful-inds))
(test-equal +nan.0 (mtx:get useful-geno 1 1))
(test-equal +nan.0 (mtx:get useful-geno 2 3))
(test-equal 1.0 (mtx:get useful-geno 0 0))
(test-equal 1.0 (mtx:get useful-geno 0 4))
(test-equal 0.0 (mtx:get useful-geno 1 3))
(test-end "usefulness")
