(define-module (mgamma useful)
  #:use-module (srfi srfi-1)
  #:use-module ((gsl vectors) #:prefix vec:)
  #:use-module ((gsl matrices) #:prefix mtx:)
  #:export (useful-individuals
            useful-snps
            useful-kinship-mtx
            useful-geno-mtx
            useful-geno-mtx-per-snps
            useful-pheno-mtx))

(define (useful-individuals pheno-mtx cvt-mtx)
  "Return a list of booleans marking individuals in PHENO-MTX.
#t is useful/indicator individuals.
#f is individuals with incomplete data and otherwise not useful."
  (let* ((inds (do ((ind 0 (1+ ind))
                    (inds '()
                          (cons (not (nan? (mtx:get pheno-mtx ind 0)))
                                inds)))
                   ((= ind (mtx:rows pheno-mtx))
                    (reverse! inds)))))
    (when (and cvt-mtx
               (positive? (mtx:cols cvt-mtx)))
      (do ((inds inds (cdr inds))
           (cvt 0 (1+ cvt)))
          ((null? inds))
        (when (zero? (mtx:get cvt-mtx cvt 0))
          (set-car! inds #f))))
    inds))

(define* (useful-snps genotypes-mtx markers pheno-mtx covariates-mtx #:key (miss-level 0.05) (maf-level 0.01))
  "Return a hash-table from SNP names to MAF values.
Only includes the useful/indicator SNPs.
Filter the SNPs for MAF being in (MAF-LEVEL, 1 - MAF-LEVEL) range and
having less than MISS-LEVEL missing values."
  (let* ((useful-inds (useful-individuals pheno-mtx covariates-mtx))
         (ind-count (length useful-inds))
         (useful-snp-table (make-hash-table))
         (mtx-rows (mtx:rows genotypes-mtx))
         (mtx-cols (mtx:columns genotypes-mtx)))
    (do ((row 0 (1+ row))
         (markers markers (cdr markers)))
        ((= row mtx-rows))
      (let* ((name (car markers))
             (maf (do ((ind 0 (1+ ind))
                       (useful-inds useful-inds (cdr useful-inds))
                       (maf 0 (if (car useful-inds)
                                  (+ maf (mtx:get genotypes-mtx row ind))
                                  maf)))
                      ((= ind mtx-cols) maf)))
             (miss-count (do ((ind 0 (1+ ind))
                              (useful-inds useful-inds (cdr useful-inds))
                              (nans 0 (if (and (car useful-inds)
                                               (nan? (mtx:get genotypes-mtx row ind)))
                                          (1+ nans)
                                          nans)))
                             ((= ind mtx-cols) nans)))
             (maf (/ maf (* 2 (- ind-count miss-count)))))
        (when (and (< (/ miss-count ind-count) miss-level)
                   (< maf-level maf (- 1 maf-level)))
          (hash-set! useful-snp-table name maf))))
    useful-snp-table))

(define (useful-kinship-mtx kinship-mtx useful-inds)
  "Only retain these individuals in the KINSHIP-MTX that are USEFUL-INDS.
Create and return a new matrix of size equal to useful individuals
number."
  (let* ((n-useful (count identity useful-inds))
         (new-mtx (mtx:alloc n-useful n-useful)))
    (do ((outer-useful-inds useful-inds (cdr outer-useful-inds))
         (kin-i 0 (1+ kin-i))
         (i 0 (if (car outer-useful-inds)
                  (1+ i)
                  i)))
        ((null? outer-useful-inds))
      (when (car outer-useful-inds)
        (do ((inner-useful-inds useful-inds (cdr inner-useful-inds))
             (kin-j 0 (1+ kin-j))
             (j 0 (if (car inner-useful-inds)
                      (1+ j)
                      j)))
            ((null? inner-useful-inds))
          (when (car inner-useful-inds)
            (mtx:set! new-mtx i j
                      (mtx:get kinship-mtx kin-i kin-j))))))
    new-mtx))

(define (useful-geno-mtx geno-mtx useful-inds)
  "Remove non-USEFUL-INDS from GENO-MTX.
Create and return a new matrix."
  (let* ((n-useful (count identity useful-inds))
         (new-mtx (mtx:alloc (mtx:rows geno-mtx) n-useful)))
    (do ((row 0 (1+ row)))
        ((= row (mtx:rows new-mtx)))
      (do ((useful-inds useful-inds (cdr useful-inds))
           (geno-i 0 (1+ geno-i))
           (i 0 (if (car useful-inds)
                    (1+ i)
                    i)))
          ((null? useful-inds))
        (when (car useful-inds)
          (mtx:set! new-mtx row i
                    (mtx:get geno-mtx row geno-i)))))
    new-mtx))

(define (useful-geno-mtx-per-snps geno-mtx markers useful-snps)
  "Only retain the USEFUL-SNPS in the GENO-MTX (ordered by MARKERS).
Return a new matrix with cleaned-up SNPs."
  (vec:with
   (tmp (mtx:columns geno-mtx) 0)
   (let* ((new-rows (hash-count (lambda (k v) v) useful-snps))
          (new-mtx (mtx:alloc new-rows (mtx:columns geno-mtx) 0))
          (i 0))
     (do ((markers markers (cdr markers))
          (row 0 (1+ row)))
         ((= i new-rows))
       (when (hash-ref useful-snps (car markers) #f)
         (mtx:row->vec! geno-mtx row tmp)
         (mtx:vec->row! tmp new-mtx i)
         (set! i (1+ i))))
     new-mtx)))

(define (useful-pheno-mtx pheno-mtx useful-inds pheno-nums)
  "Only retain USEFUL-INDividualS in PHENO-MTX.
Also only retain the PHENO-NUMS-numbered columns.
Return a new matrix with cleaned-up ones."
  (let* ((n-useful (count identity useful-inds))
         (new-mtx (mtx:alloc n-useful (length pheno-nums))))
    (do ((useful-inds useful-inds (cdr useful-inds))
         (row 0 (1+ row))
         (i 0))
        ((null? useful-inds))
      (when (car useful-inds)
        (do ((col 0 (1+ col))
             (j 0))
            ((= col (mtx:columns pheno-mtx)))
          (when (member col pheno-nums)
            (mtx:set! new-mtx i j
                      (mtx:get pheno-mtx row col))
            (set! j (1+ j))))
        (set! i (1+ i))))
    new-mtx))
