#! /usr/local/bin/guile \
--no-auto-compile -e main -s
!#

(use-modules (mgamma core))
(use-modules (ice-9 getopt-long))
(use-modules (ice-9 format))
(use-modules (srfi srfi-1))

(define regular-options
  `((help    (single-char #\h) (value #f))
    (kinship (single-char #\k) (value #f))
    (geno    (single-char #\g) (value #t) (predicate ,file-exists?))
    (pheno   (single-char #\p) (value #t) (predicate ,file-exists?))
    (lmdb    (single-char #\l) (value #t))
    (output  (single-char #\o) (value #t))))

(define (main args)
  (let* ((options (getopt-long args regular-options))
         (help (or (option-ref options 'help #f)
                   (null? (cdr args))))
         (kinship (option-ref options 'kinship #f))
         (geno.txt (option-ref options 'geno #f))
         (pheno.txt (option-ref options 'pheno #f))
         (lmdb-dir (or (option-ref options 'lmdb #f)
                       (and geno.txt
                            (string-append "/tmp/" (basename geno.txt) "-lmdb/"))))
         (output (option-ref options 'output #f))
         (help-fn (lambda ()
                    (format #t "~&~a [options...]
  -h, --help               Show this help
  -k, --kinship            Kinship matrix computation
  -g, --geno    geno.txt   Genotype file
  -p, --pheno   pheno.txt  Phenotype file
  -l, --lmdb    /lmdb/dir  LMDB directory to cache genotype/kinship data in
  -o, --output  file       Output file for kinship/LMM results

Kinship matrix

 When --kinship is provided, compute kinship matrix for individuals in dataset.
 Requires --geno and --pheno.
 Puts the results in --output file.
 Will be faster with --lmdb (see below)

LMDB genotype DB generation

 All the operations are significantly faster with genotypes in LMDB database.
 To generate the database, provide --geno and --lmdb options without --kinship.
 Provide a directory path to --lmdb, and the DB will be put there.
 In case this step is omitted, database is generated anyway on --kinship step.~%"
                            (first args)))))
    (cond
     (help
      (help-fn))
     ((and kinship geno.txt pheno.txt output)
      (let ((kinship (kmain geno.txt pheno.txt lmdb-dir)))
        (kinship->cxx.txt kinship output)))
     ((and geno.txt lmdb-dir)
      (geno.txt->lmdb geno.txt lmdb-dir))
     (else
      (format #t "~&Unrecognized combination of CLI options.~%~%")
      (help-fn)))))
