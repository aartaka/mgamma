#! /usr/local/bin/guile \
--no-auto-compile -e main -s
!#

(use-modules (mgamma core))
(use-modules (ice-9 getopt-long))
(use-modules (ice-9 format))
(use-modules (ice-9 match))
(use-modules (srfi srfi-1))
(use-modules (srfi srfi-26))

(define (ensure-options options . names)
  (let rec ((names names))
    (unless (null? names)
      (unless (option-ref options (car names) #f)
        (error (format #t "Missing mandatory --~a option~%" (car names))))
      (rec (cdr names)))))

(define regular-options
  `((help    (single-char #\h) (value #f))
    (kinship (single-char #\k) (value #t) (predicate ,file-exists?))
    (geno    (single-char #\g) (value #t) (predicate ,file-exists?))
    (pheno   (single-char #\p) (value #t) (predicate ,file-exists?))
    (anno    (single-char #\a) (value #t) (predicate ,file-exists?))
    (covar   (single-char #\c) (value #t) (predicate ,file-exists?))
    (maf                       (value #t) (predicate ,string->number))
    (map-size                  (value #t))
    (output  (single-char #\o) (value #t))))

(define (mdb-file? file)
  (or (string-suffix? ".mdb" file)
      (string-suffix? ".lmdb" file)))

(define (txt-file? file)
  (or (string-suffix? ".out" file)
      (string-suffix? ".txt" file)))

(define (convert options)
  (unless (option-ref options 'help #f)
    (ensure-options options 'output))
  (let* ((help (option-ref options 'help #f))
         (output (option-ref options 'output #f))
         (mdb-out? (and output (mdb-file? output)))
         (txt-out? (and output (txt-file? output))))
    (cond
     (help
      (format #t "Convert files from one format to another:
geno txt     -> lmdb:   [--map-size 10000000] --geno geno.txt --output data.mdb
kinship txt  -> lmdb:   [--map-size 10000000] --kinship kinship.txt --output data.mdb
kinship lmdb -> txt     [--map-size 10000000] --kinship data.mdb --output kinship.txt~%"))
     ((and mdb-out?
           (option-ref options 'geno #f))
      (geno.txt->lmdb (option-ref options 'geno #f) (dirname output)))
     ((and mdb-out?
           (option-ref options 'kinship #f))
      (kinship->lmdb (cxx.txt->kinship (option-ref options 'kinship #f))
                     (dirname output)))
     ((and txt-out?
           (option-ref options 'kinship #f))
      (let ((old-kinship-lmdb (option-ref options 'kinship #f))
            (mgamma-data-dir "/tmp/mgamma-data/"))
        (unless (file-exists? mgamma-data-dir)
          (mkdir mgamma-data-dir))
        (copy-file old-kinship-lmdb (string-append mgamma-data-dir "/data.mdb"))
        (kinship->cxx.txt (lmdb->kinship mgamma-data-dir)
                          output)
        (delete-file "/tmp/mgamma-data/data.mdb")))
     (else
      (error "Cannot convert between these formats (yet?)")))
    (when (and mdb-out?
               (not help))
      (rename-file (string-append (dirname output) "/data.mdb")
                   output))))

(define (suffixed-string->number string)
  (* (or (string->number string)
         (string->number (substring string 0 (1- (string-length string)))))
     (cond
      ((string-suffix-ci? "K" string) 1000)
      ((string-suffix-ci? "M" string) 1000000)
      ((string-suffix-ci? "G" string) 1000000000)
      (else 1))))

(define (kinship options)
  (unless (option-ref options 'help #f)
    (ensure-options options 'output 'geno 'pheno))
  (let* ((help (option-ref options 'help #f))
         (output (option-ref options 'output #f))
         (mdb-out? (and output (mdb-file? output)))
         (txt-out? (and output (txt-file? output))))
    (cond
     (help
      (format #t "Compute kinship matrix based on the genotype and phenotype files:
mgamma kinship [--maf 0.1] [--map-size 10000000] --geno geno.(lmdb|txt) --pheno pheno.txt --output kinship.(txt|mdb)~%"))
     (else
      (let ((map-size (option-ref options 'map-size #f)))
        (when map-size
          (mapsize (suffixed-string->number map-size))))
      (let* ((geno (option-ref options 'geno #f))
             (geno-txt? (txt-file? geno)))
        (match (if geno-txt?
                   (geno.txt->genotypes-mtx geno)
                   (lmdb->genotypes-mtx (dirname geno)))
          ((geno-mtx markers)
           (let* ((useful-snps (useful-snps geno-mtx markers (option-ref options 'pheno #f)
                                            #:maf-level (or (string->number (option-ref options 'maf ""))
                                                            0.01)))
                  (_ (cleanup-mtx geno-mtx))
                  (kin-mtx (kinship-mtx geno-mtx (hash-count (cut or #t <> <>) useful-snps))))
             (if txt-out?
                 (kinship->cxx.txt kin-mtx output)
                 (kinship->lmdb kin-mtx (dirname output)))))))))
    (when (and mdb-out?
               (not help))
      (rename-file (string-append (dirname output) "/data.mdb")
                   output))))

(define (main args)
  (let* ((command (if (> (length args) 1)
                      (second args)
                      #f))
         (help (or (not command)
                   (string= "--help" command)
                   (string= "help" command)
                   (string= "-h" command)))
         (command (if command
                      (string->symbol command)
                      #f))
         (options (if help
                      #f
                      (getopt-long (cdr args) regular-options))))
    (cond
     ((or help (not command))
      (unless command
        (format #t "Provide a valid command~%"))
      (format #t "mgamma is a Genome-Wide Association tool
Usage: mgamma command [options...]

Valid commands:
 convert:  Convert the text files for --kinship or --geno to LMDB versions and back.
 kinship:  Calculate kinship matrix for --geno and --pheno.
 gwa:      Run LMM9 on the --kinship matrix and --geno/--pheno files.

Use `mgamma convert --help' to get more detailed help on convert.
Same for kinship and gwa.~%"))
     (command
      (case command
        ((convert)
         (convert options))
        ((kinship)
         (kinship options))
        ((gwa)
         (format #t "Not implemented yet")))))))
