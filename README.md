# mgamma does GWA in Guile

This is a work-in-progress attempt to reproduce some of
[GEMMA](https://github.com/genetics-statistics/GEMMA) functionality and provide next generation GWA tooling for the GeneNetwork project.

IMPORTANT: WIP. YMMV.

# Usage

Get help and example usage:
``` sh
mgamma
mgamma help
mgamma -h
mgamma --help
mgamma convert --help
mgamma kinship --help
mgamma gwa --help
```

Convert a textual genotype file to mgamma lmdb format

``` sh
mgamma convert -g genofile -a annofile [--map-size 10M] -o geno.mdb
```

Use this to compute the kinship matrix

``` sh
mgamma kinship [--maf 0.1] [--map-size 10M] -g geno.mdb -p pheno.txt -o kinship.mdb
```

To compute the LMM we use

``` sh
mgamma gwa -g geno.mdb -p pheno.txt [--map-size 10M] -k kinship.mdb -o assoc.txt
```

(note the annotations are now part of the geno file and we use GEMMA's -lmm 9 by default).

TODO: When no file names are given we default to above names, so:

``` sh
mgamma kinship -o kinship.mdb
mgamma gwa -o assoc.txt
```

should compute the same kinship.mdb and assoc.txt files from default geno/pheno files.

## Example session

Using shell to get LMM params for given geno/pheno.

``` sh
# Convert geno to Mgamma LMDB-based format.
mgamma convert -g /path/to/mgamma/example/BXD_geno.txt -o /tmp/mgamma-geno/data.mdb
# Compute kinship matrix from geno&pheno.
mgamma kinship -g /tmp/mgamma-geno/data.mdb -p /path/to/mgamma/example/BXD_pheno.txt
# Compute LMM params from geno, pheno, and kinship.
mgamma gwa -g /tmp/mgamma-geno/data.mdb -p /path/to/mgamma/example/BXD_pheno.txt -k /tmp/mgamma-kin/data.mdb -o BXD.assoc.txt
```

or, in Scheme REPL (more verbose, but allows playing with data interactively).

``` scheme
(use-modules ((mgamma core) #:prefix mgamma:))
(use-modules (srfi srfi-1))
(define geno+markers (mgamma:geno.txt->lmdb "/path/to/mgamma/example/BXD_geno.txt" "/tmp/mgamma-geno/"))
(define geno-mtx (first geno+markers))
(define markers (second geno+markers))
(define pheno-mtx (pheno.txt->pheno-mtx "/path/to/mgamma/example/BXD_pheno.txt"))
(define kinship (kinship-mtx geno-mtx geno-markers (useful-snps geno-mtx geno-markers pheno-mtx #f)))
(define params (analyze geno-mtx geno-markers kinship #f pheno-mtx #f))
(snp-params->assoc.txt params "BXD.assoc.txt")
```

## LOCO

In the near future we'll support LOCO to compute kinship and SNPs by default. To skip LOCO use:

```
mgamma kinship -g geno.mdb -o kinship.mdb --no-loco
mgamma gwa -g geno.mdb -p pheno.txt -k kinship.mdb -o gwa.mdb
```

will compute the matrix without LOCO because it is missing in the kinship file.

## Covariates

To add covariates run the LMM command with a covariate file

```
mgamma gwa -g geno.mdb -p pheno.txt -c covariates.txt -k kinship.lmdb
```

## Extras

Convert a kinship file to its (old) tab delimited information

```
mgamma convert -k kinship.mdb -o kinship.out
```
