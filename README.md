# mgamma does GWA in Guile

This is a work-in-progress attempt to reproduce some of
[GEMMA](https://github.com/genetics-statistics/GEMMA) functionality and provide next generation GWA tooling for the GeneNetwork project.

IMPORTANT: WIP. YMMV.

# Usage

Get help and example usage:
```
mgamma
mgamma help
mgamma -h
mgamma --help
mgamma convert --help
mgamma kinship --help
mgamma gwa --help
```

Convert a textual genotype file to mgamma lmdb format

```
mgamma convert -g genofile -a annofile [--map-size 10M] -o geno.mdb
```

Use this to compute the kinship matrix

```
mgamma kinship [--maf 0.1] [--map-size 10M] -g geno.lmdb -p pheno.txt -o kinship.mdb
```

To compute the LMM we use

```
mgamma gwa -g geno.mdb -p pheno.txt [--map-size 10M] -k kinship.mdb -o assoc.txt
```

(note the annotations are now part of the geno file and we use GEMMA's -lmm 9 by default).

TODO: When no file names are given we default to above names, so:

```
mgamma kinship -o kinship.mdb
mgamma gwa -o assoc.txt
```

should compute the same kinship.mdb and assoc.txt files from default geno/pheno files.

## LOCO

In the near future we'll support LOCO to compute kinship and SNPs by default. To skip LOCO use:

```
mgamma kinship -g geno.lmdb -o kinship.lmdb --no-loco
mgamma gwa -g geno.lmdb -p pheno.txt -k kinship.lmdb -o gwa.lmdb
```

will compute the matrix without LOCO because it is missing in the kinship file.

## Covariates

To add covariates run the LMM command with a covariate file

```
mgamma gwa -g geno.lmdb -p pheno.txt -c covariates.txt -k kinship.lmdb
```

## Extras

Convert a kinship file to its (old) tab delimited information

```
mgamma convert -k kinship.lmdb -o kinship.out
```
