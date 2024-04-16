#! /usr/local/bin/guile \
--no-auto-compile -e main -s
!#

(use-modules (gemma core))

(define (main args)
  (write args)
  (newline))
