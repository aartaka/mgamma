#! /usr/local/bin/guile \
--no-auto-compile -e main
!#

(use-modules (gemma core))

(define (main args)
  (write args)
  (newline))
