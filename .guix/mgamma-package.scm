(define-module (mgamma-package)
 #:use-module (gnu packages guile)
 #:use-module (gnu packages guile-xyz)
 #:use-module (gnu packages databases)
 #:use-module (gnu packages bash)
 #:use-module (guix packages)
 #:use-module (guix gexp)
 #:use-module (guix utils)
 #:use-module (guix build-system guile)
 #:use-module (guix git-download)
 #:use-module ((guix licenses) #:prefix license:)
 #:use-module (guile-gsl-package)
 #:use-module (guile-lmdb-package)
 #:use-module (guile-lapack-package))

(define-public mgamma-git
  (package
    (name "mgamma-git")
    (version "0.0.1")
    (source (local-file ".."
                        "mgamma-git-checkout"
                        #:recursive? #t
                        #:select? (or (git-predicate (dirname (current-source-directory)))
                                      (const #t))))
    (build-system guile-build-system)
    (arguments
     '(#:source-directory "modules"
       #:phases (modify-phases %standard-phases
                  (add-before 'build 'substitute-guile
                    (lambda* (#:key inputs #:allow-other-keys)
                      (let ((guile-bin (string-append (assoc-ref inputs "guile")
                                                      "/bin/guile")))
                        (substitute*
                            '("bin/mgamma")
                          (("/usr/local/bin/guile")
                           guile-bin))
                        #t)))
                  (add-after 'install 'make-bin
                    (lambda* (#:key inputs outputs #:allow-other-keys)
                      (let* ((bin (string-append (assoc-ref outputs "out")
                                                 "/bin"))
                             (real (string-append bin "/mgamma"))
                             (bash (string-append (assoc-ref inputs "bash")
                                                  "/bin/bash"))
                             (self (string-append (assoc-ref outputs "out")
                                                  "/share/guile/site/3.0")))
                        (mkdir-p bin)
                        (copy-file "bin/mgamma" real)
                        (wrap-program real #:sh bash
                                      ;; For better backtraces.
                                      (list "COLUMNS" ":" '= (list "1000"))
                                      (list "GUILE_LOAD_PATH" ":" '= (list self (getenv "GUILE_LOAD_PATH"))))))))))
    (native-inputs (list guile-3.0))
    (inputs (list bash
                  guile-3.0
                  guile-json-4
                  guile-lmdb-git
                  guile-gsl-git
                  guile-lapack-git))
    (home-page "https://github.com/aartaka/guile-gemma")
    (synopsis "")
    (description "")
    (license license:gpl3+)))

mgamma-git
