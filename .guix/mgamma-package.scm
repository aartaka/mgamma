(define-module (mgamma-package)
  #:use-module (guix build-system guile)
  #:use-module (guix gexp)
  #:use-module (guix git-download)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix packages)
  #:use-module (guix utils)
  #:use-module (gnu packages base)
  #:use-module (gnu packages bash)
  #:use-module (gnu packages commencement)
  #:use-module (gnu packages databases)
  #:use-module (gnu packages guile)
  #:use-module (gnu packages guile-xyz)
  #:use-module (gnu packages maths)
  #:use-module (gnu packages pkg-config)
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
                           guile-bin)))))
                  (add-before 'build 'substitute-openblas-so
                   (lambda* (#:key inputs outputs #:allow-other-keys)
                     (let ((openblas (string-append (assoc-ref inputs "openblas")
                                                    "/lib/libopenblasp-r0.3.20.so")))
                       (substitute*
                        '("modules/mgamma/utils.scm")
                        (("libopenblas.so")
                         openblas))
                       #t)))
                  (add-before 'build 'build-libmgamma
                    (lambda* (#:key inputs outputs #:allow-other-keys)
                      (let ((extdir (string-append (assoc-ref outputs "out")
                                                   "/lib/guile/3.0/extensions")))
                        (invoke "make" "-C" "extension"
                                (string-append "CC="
                                               (assoc-ref inputs "gcc-toolchain")
                                               "/bin/gcc")
                                "libmgamma.so")
                        (mkdir-p extdir)
                        (copy-file "extension/libmgamma.so"
                                   (string-append extdir "/libmgamma.so")))))
                  (add-before 'build 'substitute-libmgamma
                    (lambda* (#:key outputs #:allow-other-keys)
                      (substitute*
                          '("modules/mgamma/utils.scm")
                        (("load-extension \"libmgamma.so\"")
                         ;; Bump the version when updating
                         (format #f "load-extension \"~a/lib/guile/3.0/extensions/libmgamma.so\""
                                 (assoc-ref outputs "out"))))))
                  (add-after 'build 'make-bin
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
                  gnu-make
                  guile-3.0
                  gcc-toolchain
                  openblas
                  guile-json-4
                  guile-lmdb-git
                  guile-gsl-git
                  guile-lapack-git
                  pkg-config))
    (home-page "https://github.com/aartaka/guile-gemma")
    (synopsis "Mgamma is a tool for Genome-Wide Association Studies (GWAS).")
    (description "This package provides a Guile @code{(mgamma core)} module and a
@code{mgamma} CLI tool. Supported use-cases and commands are:
@itemize
@item @code{convert} to convert between multiple data formats and
Mgamma-specific LMDB-based format.
@item @code{kinship} to compute kinship matrix for given genetic data.
@item @code{gwa} to compute LMM (Linear Mixed Model) parameters based
on kinship matrix and genotypes.  @end itemize")
    (license license:gpl3+)))

mgamma-git
