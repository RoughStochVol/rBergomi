(TeX-add-style-hook
 "documentation"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "amsmath"
    "amsfonts")
   (LaTeX-add-labels
    "sec:techn-docum-verbmc_b"
    "sec:background"
    "sec:implementation"))
 :latex)

