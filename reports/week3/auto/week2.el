(TeX-add-style-hook
 "week2"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper" "10pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("babel" "english") ("inputenc" "utf8") ("fontenc" "T1")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "introduction"
    "results"
    "article"
    "art10"
    "babel"
    "inputenc"
    "graphicx"
    "amsmath"
    "amssymb"
    "amsthm"
    "hyperref"
    "url"
    "listings"
    "verbatim"
    "color"
    "fontenc"
    "csvsimple"
    "array"
    "booktabs"
    "longtable")
   (TeX-add-symbols
    '("figureFromFile" 3)
    '("verbatimtable" 2)
    '("conj" 1))
   (LaTeX-add-labels
    "tab:#1"
    "fig:#1")))

