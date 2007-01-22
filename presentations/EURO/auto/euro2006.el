(TeX-add-style-hook "euro2006"
 (lambda ()
    (LaTeX-add-bibitems
     "BIS2005"
     "Bolder1999"
     "Geyer1999"
     "Jankowitsch2004"
     "Nelson1987"
     "Svensson1994")
    (LaTeX-add-labels
     "accruedinterest"
     "bondpriceeq"
     "yield"
     "bondprceq2"
     "forwrate"
     "duration"
     "spread"
     "nelson"
     "multisv"
     "multsv")
    (TeX-run-style-hooks
     "bm"
     "amsbsy"
     "amsthm"
     "amsmath"
     "amssymb"
     "fontenc"
     "T1"
     "times"
     "latex2e"
     "beamer10"
     "beamer"
     "handout"
     "mathserif"
     "10pt")))

