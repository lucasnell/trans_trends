#!/bin/bash

## This creates the PDF from the __ms.tex document.
## On macOS, you should just be able to click on it from Finder

export name="__ms"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd "$DIR"

pdflatex -interaction=nonstopmode ${name}
biber ${name}
pdflatex -interaction=nonstopmode ${name}
pdflatex -interaction=nonstopmode ${name}
rm -Rf *.aux *.bbl *.bcf *.blg *.log *.out *.run.xml *.synctex.gz



## If the proper directory exists, copy the new manuscript to a folder on Box
## so that it can be shared:

[ -d ~/"Box Sync/trans_trends/" ] && \
  cp __ms.pdf ~/"Box Sync/trans_trends/" && \
  mv ~/"Box Sync/trans_trends/__ms.pdf" ~/"Box Sync/trans_trends/manuscript.pdf"


exit
