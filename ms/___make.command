#!/bin/bash


echo -e "\n\n~~~~~~~~~~~~~~~~~~\n"
echo Creating the LaTeX document...
echo -e "\n~~~~~~~~~~~~~~~~~~\n"

echo -e "\n\n\n"
echo Should I update just locally? \(yes/no, defaults to no\)
echo -e "\n"

read local

echo -e "\n"



if [ -z "$local" ]; then
    local="no"
fi
if [ "$local" == "y" ]; then
    local="yes"
fi
if [ "$local" == "n" ]; then
    local="no"
fi

if [ "$local" != "yes" ] && [ "$local" != "no" ]; then
    echo -e "\n\n"
    echo ERROR: invalid input
    echo -e "\n\n"
    exit
fi



if [ "$local" = "yes" ]; then
    echo Okay, I\'ll just update the version in your \`trans_trends/ms\` directory.
    echo -e "\n"
fi
if [ "$local" = "no" ]; then
    echo Okay, I\'ll update the version in your \`trans_trends/ms\` directory,
    echo plus the Box version
    echo -e "\n"
fi


sleep 2





## This creates the PDF from the __ms.tex document.
## On macOS, you should just be able to click on it from Finder

export name="__ms"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd "$DIR"

pdflatex -interaction=nonstopmode ${name}
bibtex ${name}
pdflatex -interaction=nonstopmode ${name}
pdflatex -interaction=nonstopmode ${name}
rm -Rf *.aux *.bbl *.bcf *.blg *.log *.out *.run.xml *.synctex.gz



## If the proper directory exists, copy the new manuscript to a folder on Box
## so that it can be shared:



if [ "$local" = "no" ]; then
    if [ -d ~/"Box Sync/trans_trends/" ]; then
        cp __ms.pdf ~/"Box Sync/trans_trends/" && \
        mv ~/"Box Sync/trans_trends/__ms.pdf" ~/"Box Sync/trans_trends/manuscript.pdf"
    fi
fi



exit
