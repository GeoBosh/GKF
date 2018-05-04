#!/bin/bash

R CMD Sweave vignette.Rnw

pdflatex vignette.tex

if [ $# -ge 1 ]
then
   bibtex vignette
   pdflatex vignette.tex
   pdflatex vignette.tex
   echo "No display !" 
else
  evince vignette.pdf &
fi
