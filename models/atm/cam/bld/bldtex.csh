#!/bin/csh
#
#	bldtex.csh
#
#	Create postscript and pdf file from LaTex design doc.
#
#	$Id: bldtex.csh,v 1.1.6.1 2002/05/13 17:21:23 eaton Exp $
#
set files=("Requirements" "Design");
foreach file ( $files )
  echo "Build the $file file";
  latex $file.tex;
  latex $file.tex;
  ./dotoc $file.tex;
  bibtex $file.tex;
  latex $file.tex;
  dvips -o $file.ps $file.dvi;
  dvipdf $file.dvi $file.pdf;
end
