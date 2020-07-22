#!/usr/bin/gnuplot -persist
set xtics ( "Gamma" 0,  "X" 57,  "W" 86,  "L" 115,  "Gamma" 185,  "K" 259 )
unset key
nRows = real(system("awk '$1==\"kpoint\" {nRows++} END {print nRows}' bandstruct.kpoints"))
nCols = real(system("wc -c < bandstruct.eigenvals")) / (8*nRows)
formatString = system(sprintf("echo '' | awk 'END { str=\"\"; for(i=0; i<%d; i++) str = str \"%%\" \"lf\"; print str}'", nCols))
plot for [i=1:nCols] "bandstruct.eigenvals" binary format=formatString u 0:i w l
