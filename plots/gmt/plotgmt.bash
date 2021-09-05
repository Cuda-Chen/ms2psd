#!/bin/bash

PDFDATA=../../pdf_out.txt
PDFPERIODS=../../center_periods_out.txt
HIGHNOISEMODEL=../../highnoise.mod
LOWNOISEMODEL=../../lownoise.mod
COLORMAP=psdpdf.cpt

PDFGRID=pdf.grd
PLOTFILE=plot.ps

gmt begin example
    gmt xyz2grd ${PDFDATA} -R0/102/0/150 -I1/1 -G${PDFGRID} -V
    gmt grdview ${PDFGRID} -JX6i/5i -Bxa1.0+l"log10(Period)" -Bya10f5+l"Power [10log10(m**2/sec**4/Hz)] [dB]" -BWsne -C${COLORMAP} -S100 -Qs -N0 -V -Y4.0
    gmt colorbar -C${COLORMAP} -B.02 -D6.15i/2.5i/5.0i/0.25i -V
    #gmt psconvert ${PLOTFILE} -A -Tf -V
gmt end show
