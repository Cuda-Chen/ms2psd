#!/bin/bash

PDFDATA=../../pdf_out.txt
PDFPERIODS=../../center_periods_out.txt
HIGHNOISEMODEL=../../highnoise.mod
LOWNOISEMODEL=../../lownoise.mod

PDFGRID=pdf.grd
PLOTFILE=plot.ps

gmt xyz2grd ${PDFDATA} -R0/102/0/150 -I1/1 -G${PDFGRID}

gmt grdview ${PDFGRID} -JX6i/5i -Ba1.0:"log10(Period)":/a10f5:"Power [10log10(m**2/sec**4/Hz)] [dB]":Wsne -pdf testplot
