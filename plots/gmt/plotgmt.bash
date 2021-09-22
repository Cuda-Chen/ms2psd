#!/bin/bash

set -eux

PDFDATA=../../pdf_out.txt
PDFPERIODS=../../center_periods_out.txt
HIGHNOISEMODEL=../../highnoise.mod
LOWNOISEMODEL=../../lownoise.mod
COLORMAP=psdpdf.cpt

PDFGRID=pdf.nc
PLOTFILE=plot.ps

PERIOD_MAX=`cat ${PDFPERIODS} | awk 'BEGIN {max=0} {if($1>max) max=$1} END {print log(max)/log(10)}'`
PERIOD_MIN=`cat ${PDFPERIODS} | awk 'BEGIN {min=10000} {if($1<min) min=$1} END {print log(min)/log(10)}'`
PERIODS=(`cat ../../center_periods_out.txt`)
x=${PERIODS[0]}
y=${PERIODS[1]}
PERIOD_INTERVAL_LOG_SCALE=`awk -v var1=${x} -v var2=${y} 'BEGIN{print log(var2)/log(10) - log(var1)/log(10)}'`
COLORMAP_ADJUST=psdpdf_adjust.cpt

# Create netCDF grid
gmt xyz2grd ${PDFDATA} -R${PERIOD_MIN}/${PERIOD_MAX}/-200/-50 -I${PERIOD_INTERVAL_LOG_SCALE}/1 -G${PDFGRID} -ZTLA -V

# Get min/max of PDF
PDF_MIN=`gmt grdinfo pdf.nc | grep v_min | awk '{print $3}'`
PDF_MAX=`gmt grdinfo pdf.nc | grep v_min | awk '{print $5}'`

# Adjust color palette according to min/max of PDF
gmt makecpt -C${COLORMAP} -T${PDF_MIN}/${PDF_MAX} > ${COLORMAP_ADJUST}

gmt begin example
    # Plot PDF
    gmt grdview ${PDFGRID} -JX6i/5i -Bxa1.0+l"log10(Period)" -Bya10f5+l"Power [10log10(m**2/sec**4/Hz)] [dB]" -BWsne -C${COLORMAP_ADJUST} -S100 -Qs -N0 -V -Y4.0
    #gmt colorbar -C${COLORMAP_ADJUST} -B.02 -D6.15i/2.5i/5.0i/0.25i -V
    gmt colorbar -C${COLORMAP_ADJUST} -Li0.05 -D6.15i/2.5i/5.0i/0.25i -V
    #gmt psconvert ${PLOTFILE} -A -Tf -V
gmt end show
