#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#  
#  This script creates plots for find_distortion.c mode 3 and 4
#
#  Oct 28, 2008 Edward J. Los - Initial version
#  Dec 15, 2008 Edward J. Los - Use 50 x 50 grids
# 
#  Usage ./plot_distortion.csh  <basename> <series>
#             <basename> is the file prefix
#             <series> is the series, or a string of multiple series
#
set basename = $1
set serieslist = $2 
date
source ${DASCH_SCRIPTS}/initgmt
foreach series ($2)
  cd ${DASCH_MATCH}
  pwd
  echo "Processing series ${series}"

  column -b -i ${basename}_${series}_count.db dix diy count > ${basename}_${series}_count.txt
  minmax ${basename}_${series}_count.txt
  xyz2grd ${basename}_${series}_count.txt -R0/50/0/50 -I1.0/1.0 -G${basename}_${series}_count.grd
  grdview ${basename}_${series}_count.grd -JX4i/4i  -JZ2i -E225/30 -B10:"X":/10:"Y":/50:"Z":WSZne >! ${basename}_${series}_count.ps
  
  column -b -i ${basename}_${series}_median.db dix diy dxmedian > ${basename}_${series}_dxmedian.txt
  minmax ${basename}_${series}_dxmedian.txt
  xyz2grd ${basename}_${series}_dxmedian.txt -R0/50/0/50 -I1.0/1.0 -G${basename}_${series}_dxmedian.grd
  grdview ${basename}_${series}_dxmedian.grd -JX4i/4i  -JZ2i -E225/30 -B10:"X":/10:"Y":/1:"Z":WSZne >! ${basename}_${series}_dxmedian.ps
  
  column -b -i ${basename}_${series}_median.db dix diy dymedian > ${basename}_${series}_dymedian.txt
  minmax ${basename}_${series}_dymedian.txt
  xyz2grd ${basename}_${series}_dymedian.txt -R0/50/0/50 -I1.0/1.0 -G${basename}_${series}_dymedian.grd
  grdview ${basename}_${series}_dymedian.grd -JX4i/4i  -JZ2i -E225/30 -B10:"X":/10:"Y":/1:"Z":WSZne >! ${basename}_${series}_dymedian.ps
  
  
  row '(dxrms < 90.0)' < ${basename}_${series}_rms.db |column -b  dix diy dxrms > ${basename}_${series}_dxrms.txt
  minmax ${basename}_${series}_dxrms.txt
  xyz2grd ${basename}_${series}_dxrms.txt -R0/50/0/50 -I1.0/1.0 -G${basename}_${series}_dxrms.grd
  grdview ${basename}_${series}_dxrms.grd -JX4i/4i  -JZ2i -E225/30 -B10:"X":/10:"Y":/1:"Z":WSZne >! ${basename}_${series}_dxrms.ps

  row '(dyrms < 90.0)' < ${basename}_${series}_rms.db |column -b  dix diy dyrms > ${basename}_${series}_dyrms.txt
  minmax ${basename}_${series}_dyrms.txt
  xyz2grd ${basename}_${series}_dyrms.txt -R0/50/0/50 -I1.0/1.0 -G${basename}_${series}_dyrms.grd
  grdview ${basename}_${series}_dyrms.grd -JX4i/4i  -JZ2i -E225/30 -B10:"X":/10:"Y":/1:"Z":WSZne >! ${basename}_${series}_dyrms.ps


  rm ${basename}_${series}_count.txt
  rm ${basename}_${series}_dxmedian.txt
  rm ${basename}_${series}_dxrms.txt
  rm ${basename}_${series}_dymedian.txt
  rm ${basename}_${series}_dyrms.txt
  rm ${basename}_${series}_count.grd
  rm ${basename}_${series}_dxmedian.grd
  rm ${basename}_${series}_dxrms.grd
  rm ${basename}_${series}_dymedian.grd
  rm ${basename}_${series}_dyrms.grd

end

date
