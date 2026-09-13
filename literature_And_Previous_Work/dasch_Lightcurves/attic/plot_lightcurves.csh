#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
####################### lc_scatter_flags ############################################################################
# Silas Laycock 2007
# This script generates a summary table containing lightcurve parameters for all lightcurves in the local directory.
# There is one row per star, and computed quantities include Median magnitude, RMS, total number of points in lightcurve, 
# number of good points, brightest and faintest magnitude recorded for that star... etc.
# A Quality criterion/threshold is applied to the points included in these statistics. This is done so that unreliable measurements 
# do not contaminate the values. Typical criteria are related to distance from center of plate, photometric fitting error etc...
# In addition to the summary table, lightcurves are written out in a plottable (selected columns, simple ascii) format with 
# a special plot-symbol column that denotes different classes of data by the use of different plotting symbols.
# The symbol codes are taken from the WIP (pgplot) manual.
# For each lightcurve a plotting template is also written, these are WIP macros and can be plotted by executing 'wip lc_OBJECT.wip' 
# If WIP is installed on the system just uncomment the line near the end of this script and the plotting will be done at run-time.
#
# Modified Sept 13th 2007. (Version "local7" to match extract_lightcurve_flags_local7.csh)
# Corrections made to compound IF statement used for assigning plot symbols.
# New flags have been introduced upstream in the local-calib and lightcurve extraction scripts 
# for identifying detected/non-detected stars (detected=0/1) and tracking local-calibration status (cal_local=1 done or 0=not-done)
# These flags are utilized by this version of the script. Purpose is unambiguous seperation of good, 'bad' & upper-limit points.
# Rationalized the Quality criterion into string variable "CONDITIONS"
# Simplified the long compound IF statement by using the new flags.
# Modified Sept 17th 2007
# Added long planned features: 	line on lower panel to show number of neighboring stars used for local calibration
# 				 color-coded blocks to indicate plate-series.
######################################################################################################################
#
# Sep 11, 2007: Edward Los  Add ( $#rms_iso == 0) error check; select first rms entry in the statstable output
# Sep 20, 2007: Edward Los  Merge with lc_scatter_flags_local7.csh
# Nov  6, 2007  Edward J. Los - Remove directory references and use environment variables instead
# Feb 25, 2008: Sumin Tang - Add sigma-clipping for DASCH median, and increase the threshold for median calculation from 3 to 5
# Mar 18, 2008: Edward J. Los - Most of the contents of this file have moved to extract_lightcurves.c.  Summary files are
#                               generated through extract_goodpoints.csh
# Jun 14, 2008: Edward J. Los - Add rms error bar plot.
# Aug 22, 2008: Edward J. Los - Add ERROR_BAR_FACTOR, the factor by which magcal_local_rms and magcal_iso_rms has been
#                               reduced so that the clipped median RMS is equal to the clipped zero-based RMS of 
#                               magcal_local_rms for good stars.
# May  4, 2008: Edward J. Los - Change name to plot_lightcurves
# Jun 19, 2009: Edward J. Los - Add GSC catalog information to the lightcurves
# Jun 22, 2009: Edward J. Los - Add JDBEGIN and JDEND to correct skew between upper and lower graph labels.
#
#               NOTE: format of the .db and .col files is identical and defined in extract_lightcurves.c
#

set scripts = $DASCH_SCRIPTS
foreach lc ( lc_*.db )
 set name = `echo $lc:r | sed 's/+/\\+/g'`
 set colfile = ${lc:r}.col
 set txtfile = ${lc:r}.txt


 column -b -i $lc > $colfile

 set raw_median_local  = `cat $txtfile | gawk 'NR==2{print  $1}'`
 set raw_rms_local     = `cat $txtfile | gawk 'NR==2{print  $2}'`
 set clip_median_local = `cat $txtfile | gawk 'NR==2{print  $3}'`
 set clip_rms_local    = `cat $txtfile | gawk 'NR==2{print  $4}'`
 set npoints           = `cat $txtfile | gawk 'NR==2{print  $5}'`
 set ngood_points      = `cat $txtfile | gawk 'NR==2{print  $6}'`
 set yrbegin           = `cat $txtfile | gawk 'NR==2{print  $7}'`
 set yrend             = `cat $txtfile | gawk 'NR==2{print  $8}'`
 set error_bar_factor  = `cat $txtfile | gawk 'NR==2{print  $9}'`
 set jdbegin           = `cat $txtfile | gawk 'NR==2{print $10}'`
 set jdend             = `cat $txtfile | gawk 'NR==2{print $11}'`


 set REF     = `cat $txtfile | gawk 'NR==2{print $12}'`
 set ra      = `cat $txtfile | gawk 'NR==2{print $13}'`
 set dec     = `cat $txtfile | gawk 'NR==2{print $14}'`
 set Stdmag  = `cat $txtfile | gawk 'NR==2{print $15}'`
 set color   = `cat $txtfile | gawk 'NR==2{print $16}'`
 set VFlag   = `cat $txtfile | gawk 'NR==2{print $17}'`
 set MAGFlag = `cat $txtfile | gawk 'NR==2{print $18}'`
 set RaPM    = `cat $txtfile | gawk 'NR==2{print $19}'`
 set DecPM   = `cat $txtfile | gawk 'NR==2{print $20}'`
 set minMAG  = `cat $txtfile | gawk 'NR==2{print $21}'`
 set maxMAG  = `cat $txtfile | gawk 'NR==2{print $22}'`

 

 cat ${scripts}/plot_lc_template_flags_local.wip | sed "s/OUTPUT/${lc:r}/g" | sed "s/FILE/${lc:r}.col/g" | sed "s/CLIP_MEDIAN_LOCAL/${clip_median_local}/g" | sed "s/RAW_MEDIAN_LOCAL/${raw_median_local}/g" | sed "s/NPOINTS/${npoints}/g" | sed "s/NGPOINTS/${ngood_points}/g" | sed "s/CLIP_RMS_LOCAL/${clip_rms_local}/g" | sed "s/RAW_RMS_LOCAL/${raw_rms_local}/g" | sed "s/YRBEGIN/${yrbegin}/g" | sed "s/YREND/${yrend}/g" | sed "s/ERROR_BAR_FACTOR/${error_bar_factor}/g" | sed "s/JDBEGIN/${jdbegin}/g" | sed "s/JDEND/${jdend}/g" | sed "s/REFNO/${REF}/g" | sed "s/RASC/${ra}/g" | sed "s/DECLINATION/${dec}/g" | sed "s/STDMAG/${Stdmag}/g" | sed "s/COLOR/${color}/g" | sed "s/VFLAG/${VFlag}/g" | sed "s/MAGFLAG/${MAGFlag}/g" | sed "s/RAPM/${RaPM}/g" | sed "s/DECPM/${DecPM}/g"  >! ${lc:r}.wip

 cat ${scripts}/plot_lc_template_flags_rms.wip   | sed "s/OUTPUT/${lc:r}_rms/g" | sed "s/FILE/${lc:r}.col/g" | sed "s/CLIP_MEDIAN_LOCAL/${clip_median_local}/g" | sed "s/RAW_MEDIAN_LOCAL/${raw_median_local}/g" | sed "s/NPOINTS/${npoints}/g" | sed "s/NGPOINTS/${ngood_points}/g" | sed "s/CLIP_RMS_LOCAL/${clip_rms_local}/g" | sed "s/RAW_RMS_LOCAL/${raw_rms_local}/g" | sed "s/YRBEGIN/${yrbegin}/g" | sed "s/YREND/${yrend}/g" | sed "s/ERROR_BAR_FACTOR/${error_bar_factor}/g" | sed "s/JDBEGIN/${jdbegin}/g" | sed "s/JDEND/${jdend}/g" | sed "s/REFNO/${REF}/g" | sed "s/RASC/${ra}/g" | sed "s/DECLINATION/${dec}/g" | sed "s/STDMAG/${Stdmag}/g" | sed "s/COLOR/${color}/g" | sed "s/VFLAG/${VFlag}/g" | sed "s/MAGFLAG/${MAGFlag}/g" | sed "s/RAPM/${RaPM}/g" | sed "s/DECPM/${DecPM}/g"  | sed "s/MINMAG/${minMAG}/g"  | sed "s/MAXMAG/${maxMAG}/g" >! ${lc:r}_rms.wip

# If WIP is installed on the system, you can uncomment the following line to produce plots.

wip ${lc:r}.wip
wip ${lc:r}_rms.wip



end

# These commands gather all the lightcurve plots into a single multi-page postscript file.
#/home/silas/pscombine.csh lc*.ps 
#mv newfile.ps  lightcurves.ps


