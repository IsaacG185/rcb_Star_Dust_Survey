#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# This small script runs the local re-calibration routine (map method) on a
# whole set of plates. Local recalibration adjusts each star's magnitude by the
# mean calibration residual of its neighbors. Starbase commands are used to
# properly format the input and output.
#
# Inputs:
#
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.db
#
# Outputs:
#
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.local.db
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.local.db.REF-i
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.local.db.dec-i
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.local.db.ra-b
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_dmagcor.db
#
# Database updates: none.
#
# Invokes:
#
# - contour_plot
#
# External dependencies:
#
# - Starbase (index)
#
# Variables used:
#
# - $DASCH_INGEST

solutionNumber=$1
imagelist=$2

if [ $# -eq 2 ] ; then
  qualifier=""
elif [ $# -eq 3 ] ; then
  qualifier="_$3"
else
  echo "usage: $0 <solutionNumber> <list> [qualifier]"
  echo "  solutionNumber is the multiple exposure iteration"
  echo "  list      is the list of mosaics to process"
  echo "  qualifier is a string append to output filenames"
  exit 1
fi

if [ $solutionNumber -eq 0 ] ; then
  solutionString=""
else
  solutionString="_s$solutionNumber"
fi

echo "Processing $imagelist"

for listentry in $(cat $imagelist) ; do
  mosaicsize=$(getlocation -a $listentry -e $solutionNumber)
  width=$(echo "$mosaicsize" |awk '{print $1}')
  height=$(echo "$mosaicsize" |awk '{print $2}')
  platescale=$(echo "$mosaicsize" |awk '{print $3}')

  if [ -z "$platescale" ] ; then
    # echo "ERROR getlocation failed with args  -a $listentry -e $solutionNumber)
    continue
  fi

  cd $DASCH_INGEST
  dbbase="${listentry}${solutionString}${qualifier}"
  db="${dbbase}.out.db"
  localdb="${dbbase}.out.local.db"
  rm -f ${localdb}

  if [ -e $db ] ; then
    date

    # old code 1 was here.

    rm -f ${dbbase}_dmagcor.grid
    rm -f ${dbbase}_dmagcor.db
    rm -f ${dbbase}_dmagcor_copy.db

    set -x
    contour_plot \
      -i $db \
      -lx X \
      -ly Y \
      -c 5 \
      -p ${dbbase}_dmag.ps/ps \
      -r 3 \
      -s 0.5 \
      -o ${localdb} \
      -g ${dbbase}_dmagcor.db \
      -a $platescale \
      -mx $width \
      -my $height

    # If contour_plot emits no rows, subsequent processing will
    # fail. In that case the local.db file has two lines: the
    # column names and the delimiter hyphens

    if [ ! -f "$localdb" ] ; then
      echo >&2 "error: no such contour_plot output file $localdb"
      exit 32
    fi

    if [ "$(wc -l "$localdb" |awk '{print $1}')" -lt 3 ] ; then
      echo >&2 "error: no data in contour_plot output file $localdb"
      exit 32
    fi

    # old code 2 was here

    # Sort and index the table for later use.
    index -mb ${localdb} ra
    index -mi ${localdb} dec
    index -mi ${localdb} REF
  else
    echo "Database $db does not exist"
  fi
done

# old code 1:
#
# First select the good stars only. These will be used to define the calibration map and will be corrected.
# We choose only those stars brighter than their local limiting magnitude. Other selections could aso be used.
# row -i $db '(magcal_iso<(limiting_mag-0.5))' >! goodstars${listentry}${solutionString}${qualifier}.tmp
# NOTE: The above 0.5 magnitude value should agree with the definition of MAX_LIMITING_MAG in computemag.c
#
# Here we select all the other stars that didn't make the cut, we want to keep them, but not try to recalibrate them.
# row -i $db '(magcal_iso>=limiting_mag)' >! otherstars${listentry}${solutionString}${qualifier}.tmp
#
# Write a list of the selected stars' RA, DEC and CALIBRATION_RESIDUAL
# NOTE: THE COORDS SHOULD IDEALLY BE IN X,Y FORM, RATHER THAN RA,DEC
# Here we are using the magnitude quantity dmag_iso, but it could equally well be dmag_aper, dmag_kron, dmag_area .... etc
# column -i goodstars${listentry}${solutionString}${qualifier}.tmp -b REF ra dec dmag_iso >! dmag${listentry}${solutionString}${qualifier}.tmp
# column -i goodstars${listentry}${solutionString}${qualifier}.tmp -b REF X_IMAGE Y_IMAGE dmag_iso magcal_iso >! dmag${listentry}${solutionString}${qualifier}.tmp


# old code 2:
#
# Old code to run the local re-calibration routine
#$DASCH_SCRIPTS/contour_plot.csh dmag${listentry}${solutionString}${qualifier}.tmp 50 50 RA Dec 5 ${dbbase}_dmag.ps/ps 3 0.5 ${dbbase}_dmagcor.col
#
# Format the output of the routine, for starbase.
# Column REF2 == REF and can be used for checking that no miss-matches occured in seperating and subsequently recombining the data.
#echo "magcor_local magcal_local_error npoints_local local_bin" |sed 's/ /\t/g' >! dmagcor${listentry}${solutionString}${qualifier}.tmp
#echo "------------ ------------------ ------------- ---------" |sed 's/ /\t/g' >> dmagcor${listentry}${solutionString}${qualifier}.tmp
#cat ${dbbase}_dmagcor.col |awk '{OFS="\t"}{print $5,$6,$7,$8}' >> dmagcor${listentry}${solutionString}${qualifier}.tmp
#
# Compute the recalibrated magnitudes. This step is necessary because the fortran calibration routine works in dmag units, for generality.
# magcor_local -- the local calibration offset to be applied (shorthand for: magnitude_correction_local)
# magcal_local -- Corrected magnitude (shorthand for: magnitude_calibrated_local)
# dmag_local   -- Residual difference between primary calibrator and locally calibrated magnitude
#                It should be clear that [(dmag_iso - dmag_local) / dmag_iso] is the fractional improvement attained.
# cal_local   -- This flag indicates whether local calibration was performed on this star. 1=YES, 0=NO
#
#paste goodstars${listentry}${solutionString}${qualifier}.tmp dmagcor${listentry}${solutionString}${qualifier}.tmp |column -a magcal_local dmag_local cal_local |compute 'magcal_local = magcal_iso - magcor_local; dmag_local = dmag_iso - magcor_local; cal_local=1' >! ${localdb}
#
# For the "other" stars that were too faint to be used (or otherwise unsuitable), initiate 'placeholder' local columns and assign flag values to them.
# Then combine these stars with the "good" stars so there is a single merged output table with proper values for all entries.
# The same columns must be present, and IN THE SAME ORDER, in both the "good" and "other" tables for correct merging!
# In particular for the "other" stars make a sensible choice for the placeholder values such as:  magcal_local=magcal_iso (or 99)
# HERE I CHOOSE TO SET THE UNCORRECTED LOCAL VALUES TO THE INPUT VALUES
# THE FLAG IN THIS CASE IS cal_local=0 (means point was NOT calibrated)
#
# Note: all other stars, including all of the unmatched sextractor images, will be handled in recover_points.c
#
# column -i otherstars${listentry}${solutionString}${qualifier}.tmp -a magcor_local magcal_local_error npoints_local magcal_local dmag_local cal_local |compute 'magcor_local=0; magcal_local_error=magcal_iso_rms; npoints_local=0; magcal_local=magcal_iso; dmag_local=dmag_iso; cal_local=0' |column -b >> ${localdb}


# Modified: Sept 12, 2007 Silas Laycock.
# Introduced a 'flag_local' column to indicate corrected/uncorrected stars
# Added several explanatory comments.
#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#  Dec 11, 2007 Edward J. Los - Reject stars if FLAGS > 256
#  Dec 20, 2007 Edward J. Los - Reject stars if FLAGS >= 65536
#  Dec 29, 2007 Edward J. Los - Support recalibrated GSC magnitudes.
#  Jan 22, 2008 Edward J. Los - Back out the FLAGS test
#                               Replace magcal_iso_error with magcal_iso_rms
#                               Replace contour_plot.csh with contour_plot.c so we can generate
#                                 a calibration grid file to apply to stars which did not make the cut.
#                               Use X_IMAGE and Y_IMAGE for binning.
#                               Skip all stars greater than the limiting magnitude.
#                               Add the local calibration bin number for cross-checking
#  Apr 19, 2008 Edward J. Los - Remove stale files
#  Jul  8, 2008 Edward J. Los - Use only stars that are 0.5 brighter than the limiting magnitude for
#                               local calibration.
#  Feb 16, 2009 Edward J. Los - Change .grid to .db
#  Sep 29, 2009 Edward J. Los - Add multiple solution support
#  Mar 11, 2011 Edward J. Los - Perform all filtering within the contour_plot routine
#  Jul 25, 2017 Edward J. Los - Remove the nx and ny parameters and make the number of smoothing bins dependent on mx and my
