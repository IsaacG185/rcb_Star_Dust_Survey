#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Look for close multiple exposures. No-op for solnum > 0.
#
# Inputs:
#
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww.db
# - (doesn't need FITS, I think?)
#
# Outputs:
#
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww_background.db
#
# Database updates:
#
# - Set multiplicity flags on `FitWCS` in `scanner.mosaics`
#
# Invokes:
#
# - getlocation
# - getdirectory
# - search_close
#
# Variables used:
#
# - $DASCH_MATCH

if [ $# -ne 2 ] ; then
  echo "usage: $0 <solutionNumber> <list>"
  exit 1
fi

solutionNumber=$1
imagelist=$2

binning=16
doplots=0
dobackground=1
plotstr=""
calibrationSquare=" "
analysis_threshold=2.2

if [ $doplots -eq 1 ] ; then
  plotstr=" -q -p "
fi

if [ $solutionNumber != 0 ] ; then
  # We are done. This code is independent of the solution number
  exit
fi

echo "Processing $imagelist"

for plate in $(cat $imagelist) ; do
  backgroundname=${DASCH_MATCH}/${plate}_background.db
  rm -f $backgroundname

  backgroundtmp=${DASCH_MATCH}/${plate}_background.tmp
  rm -f $backgroundtmp

  if [ $dobackground -eq 1 ] ; then
    mosaicDirectory=$(getdirectory $plate)
    calibrationSquare=" -b 1 -d $backgroundtmp ${mosaicDirectory}/background.reg"
  fi

  series=$(getlocation -s $plate |awk '{print $1}')

  outputname=${DASCH_MATCH}/${plate}_close.db
  rm -f $outputname

  gmt0name=${DASCH_MATCH}/${plate}_close_gmt0.txt
  gmt1name=${DASCH_MATCH}/${plate}_close_gmt1.txt
  gmt2name=${DASCH_MATCH}/${plate}_close_gmt2.txt
  gmt3name=${DASCH_MATCH}/${plate}_close_gmt3.txt
  gmt4name=${DASCH_MATCH}/${plate}_close_gmt4.txt
  emulsionname=${DASCH_MATCH}/${plate}_close_plot.db

  rm -f $gmt0name
  rm -f $gmt1name
  rm -f $gmt2name
  rm -f $gmt3name
  rm -f $gmt4name
  rm -f $emulsionname

  if [ -e ${DASCH_MATCH}/${plate}_tnx.db ] ; then
    inputname=${DASCH_MATCH}/${plate}_tnx.db
  elif [ -e ${DASCH_MATCH}/${plate}.db ] ; then
    inputname=${DASCH_MATCH}/${plate}.db
  else
    echo "ERROR: search_close found no sextractor file for $plate"
    continue
  fi

  mosaicsize=$(getlocation -a $plate)
  width=$(echo "$mosaicsize" |awk '{print $1}')
  height=$(echo "$mosaicsize" |awk '{print $2}')
  platescale=$(echo "$mosaicsize" |awk '{print $3}')
  if [ -z "$platescale" ] ; then
    # echo "ERROR getlocation failed with args -a $image"
    continue
  fi

  if [ -e $inputname ] ; then
    set -x

    search_close $plotstr \
      -t $analysis_threshold \
      -u $calibrationSquare \
      -w $width \
      -s $platescale \
      -h $height \
      -r ${plate} \
      -i $inputname \
      -o $outputname

    set +x

    if [ -e $backgroundtmp ] ; then
        sorttable -n NUMBER <$backgroundtmp >$backgroundname
        rm $backgroundtmp
    fi
  else
      echo "ERROR: Cannot find input file $inputname"
  fi
done

#  Dec 20, 2011 Edward J. Los - Adapted from filter_wedge.csh
#  Jun  1, 2012 Edward J. Los - Add background search of calibration squares
#  Jul  3, 2012 Edward J. Los - Integrate high background object support with rest of pipeline
#  Jul 23, 2012 Edward J. Los - Skip this for solutionNumber != 0;
