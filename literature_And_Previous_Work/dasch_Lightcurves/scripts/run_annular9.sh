#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Run the `annular9.m` script.
#
# It is possible that not all spatial bins will calibrate successfully. That is OK.
#
# Inputs:
#
# - $DASCH_CATALOGBIN/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_a{1..9}.colorterm.txt
# - $DASCH_CATALOGBIN/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_a{1..9}.db
#
# Outputs:
#
# - $DASCH_BINOUTPUT/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_a{1..9}.grid
# - $DASCH_BINOUTPUT/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_a{1..9}.out
# - $DASCH_BINOUTPUT/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_a{1..9}.para
# - $DASCH_BINOUTPUT/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.annular9_marker.txt
# - ~DELETES~ $DASCH_CATALOGBIN/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.colorterm_marker.txt
#
# Database updates: none.
#
# Invokes:
#
# - annular9.m
#
# External dependencies:
#
# - Octave
#
# Variables used:
#
# - $DASCH_BINOUTPUT
# - $DASCH_CATALOGBIN
# - $DASCH_NUMBINS (needed by annular9.m)
# - $DASCH_PLOT (needed by annular9.m)
# - $DASCH_SCRATCH
# - $DASCH_SCRIPTS

solutionNumber=$1
imagelist=$2

if [ $# -eq 2 ] ; then
  qualifier=""
  magnitudefile=""
  qualifier=""
  qualifier2=""
elif [ $# -eq 3 ]; then
  qualifier=$3
  qualifier2="_$3"
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

for plate in $(cat $imagelist) ; do
  catalog=$(ls ${DASCH_CATALOGBIN}/${plate}${solutionString}${qualifier2}_a?.db)
  markerfile=${DASCH_BINOUTPUT}/${plate}${solutionString}${qualifier2}.annular9_marker.txt
  platelist=${DASCH_SCRATCH}/${plate}.xxxx
  echo "$plate" >$platelist

  if [ -n "$catalog" ] ; then
    rm -f $markerfile
    octave $DASCH_SCRIPTS/annular9.m $solutionNumber $platelist $qualifier

    if [ ! -e $markerfile ] ; then
      echo "ERROR no marker file $markerfile"
    fi
  fi

  rm -f $platelist
done

#  Apr 13, 2011 Edward J. Los - Initial version
#  Apr 15, 2014 Edward J. Los - Mark miscellaneous errors deferred with ERRORD
