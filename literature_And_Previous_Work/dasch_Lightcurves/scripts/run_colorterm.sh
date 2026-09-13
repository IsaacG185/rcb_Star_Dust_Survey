#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Run the `colorterm.m` script.
#
# Inputs:
#
# - $DASCH_CATALOGBIN/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_a{1..9}.db
#
# Outputs:
#
# - $DASCH_CATALOGBIN/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_a{1..9}.colorterm.txt
# - $DASCH_CATALOGBIN/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.colorterm_marker.txt
#
# Database updates: none.
#
# Invokes:
#
# - colorterm.m
#
# External dependencies:
#
# - Octave
#
# Variables used:
#
# - $DASCH_CATALOGALL
# - $DASCH_CATALOGBIN
# - $DASCH_PLOT (needed by colorterm.m)
# - $DASCH_SCRATCH
# - $DASCH_SCRIPTS

solutionNumber=$1
imagelist=$2

if [ $# -eq 2 ] ; then
  qualifier=""
  magnitudefile=""
  qualifier=""
  qualifier2=""
elif [ $# -eq 3 ] ; then
  qualifier=$3
  qualifier2="_$3"
else
  echo "usage: $0 <solutionNumber> <list> [qualifier]"
  echo "  solutionNumber is the multiple exposure iteration"
  echo "  list      is the list of mosaics to process"
  echo "  qualifier is a string append to output filenames"
  exit
fi

if [ $solutionNumber -eq 0 ] ; then
  solutionString=""
else
  solutionString="_s$solutionNumber"
fi

echo "Processing $imagelist"

for plate in $(cat $imagelist) ; do
  inputfile=${DASCH_CATALOGALL}/${plate}${solutionString}${qualifier2}.db
  markerfile=${DASCH_CATALOGBIN}/${plate}${solutionString}${qualifier2}.colorterm_marker.txt
  platelist=${DASCH_SCRATCH}/${plate}.xxxx
  echo "$plate" >$platelist

  if [ -e $inputfile ] ; then
    rm -f $markerfile
    octave $DASCH_SCRIPTS/colorterm.m $solutionNumber $platelist $qualifier

    if [ ! -e $markerfile ] ; then
      echo "ERROR no marker file $markerfile"
    fi
  fi

  rm -f $platelist
done

#  Sep 13, 2010 Edward J. Los - Initial version
#  Oct  4, 2010 Edward J. Los - Correct for lists with more than one plate
#  Apr  7, 2011 Edward J. Los - Correct handling of kepler and apass files
#  Apr 15, 2014 Edward J. Los - Mark miscellaneous errors deferred with ERRORD
