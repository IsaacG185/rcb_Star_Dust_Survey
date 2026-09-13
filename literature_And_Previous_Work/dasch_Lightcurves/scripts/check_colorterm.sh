#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# This file looks for the colorterm marker
# `${DASCH_CATALOGBIN}/${plate}${solutionString}.colorterm_marker.txt` and sets
# or clears the ColortermCrash bit in the FitWCS column of the mosaics table of
# the scanner database.
#
# Inputs:
#
# - $DASCH_CATALOGBIN/${series}${id}_${mosnum}_01${rot}ww${sol_tag}.colorterm_marker.txt
#
# Outputs:
#
# - None
#
# Database updates:
#
# - Potentially updates `FitWCS` in `scanner.mosaics` as per above.
#
# Invokes:
#
# - setFitWCS

if [ $# -eq 2 ] ; then
  qualifier=
elif [ $# -eq 3 ] ; then
  qualifier="$3"
else
  echo "usage: $0 <solutionNumber> <list> [qualifier]"
  echo "  if qualifier is non-empty, this command is a no-op"
  exit 1
fi

solutionNumber=$1
imagelist=$2

if [ $solutionNumber -eq 0 ] ; then
  solutionString=""
else
  solutionString="_s$solutionNumber"
fi

if [ -n "$qualifier" ] ; then
  echo "check_colorterm is a no-op for qualifier $qualifier"
  exit 0
fi

echo "Processing $imagelist"

for plate in $(cat $imagelist) ; do
  inputfile=${DASCH_CATALOGALL}/${plate}${solutionString}.db
  markerfile=${DASCH_CATALOGBIN}/${plate}${solutionString}.colorterm_marker.txt

  if [ -e $inputfile ] ; then
    if [ -e $markerfile ] ; then
      setFitWCS -e $solutionNumber -b ColortermCrash -p ${markerfile} -c
    else
      setFitWCS -e $solutionNumber -b ColortermCrash -p ${markerfile} -s
    fi
  fi
done

#  Dec 16, 2009 Edward J. Los - Initial version
#  Dec 17, 2009 Edward J. Los - Set/clear flag only if colorterm input file exists
