#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Perform magnitude-dependent local correction after run_local_calibration.sh is
# executed.
#
# This calibration might not succeed. As far as the pipeline is concerned, that
# is undesirable, but OK.
#
# Inputs:
#
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.local.db
#
# Outputs:
#
# - $DASCH_INGEST/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_magdepcalibrate.db
#
# Database updates: none.
#
# Invokes:
#
# - magdepcalibrate
#
# Variables used:
#
# - $DASCH_INGEST

solutionNumber=$1
imagelist=$2

binning=16
doplots=0
plotstr=""
useRaDec=0

if [ $# -eq 2 ] ; then
  qualifier=""
  qualcmd=""
elif [ $# -eq 3 ] ; then
  qualifier="_$3"
  qualcmd="-q $qualifier"
else
  echo "usage: $0 <solutionNumber> <list> [qualifier]"
  exit 1
fi

if [ $solutionNumber -eq 0 ] ; then
  solutionString=""
else
  solutionString="_s$solutionNumber"
fi

echo "Processing $imagelist"

for plate in $(cat $imagelist) ; do
  inputname=${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.local.db
  calibratename=${DASCH_INGEST}/${plate}${solutionString}${qualifier}_magdepcalibrate.db
  plot1name=${DASCH_INGEST}/${plate}${solutionString}${qualifier}_magdep1.ps
  plot2name=${DASCH_INGEST}/${plate}${solutionString}${qualifier}_magdep2.ps

  rm -f $calibratename
  rm -f $plot1name
  rm -f $plot2name

  if [ -e $inputname ] ; then
    set -x
    magdepcalibrate \
      $qualcmd \
      -r ${plate}${solutionString} \
      -i $inputname \
      -c $calibratename
    set +x

    if [ ! -f $calibratename ] ; then
      echo "WARNING: calibration failed for $plate; processing will continue"
    fi
  else
    echo "ERRORD: Cannot find input file $inputname"
  fi
done

#  Oct  4, 2011 Edward J. Los - modelled after Sumin Tang's
#                               bright_dubious_variales.m  of Sep 24, 2011
#  (Matlab) ported to Octave as calibratebright.m Apr 15, 2014 Edward J. Los -
#  Mark miscellaneous errors deferred with ERRORD
