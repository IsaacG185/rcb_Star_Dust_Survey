#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# This file performs a statistical analysis of Sextractor parameters and flags
# objects presumed to be not stars so that they can be removed later in the
# pipeline
#
# Inputs:
#
# - ${DASCH_MATCH}/match_${name}${sol_tag}_tnx_u.db - All objects matched with the GSC catalog
# - ${DASCH_MATCH}/${name}${sol_tag}_tnx.db - All objects found by sextractor.
#
# Outputs:
#
# - ~UPDATE~ ${DASCH_MATCH}/match_${name}${sol_tag}_tnx${qual}_u.db
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_defect.db
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_summary_defect.db
#
# Database updates: none.
#
# Invokes:
#
# - filter_defect
# - getlocation
#
# External dependencies:
#
# - Octave
# - Starbase (column, jointable, sorttable)
#
# Variables used:
#
# - $DASCH_MATCH
# - $DASCH_SCRIPTS

local=$(pwd)
totalcounter=0
goodcounter=0
missingcounter=0
doplots=0
plotstr=""
solutionNumber=$1
imagelist=$2

if [ $# -eq 2 ] ; then
  qualifier=""
  qualcmd=""
elif [ $# -eq 3 ] ; then
  qualifier="_$3"
  qualcmd="-q $qualifier"
else
  echo "usage: $0 <solutionNumber> <list> [qualifier]"
  echo "  solutionNumber is the multiple exposure iteration"
  echo "  list      is the list of mosaics to process"
  echo "  rms    catalog error. rms is in magnitudes"
  echo "  qualifier is a string appended to output filenames"
  exit 1
fi

if [ $solutionNumber -eq 0 ] ; then
    solutionString=""
else
    solutionString="_s$solutionNumber"
fi

if [ $doplots -eq 1 ] ; then
    plotstr=" -g "
fi

echo "$(date) Processing ${imagelist}${qualifier}"

for name in $(cat $imagelist) ; do
  mosaicsize=$(getlocation -a $name -e $solutionNumber)
  width=$(echo "$mosaicsize" |awk '{print $1}')
  height=$(echo "$mosaicsize" |awk '{print $2}')
  if [ -z "$height" ] ; then
    echo >&2 "ERROR getlocation failed with args -e $solutionNumber -a $name"
    exit 1
  fi

  matchfilename="${DASCH_MATCH}/match_${name}${solutionString}_tnx${qualifier}_u.db"
  sextractorfilename="${DASCH_MATCH}/${name}${solutionString}_tnx.db"
  outputfilename="${DASCH_MATCH}/${name}${solutionString}${qualifier}_defect.db"
  summaryfilename="${DASCH_MATCH}/${name}${solutionString}${qualifier}_summary_defect.db"

  if [ ! -e $matchfilename ] ; then
    echo >&2 "fatal error: no such input $matchfilename"
    exit 1
  fi

  rm -f $outputfilename
  rm -f $summaryfilename

  ((totalcounter=$totalcounter + 1))

  echo "Processing ${name}${solutionString} at $(date)"

  cd $DASCH_MATCH

  if [ -e $matchfilename ] ; then
    if [ -e $sextractorfilename ] ; then
        ((goodcounter=$goodcounter + 1))

        set -x
        filter_defect \
          $plotstr \
          $qualcmd \
          -p ${name}${solutionString} \
          -w $width \
          -h $height \
          -m $matchfilename \
          -s $sextractorfilename \
          -o $outputfilename \
          -r $summaryfilename
        set +x

        if [ -e $outputfilename ] ; then
          matchheader=$(column -h -i $matchfilename |head -n 1 |sed 's/\t/ /g')
          tempname2=${DASCH_MATCH}/${name}${solutionString}${qualifier}_defect.tmp2
          tempname3=${DASCH_MATCH}/${name}${solutionString}${qualifier}_defect.tmp3
          tempname4=${DASCH_MATCH}/${name}${solutionString}${qualifier}_defect.tmp4

          index -n -mb $matchfilename NUMBER
          column -i $outputfilename NUMBER AFLAGS >$tempname2
          index  -n -mb $tempname2 NUMBER
          jointable -j NUMBER -n $matchfilename $tempname2 >$tempname3
          column -a -i $tempname3 AFLAGS |compute 'AFLAGS = AFLAGS_2' |column $matchheader >$tempname4
          mv $tempname4 $matchfilename
          rm ${matchfilename}.NUMBER.b
          rm ${tempname2}.NUMBER.b
          rm -f $tempname2
          rm -f $tempname3
          rm -f $tempname4
        fi
    else
      echo "ERROR: No sextractor file found: $sextractorfilename"
      ((missingcounter=$missingcounter + 1))
    fi
  else
     ((missingcounter=$missingcounter + 1))
     echo "ERROR: No match file found: $matchfilename"
  fi
done

echo "Total Plates: $totalcounter Found: $goodcounter Missing: $missingcounter"
echo "$(date) Completed ${imagelist}${qualifier}"

#  Jun  3, 2008 Edward J.Los - Initial version
#  Jan 27, 2009 Edward J. Los - split FLAGS into AFLAGS and BFLAGS
#  Aug 20, 2009 Edward J. Los - Add Kepler catalog support
#  Sep 29, 2009 Edward J. Los - add multiple exposure support
#  Nov 20, 2010 Edward J. Los - Set directory to the match directory to avoid index cross-link issues.
