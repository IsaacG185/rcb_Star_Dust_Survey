#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Apply calibrations to catalogs.
#
# Inputs:
#
# - ${DASCH_INGEST}/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.local.db - Best quality locally calibrated stars
# - ${DASCH_INGEST}/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_dmagcor.db - used to apply local calibration
# - ${DASCH_INGEST}/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.out.spatial_bins.db - defines the limits
# - ${DASCH_MATCH}/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_tnx_u.db - All objects matched with the GSC catalog
# - ${DASCH_MATCH}/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.db - All objects found by sextractor.
# - ${DASCH_MATCH}/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_defect.db - Defect filter output.
# - ${DASCH_CATALOGBIN}/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_a${bin}.grid - used to convert MAG_ISO to magcal_iso
#
# Outputs:
#
# - ${DASCH_INGEST}/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_allobjects.db
# - ~UPDATE~ ${DASCH_INGEST}/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_dmagcor.db
# - ${DASCH_INGEST}/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_dmagcor_copy.db (backup of original dmagcor.db file)
#
# Database updates: none.
#
# Invokes:
#
# - getlocation
# - recover_points
#
# Variables used:
#
# - $DASCH_BINOUTPUT
# - $DASCH_CATALOGALL
# - $DASCH_INGEST
# - $DASCH_MATCH
# - $DASCH_NUMBINS

totalcounter=0
goodcounter=0
missingcounter=0

solutionNumber=$1
imagelist=$2

# This script used to require a `catrms=$MAGS` argument, but with the current
# "RMS Algorithm E" that parameter is not used. If we end up needing it again in
# the future, we should have this script automatically set the parameter based
# on the qualifier, I think.

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
  echo "  qualifier is a string appended to output filenames"
  exit 1
fi

if [ $solutionNumber -eq 0 ] ; then
  solutionString=""
else
  solutionString="_s$solutionNumber"
fi

echo "$(date) Processing ${imagelist}${qualifier}"

for name in $(cat $imagelist) ; do
  mosaicsize=$(getlocation -a $name -e $solutionNumber)
  width=$(echo "$mosaicsize" |awk '{print $1}')
  height=$(echo "$mosaicsize" |awk '{print $2}')
  extinction_file=$(echo ${DASCH_CATALOGALL}/${name}${solutionString}${qualifier}_extinction.db)
  script_output=$(getlocation $name -e $solutionNumber)
  plate_scale=$(echo "$script_output" |awk '{print $4}')

  if [ -z "$plate_scale" ] ; then
    # echo "ERROR getlocation failed with args  $name -e $solutionNumber)
    continue
  fi

  localfilename="${DASCH_INGEST}/${name}${solutionString}${qualifier}.out.local.db"
  dmagcorfilename="${DASCH_INGEST}/${name}${solutionString}${qualifier}_dmagcor.db"
  dmagcorcopy="${DASCH_INGEST}/${name}${solutionString}${qualifier}_dmagcor_copy.db"

  if [ -e $dmagcorfilename -a ! -e $dmagcorcopy ] ; then
    echo "cp $dmagcorfilename $dmagcorcopy"
    cp $dmagcorfilename $dmagcorcopy
  fi

  if [ "$qualifier" == "_kepler" -o "$qualifier" == "_apass" -o "$qualifier" == "_gaia" -o "$qualifier" == "_atlas" -o "$qualifier" == "_experimental" ] ; then
    matchfilename="${DASCH_MATCH}/match_${name}${solutionString}_tnx${qualifier}_u.db"
    nomtchfilename="${DASCH_MATCH}/match_${name}${solutionString}_tnx${qualifier}_n.db"
    blendfilename="${DASCH_MATCH}/match_${name}${solutionString}_tnx${qualifier}_b.db"
    dradfilename="${DASCH_MATCH}/${name}${solutionString}${qualifier}_drad.db"
  else
    matchfilename="${DASCH_MATCH}/match_${name}${solutionString}_tnx_u.db"
    nomtchfilename="${DASCH_MATCH}/match_${name}${solutionString}_tnx_n.db"
    blendfilename="${DASCH_MATCH}/match_${name}${solutionString}_tnx_b.db"
    dradfilename="${DASCH_MATCH}/${name}${solutionString}_drad.db"
  fi

  sextractorfilename="${DASCH_MATCH}/${name}${solutionString}_tnx.db"
  outputfilename="${DASCH_INGEST}/${name}${solutionString}${qualifier}_allobjects.db"

  rm -f $outputfilename

  if [ $solutionNumber -eq 0 ] ; then
    rm -f ${DASCH_INGEST}/${name}${solutionString}_s0${qualifier}_allobjects.db
  fi

  ((totalcounter+=1))

  if [ -e $localfilename ] ; then
    echo "Processing ${name}${solutionString} at $(date)"
    plate=$(echo "${name}${solutionString}" |awk -F"_" '{print $1}')

    cd $DASCH_INGEST

    if [ ! -e $dradfilename  ] ; then
      echo "ERROR: drad file $dradfilename does not exist"
      continue
    fi

    if [ -e $matchfilename ] ; then
      if [ -e $sextractorfilename ] ; then
        ((goodcounter+=1))

        # Create lowess grid calibration files
        grid=$(ls ${DASCH_BINOUTPUT}/${name}${solutionString}${qualifier}_a?.grid)

        if [ -n "$grid" ] ; then
          for grid in ${DASCH_BINOUTPUT}/${name}${solutionString}${qualifier}_a?.grid ; do
              tmpgridfile="${DASCH_INGEST}/$(basename $grid .grid)_grid.tm"
              echo "maggrid isogrid griderr flaggrid" >$tmpgridfile
              echo "------- ------- ------- --------" >>$tmpgridfile
              cat $grid >>$tmpgridfile
              awk '{OFS="\t"}{print $1,$2,$3,$4}' <$tmpgridfile >${tmpgridfile}p
          done

          # Now process the results
          set -x
          recover_points \
            -a $solutionNumber \
            $qualcmd \
            -c catrms=0 \
            -p ${name}${solutionString} \
            -e $extinction_file \
            -w $width \
            -h $height \
            -l $localfilename \
            -m $matchfilename \
            -n $nomtchfilename \
            -d $blendfilename \
            -r $plate_scale \
            -s $sextractorfilename \
            -o $outputfilename \
            -b $DASCH_NUMBINS
          set +x

          for grid in ${DASCH_BINOUTPUT}/${name}${solutionString}${qualifier}_a?.grid ; do
            tmpgridfile="${DASCH_INGEST}/$(basename $grid .grid)_grid.tm"
            rm $tmpgridfile
            rm ${tmpgridfile}p
          done
        else
          echo "WARNING: No calibration files found for ${name}${solutionString}"
        fi
      else
        echo "ERROR: No sextractor file found: $sextractorfilename"
        ((missingcounter+=1))
      fi
    else
      ((missingcounter+=1))
      echo "ERROR: No match file found: $matchfilename"
    fi
  else
    echo "No file found: $localfilename"
    ((missingcounter+=1))
  fi
done

echo "Total Plates: $totalcounter Found: $goodcounter Missing: $missingcounter"
echo "$(date) Completed ${imagelist}${qualifier}"

#  Feb  1, 2008 Edward J. Los - Initial version
#  Feb 17, 2008 Edward J. Los - Convert missing file error to a warning
#  Feb 18, 2008 Edward J. Los - Currect pass 2 by adding a qualifier
#                             - Correct Julian Date
#                             - Add correction RMS for the catalog error
#  Feb 25, 2008 Edward J. Los - Handle single bin case by rejecting any points
#                               outside of the spatial limit of the calibrated data.
#  Mar  5, 2008 Edward J. Los - Add the number of bins as a parameter
#                               Provide pointer to an extinction file
#  Mar 14, 2008 Edward J. Los - Remove any possible stale output files from previous
#                               executions of the pipeline.
#  Apr  1, 2008 Edward J. Los - Move match file to $DASCH_MATCH
#  Aug 11, 2008 Edward J. Los - Add magnitude adjustment for blended stars
#                               Add the plate scale parameter
#  Aug 20, 2009 Edward J. Los - Add Kepler Input Catalog support
#  Sep 29, 2009 Edward J. Los - add multiple exposure support
#  Nov 16, 2009 Edward J. Los - Move Julian Date calculation to within recover_points
#  Nov 25, 2009 Edward J. Los - Check for drad file before calling recover_points
#  May 23, 2010 Edward J. Los - The kepler pass needs only the kepler drad file
#  Mar 23, 2011 Edward J. Los - Add apass catalog support
#  Jul 30, 2012 Edward J. Los - Add experimental catalog support
#  Aug 16, 2012 Edward J. Los - Make a copy of dmagcor.db before recover_points modifies it.
#  Jun 27, 2010 Edward J. Los - Prevent use of a stale multiple exposure source list.
#  May 28, 2018 Edward J. Los - Add gaia catalog support
#  Oct 28, 2018 Edward J. Los - Add atlas refcat2 catalog support
