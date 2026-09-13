#!/bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Register raw mosaics automatically using Astrometry.Net.
#
# TO CHECK: we only ever use this with solnum = 0, I think;
# `find_astrometry2.sh` handles solnum > 0.
#
# Inputs:
#
# - $RAID/ExposureData/Mosaics/${series}/${id}_${mosnum}/${series}${id}_${mosnum}_16.fit
#
# Outputs:
#
# - $DASCH_ASTROMETRY/${series}${id}_${mosnum}_16_s0.corr.fits
# - $DASCH_ASTROMETRY/${series}${id}_${mosnum}_16.solve
# - $DASCH_ASTROMETRY/${series}${id}_${mosnum}_16.wcs
# - $DASCH_ASTROMETRY/${series}${id}_${mosnum}_01.db
#
# No database updates.
#
# Invokes:
#
# - db2xyls
# - getdirectory
# - getlocation
# - getseries
#
# External dependencies:
#
# - Astrometry.net (solve-field)
# - Starbase (column)
# - listhead
#
# Variables used:
#
# - $DASCH_ASTROMETRY
# - $DASCH_SCRIPTS

imagelist=$1

binning=16
debugMode=0
fractionFLAG="-f 0.01"
kronFLAG=""
configSuffix=1

echo "Processing $imagelist"

for plate in $(cat $imagelist) ; do
  solvedFlag=0
  thumbnail=$(echo "$plate.fit" |sed -e 's/_01\./_16\./g' -e 's/_01w/_16w/g' -e 's/_01r/_16r/g')
  directory=$(getdirectory $plate)
  filename=${directory}/${thumbnail}
  rootname=$(echo "$thumbnail" |sed 's/\.fit$//')
  scaletext=$(getseries $plate -f 1.1 -b $binning)
  minscale=$(echo $scaletext |awk '{print $1}')
  maxscale=$(echo $scaletext |awk '{print $3}')

  if [ -e ${DASCH_ASTROMETRY}/${plate}.db ] ; then
    echo "ERRORD: find_astrometry ${DASCH_ASTROMETRY}/${plate}.db already exists"
    continue
  fi

  rm -f ${DASCH_ASTROMETRY}/fail_${plate}.txt

  if [ -e $filename ] ; then
    xylsfile=$(echo "$filename" |sed -e 's/\.fit$/\.xyls/g')
    xylsfile1=${xylsfile}first
    xylsfile2=${xylsfile}second
    resultfile=$(echo "$filename" |sed -e 's/\.fit$/\.solve/g')

    if [ $debugMode == 1 ] ; then
      debugfile=$(echo "$filename" |sed -e 's/\.fit$/\.debug/g')
      debugflag="-d $debugfile"
    else
      debugflag=""
    fi

    width=$(listhead $filename |grep "NAXIS1" |awk '{print $3}')
    height=$(listhead $filename |grep "NAXIS2" |awk '{print $3}')

    set -x

    sex \
      -c $DASCH_SCRIPTS/Sextractor/astrometry${configSuffix}.config \
      -CATALOG_NAME $xylsfile \
      -PARAMETERS_NAME $DASCH_SCRIPTS/Sextractor/astrometry.param \
      $filename \
      -FILTER_NAME $DASCH_SCRIPTS/Sextractor/default.conv

    sextotable <$xylsfile |sed -e 's/ALPHA_J2000/ra/g' -e 's/DELTA_J2000/dec/g' >$xylsfile1

    set +x

    # cases 1-5 are old version of astrometry.net
    # cases 6 7 8 are for astrometry.net-0.38/

    for casenumber in 6 7 8 1 2 ; do
      if [ $casenumber == 1 ] ; then
        filterFLAG="-e"
        ellipticityFLAG=""
        maxdepth="--depth 1-120"
        maxstars=240
        cpulimit=""
      elif [ $casenumber == 2 ] ; then
        filterFLAG="-e"
        ellipticityFLAG="-l 0.3"
        maxdepth="--depth 1-120"
        maxstars=240
        cpulimit=""
      elif [ $casenumber == 6 ] ; then
        filterFLAG="-e"
        ellipticityFLAG=""
        maxdepth=""
        maxstars=480
        cpulimit="--cpulimit 300"
      elif [ $casenumber == 7 ] ; then
        filterFLAG="-e"
        ellipticityFLAG="-l 0.3"
        maxdepth=""
        maxstars=480
        cpulimit="--cpulimit 300"
      elif [ $casenumber == 8 ] ; then
        filterFLAG=""
        ellipticityFLAG=""
        maxdepth=""
        maxstars=480
        cpulimit="--cpulimit 300"
      fi

      set -x

      db2xyls \
        -i $xylsfile1 \
        -o $xylsfile2 \
        -w $width \
        -h $height \
        -n $maxstars \
        $fractionFLAG \
        $ellipticityFLAG \
        $kronFLAG \
        $filterFLAG \
        $debugflag

      set +x

      cp $xylsfile2 $xylsfile

      if [ $debugMode == 1 ] ; then
        column -a -b <$debugfile >${DASCH_ASTROMETRY}/${rootname}.debug
        rm $debugfile
      fi

      date
      set -x

      solve-field \
        $xylsfile \
        --no-tweak \
        --scale-units arcsecperpix \
        --scale-low $minscale \
        --scale-high $maxscale \
        --x-column X_IMAGE \
        --y-column Y_IMAGE \
        --sort-column MAG_ISO \
        --sort-ascending \
        --overwrite \
        --corr ${DASCH_ASTROMETRY}/${rootname}_s0.corr.fits \
        --width $width \
        --height $height \
        $cpulimit \
        $maxdepth \
        --no-plots \
        >&$resultfile

      set +x
      date

      centerra=$(cat $resultfile |grep "Field center: (RA,Dec)" |awk '{print $5}' |sed -e 's/[(),]//g')
      centerdec=$(cat $resultfile |grep "Field center: (RA,Dec)" |awk '{print $6}' |sed -e 's/[(),]//g')

      if [ -z "$centerdec" ] ; then
        cp $resultfile ${DASCH_ASTROMETRY}/fail_${plate}.txt
      else
        solvedFlag=1
        break
      fi
    done

    if [ $solvedFlag == 0 ] ; then
      echo "No result for $plate"
      exit 30  # signal a "managed failure" to the job management system
    fi

    rm -f ${DASCH_ASTROMETRY}/fail_${plate}.txt
    solutiontext=$(grep "Solution found" $resultfile)
    echo "$solutiontext for case $casenumber"

    listhead ${directory}/${rootname}.wcs >${rootname}.gethead
    naxis1=$(cat ${rootname}.gethead |grep 'IMAGEW' |awk '{print $3}')
    naxis2=$(cat ${rootname}.gethead |grep 'IMAGEH' |awk '{print $3}')
    ctype1=$(cat ${rootname}.gethead |grep 'CTYPE1' |awk '{print $3}' |sed 's/\-SIP//g' |sed "s/'//g")
    ctype2=$(cat ${rootname}.gethead |grep 'CTYPE2' |awk '{print $3}' |sed 's/\-SIP//g' |sed "s/'//g")
    crval1=$(cat ${rootname}.gethead |grep 'CRVAL1' |awk '{print $3}')
    crval2=$(cat ${rootname}.gethead |grep 'CRVAL2' |awk '{print $3}')
    crpix1=$(cat ${rootname}.gethead |grep 'CRPIX1' |awk '{print $3}')
    crpix2=$(cat ${rootname}.gethead |grep 'CRPIX2' |awk '{print $3}')
    cd1_1=$(cat ${rootname}.gethead |grep 'CD1_1' |awk '{print $3}')
    cd1_2=$(cat ${rootname}.gethead |grep 'CD1_2' |awk '{print $3}')
    cd2_1=$(cat ${rootname}.gethead |grep 'CD2_1' |awk '{print $3}')
    cd2_2=$(cat ${rootname}.gethead |grep 'CD2_2' |awk '{print $3}')
    scalex=0 # newer Astrometry.Net does not print this info
    scaley=0

    echo "plate ra dec scalex scaley binning naxis1 naxis2 ctype1 ctype2 crval1 crval2 crpix1 crpix2 cd1_1 cd1_2 cd2_1 cd2_2" >${DASCH_ASTROMETRY}/${plate}.tmp
    echo "----- -- --- ------ ------ ------- ------ ------ ------ ------ ------ ------ ------ ------ ----- ----- ----- -----" >>${DASCH_ASTROMETRY}/${plate}.tmp
    echo "$plate $centerra $centerdec $scalex $scaley $binning $naxis1 $naxis2 $ctype1 $ctype2 $crval1 $crval2 $crpix1 $crpix2 $cd1_1 $cd1_2 $cd2_1 $cd2_2" >>${DASCH_ASTROMETRY}/${plate}.tmp
    cat ${DASCH_ASTROMETRY}/${plate}.tmp |awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' >${DASCH_ASTROMETRY}/${plate}.db
    rm ${DASCH_ASTROMETRY}/${plate}.tmp
    rm ${rootname}.gethead
    echo "Solved: ra $centerra dec $centerdec for $plate"

    mv $resultfile ${DASCH_ASTROMETRY}/${rootname}.solve
    echo $casenumber >${DASCH_ASTROMETRY}/${rootname}.casenumber.txt

    if [ -f ${rootname}.wcs ] ; then
      mv ${rootname}.wcs ${DASCH_ASTROMETRY}/${rootname}.wcs
    fi

    if [ -f ${directory}/${rootname}.wcs ] ; then
      mv ${directory}/${rootname}.wcs ${DASCH_ASTROMETRY}/${rootname}.wcs
    fi

    rm -f ${rootname}.axy
    rm -f ${rootname}.match
    rm -f ${rootname}.rdls
    rm -f ${rootname}.solved
    rm -f ${directory}/${rootname}.axy
    rm -f ${directory}/${rootname}.match
    rm -f ${directory}/${rootname}.rdls
    rm -f ${directory}/${rootname}.solved
    rm -f ${directory}/${rootname}-indx.xyls
    rm -f ${directory}/${rootname}.corr
    rm -f ${rootname}-indx.xyls
    rm -f $resultfile
    rm -f $xylsfile
    rm -f $xylsfile1
    rm -f $xylsfile2
  else
    echo "ERROR: $filename not found"
  fi
done

#  May  2, 2008 Edward J. Los - initial version
#  May 20, 2008 Edward J. Los - handle multiple cases
#  May 26, 2008 Edward J. Los - correct solvedFlag handling
#  Jun 27, 2008 Edward J. Los - replace gethead with listhead
#  Apr 21, 2008 Edward J. Los - remove fail* file if one of the solutions succeeds
#  Mar 29, 2010 Edward J. Los - do not over-write an existing solution
#  May 17, 2011 Edward J. Los - Modify for astrometry.net-0.38   Try using the internal timer
#  Jul  5, 2011 Edward J. Los - Extend cpulimit to 42 minutes except for cases 6,7, and 8
#  Apr 15, 2014 Edward J. Los - Mark miscellaneous errors deferred with ERRORD
#  Apr  8, 2016 Edward J. Los - Set the protection of the thumbnail to read only
#  Sep 12, 2016 Edward J. Los - Back out the protection change
