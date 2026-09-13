#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Register images automatically using Astrometry.Net, in the mode for
# multiple-exposure processing.
#
# Currently, a no-op if solnum = 0 or catalog != gsc232.
#
# Inputs:
#
# - $DASCH_ASTROMETRY/${series}${id}_${mosnum}_01ww_s${solnum}_none.db
#
# Outputs:
#
# - $DASCH_ASTROMETRY/${series}${id}_${mosnum}_01ww_s${solnum}.db
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

if [ -n "$qualifier" ] ; then
  echo "find_astrometry2 is a no-op for qualifier $qualifier"
  exit 0
fi

binning=16
((prevSolutionNumber=$solutionNumber - 1))

if [ $prevSolutionNumber -lt 0 ] ; then
  echo "find_astrometry2 ignoring solutionNumber $solutionNumber"
  #echo "ERROR: exposure must be 1 or greater"
  exit 0
fi

debugMode=0
fractionFLAG="-f 0.01"
kronFLAG=""

echo "Processing $imagelist"

for plate in $(cat ${imagelist}) ; do
    solvedFlag=0
    rootname=$(echo "$plate" |sed -e 's/_01\./_16\./g' -e 's/_01w/_16w/g' -e 's/_01r/_16r/g')
    mosid=$(echo "$rootname" |cut -d_ -f1,2)
    scaletext=$(getseries $plate -f 1.1 -b $binning)
    minscale=$(echo $scaletext |awk '{print $1}')
    maxscale=$(echo $scaletext |awk '{print $3}')

    if [ -e ${DASCH_ASTROMETRY}/${plate}_s${solutionNumber}.db  ] ; then
        echo "ERROR: find_astrometry2 ${DASCH_ASTROMETRY}/${plate}_s${solutionNumber}.db already exists"
        continue
    fi

    rm -f ${DASCH_ASTROMETRY}/fail_${plate}_s${solutionNumber}.txt

    xylsfile=${DASCH_ASTROMETRY}/${rootname}_s${solutionNumber}.xyls
    xylsfile1=${DASCH_ASTROMETRY}/${plate}_s${solutionNumber}_none.db

    if [ -e $xylsfile1 ] ; then
        echo "xylsfile1 is $xylsfile1"
        xylsfile2=${xylsfile}second
        resultfile=$(echo "$xylsfile" |sed -e 's/\.xyls$/\.solve/g')

        if [ $debugMode -eq 1 ] ; then
            debugfile=$(echo "$xylsfile" |sed -e 's/\.xyls$/\.debugdb/g')
            debugflag="-d $debugfile"
        else
            debugflag=""
        fi

        mosaicsize=$(getlocation -a $plate)
        width=$(echo "$mosaicsize" |awk '{print $1}')
        height=$(echo "$mosaicsize" |awk '{print $2}')
        ((width=$width / $binning))
        ((height=$height / $binning))

        # cases 1-5 are old version of astrometry.net
        # cases 6 7 8 are for astrometry.net-0.38/

        for casenumber in 6 7 8 1 2 ; do
            if [ $casenumber -eq 1 ] ; then
                filterFLAG="-e"
                ellipticityFLAG=""
                maxdepth="--depth 1-120"
                maxstars=240
                cpulimit=""
            fi
            if [ $casenumber -eq 2 ] ; then
                filterFLAG="-e"
                ellipticityFLAG="-l 0.3"
                maxdepth="--depth 1-120"
                maxstars=240
                cpulimit=""
            fi
            if [ $casenumber -eq 6 ] ; then
                filterFLAG="-e"
                ellipticityFLAG=""
                maxdepth=""
                maxstars=480
                cpulimit="--cpulimit 300"
            fi
            if [ $casenumber -eq 7 ] ; then
                filterFLAG="-e"
                ellipticityFLAG="-l 0.3"
                maxdepth=""
                maxstars=480
                cpulimit="--cpulimit 300"
            fi
            if [ $casenumber -eq 8 ] ; then
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
                $fractionFLAG $ellipticityFLAG $kronFLAG $filterFLAG $debugflag
            set +x

            cp $xylsfile2 $xylsfile

            if [ $debugMode -eq 1 ] ; then
                column -a -b <$debugfile >${DASCH_ASTROMETRY}/${rootname}_s${solutionNumber}.debug
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
                --corr ${DASCH_ASTROMETRY}/${mosid}_16_s${solutionNumber}.corr.fits \
                --width $width \
                --height $height \
                $cpulimit \
                $maxdepth \
                --no-plots >& $resultfile
            set +x
            date

            centerra=$(cat $resultfile |grep "Field center: (RA,Dec)" |awk '{print $5}' |sed -e 's/[(),]//g')
            centerdec=$(cat $resultfile |grep "Field center: (RA,Dec)" |awk '{print $6}' |sed -e 's/[(),]//g')

            if [ -z "$centerdec" ] ; then
                cp $resultfile ${DASCH_ASTROMETRY}/fail_${plate}_s${solutionNumber}.txt
            else
                solvedFlag=1
                break
            fi
        done

        if [ $solvedFlag -eq 0 ] ; then
            setFitWCS -e $prevSolutionNumber -b Astrometry2Succeeded -p $plate -c
            echo "No result for $plate"
            rm -f $xylsfile1
            exit 30  # signal a "managed failure" to the job management system
        fi

        rm -f ${DASCH_ASTROMETRY}/fail_${plate}_s${solutionNumber}.txt

        setFitWCS -e $prevSolutionNumber -b Astrometry2Succeeded -p $plate -s
        solutiontext=$(grep "Solution found"  $resultfile)
        echo "$solutiontext for case $casenumber"

        listhead ${DASCH_ASTROMETRY}/${rootname}_s${solutionNumber}.wcs >${rootname}_s${solutionNumber}.gethead
        naxis1=$(cat ${rootname}_s${solutionNumber}.gethead |grep 'IMAGEW' |awk '{print $3}')
        naxis2=$(cat ${rootname}_s${solutionNumber}.gethead |grep 'IMAGEH' |awk '{print $3}')
        ctype1=$(cat ${rootname}_s${solutionNumber}.gethead |grep 'CTYPE1' |awk '{print $3}' |sed 's/\-SIP//g' |sed "s/'//g")
        ctype2=$(cat ${rootname}_s${solutionNumber}.gethead |grep 'CTYPE2' |awk '{print $3}' |sed 's/\-SIP//g' |sed "s/'//g")
        crval1=$(cat ${rootname}_s${solutionNumber}.gethead |grep 'CRVAL1' |awk '{print $3}')
        crval2=$(cat ${rootname}_s${solutionNumber}.gethead |grep 'CRVAL2' |awk '{print $3}')
        crpix1=$(cat ${rootname}_s${solutionNumber}.gethead |grep 'CRPIX1' |awk '{print $3}')
        crpix2=$(cat ${rootname}_s${solutionNumber}.gethead |grep 'CRPIX2' |awk '{print $3}')
        cd1_1=$(cat ${rootname}_s${solutionNumber}.gethead |grep 'CD1_1'  |awk '{print $3}')
        cd1_2=$(cat ${rootname}_s${solutionNumber}.gethead |grep 'CD1_2'  |awk '{print $3}')
        cd2_1=$(cat ${rootname}_s${solutionNumber}.gethead |grep 'CD2_1'  |awk '{print $3}')
        cd2_2=$(cat ${rootname}_s${solutionNumber}.gethead |grep 'CD2_2'  |awk '{print $3}')
        scalex=0 # newer Astrometry.Net does not print this info
        scaley=0

        echo "plate ra dec scalex scaley binning naxis1 naxis2 ctype1 ctype2 crval1 crval2 crpix1 crpix2 cd1_1 cd1_2 cd2_1 cd2_2" >${DASCH_ASTROMETRY}/${plate}_s${solutionNumber}.tmp
        echo "----- -- --- ------ ------ ------- ------ ------ ------ ------ ------ ------ ------ ------ ----- ----- ----- -----" >>${DASCH_ASTROMETRY}/${plate}_s${solutionNumber}.tmp
        echo "$plate $centerra $centerdec $scalex $scaley $binning $naxis1 $naxis2 $ctype1 $ctype2 $crval1 $crval2 $crpix1 $crpix2 $cd1_1 $cd1_2 $cd2_1 $cd2_2" >>${DASCH_ASTROMETRY}/${plate}_s${solutionNumber}.tmp

        cat ${DASCH_ASTROMETRY}/${plate}_s${solutionNumber}.tmp |awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' >${DASCH_ASTROMETRY}/${plate}_s${solutionNumber}.db
        rm ${DASCH_ASTROMETRY}/${plate}_s${solutionNumber}.tmp
        rm ${rootname}_s${solutionNumber}.gethead
        echo "Solved: ra $centerra dec $centerdec for $plate"
        echo $casenumber >${DASCH_ASTROMETRY}/${rootname}_s${solutionNumber}.casenumber.txt

        rm -f ${rootname}_s${solutionNumber}.axy
        rm -f ${rootname}_s${solutionNumber}.match
        rm -f ${rootname}_s${solutionNumber}.rdls
        rm -f ${rootname}_s${solutionNumber}.corr
        rm -f ${rootname}_s${solutionNumber}-indx.xyls
        rm -f $resultfile
        rm -f $xylsfile
        rm -f $xylsfile2
    else
        echo "WARNING $xylsfile1 not found"
    fi
done

#  Feb 10, 2009 Edward J. Los - adapted from find_astrometry.csh
#  Apr 20, 2009 Edward J. Los - revise to work as a subroutine for filter_multiple.csh
#  May  1, 2009 Edward J. Los - do not delete the input "NONE" listing if there was a success
#  Aug 28, 2009 Edward J. Los - change _e to _s to reflect difference between solutionNumber and exposureNumber
#  Nov  3, 2009 Edward J. Los - Make this script independent of filter_multiple.csh
#  Mar 29, 2010 Edward J. Los - Do not over-write an existing solution
#  May 17, 2011 Edward J. Los - Modify for astrometry.net-0.38   Try using the internal timer
#  Jul  5, 2011 Edward J. Los - Extend cpulimit to 42 minutes except for cases 6,7, and 8
#  Mar 19, 2013 Edward J. Los - Do not run this script if we are using find_astrometry01.csh
#  Dec  1, 2014 Edward J. Los - Support the Astrometry2Succeeded bit in FitWCS of the mosaics table
