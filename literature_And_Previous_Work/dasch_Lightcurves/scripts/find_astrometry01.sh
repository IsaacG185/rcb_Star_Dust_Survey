#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Register images automatically using Astrometry.Net. This script is part of the
# "prep-WW" pipeline stage that derives the "TNX" distortion calibrations.
#
# TODO: annotate inputs/outputs; maybe also merge with other
# `find_astrometry*.sh` scripts.
#
# Invokes:
#
# - db2xyls
# - getdirectory
# - getlocation
# - getseries
# - setScampRMS
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
# - $DASCH_MATCH
# - $DASCH_SCRIPTS

solutionNumber=$1
imagelist=$2

binning=1
debugMode=0
fractionFLAG="-f 0.01"
kronFLAG=""
configSuffix=1
((prevSolutionNumber=$solutionNumber - 1))

echo "Processing ${imagelist}"

if [ $solutionNumber -eq 0 ] ; then
    solutionString=""
else
    solutionString="_s$solutionNumber"
fi
#echo "solutionNumber $solutionNumber solutionString $solutionString"

for plate in $(cat $imagelist) ; do
    #exit quietly if not a patrol plate
    getlocation -n $plate

    if [ $? -ne 0 ] ; then
        #echo "find_astrometry01.csh not patrol: $plate"
        continue
    fi

    solvedFlag=0
    mosaicname=$(echo "${plate}${solutionString}.fit")
    directory=$(getdirectory $plate)
    filename=${directory}/${mosaicname}
    rootname=$(echo "$mosaicname" |sed -e 's/.fit//g')
    resultroot=$(echo $rootname |sed -e 's/r90ww//g' -e 's/r180ww//g' -e 's/r270ww//g' -e 's/ww//g' |awk '{print $1}')

    scaletext=$(getseries $plate -f 1.1 -b $binning)
    minscale=$(echo $scaletext |awk '{print $1}')
    maxscale=$(echo $scaletext |awk '{print $3}')

    if [ ! -e $filename ] ; then
        filename=${directory}/${rootname}_tnx.fit
    fi

    echo "mosaicname $mosaicname rootname $rootname resultroot $resultroot"

    if [ -e ${DASCH_ASTROMETRY}/${plate}${solutionString}_full.db ] ; then
        echo "ERROR: find_astrometry01 ${DASCH_ASTROMETRY}/${plate}${solutionString}_full.db already exists"
        continue
    fi

    if [ -e ${DASCH_ASTROMETRY}/${plate}${solutionString}_full_scamp.object ] ; then
        echo "ERROR: find_astrometry01 ${DASCH_ASTROMETRY}/${plate}${solutionString}_full_scamp.object already exists"
        continue
    fi

    if [ -e ${DASCH_ASTROMETRY}/${plate}${solutionString}_full_scamp.ref ] ; then
        echo "ERROR: find_astrometry01 ${DASCH_ASTROMETRY}/${plate}${solutionString}_full_scamp.ref already exists"
        continue
    fi

    if [ -e ${DASCH_ASTROMETRY}/${plate}${solutionString}_full_scamp.conf ] ; then
        echo "ERROR: find_astrometry01 ${DASCH_ASTROMETRY}/${plate}${solutionString}_full_scamp.conf already exists"
        continue
    fi

    # -f mode: "check and set new algorithm failure"
    setScampRMS -f -e $solutionNumber -p $plate
    if [ $? -ne 0 ] ; then
        continue
    fi

    if [ $solutionNumber -gt 0 ] ; then
        filename=${DASCH_ASTROMETRY}/${rootname}_none.db
        directory="$DASCH_ASTROMETRY"
    fi

    rm -f ${DASCH_ASTROMETRY}/fail_${plate}${solutionString}_full.txt

    if [ -e $filename ] ; then
        sextractorfile=${DASCH_MATCH}/${plate}${solutionString}_tnx.db
        if [ ! -e $sextractorfile ] ; then
            sextractorfile=${DASCH_MATCH}/${plate}${solutionString}.db
            if [ ! -e $sextractorfile ] ; then
                sextractorfile=${DASCH_MATCH}/${resultroot}.db
            fi
        fi
        echo "sextractorfile $sextractorfile"

        maskfile=${DASCH_MATCH}/${plate}${solutionString}_close.db
        if [ ! -e $maskfile ] ; then
            maskfile=${DASCH_MATCH}/${resultroot}_close.db
        fi
        echo "maskfile $maskfile"

        if [ -e $maskfile ] ; then
            maskflag="-m 0 $maskfile"
        else
            maskflag=""
        fi

        if [ $solutionNumber -eq 0 ] ; then
            xylsfile=$(echo "$filename" |sed -e 's/\.fit$/\_full.xyls/' |sed 's/_tnx//g')
            xylsfile1=${xylsfile}first
            xylsfile2=${xylsfile}second
            resultfile=$(echo "$filename" |sed -e 's/\.fit$/\_full.solve/' |sed 's/_tnx//g')

            scampobjectfile=$(echo "$filename" |sed -e 's/\.fit$/\_full_scamp.object/' |sed 's/_tnx//g')
            scampreffile=$(echo "$filename" |sed -e 's/\.fit$/\_full_scamp.ref/' |sed 's/_tnx//g')
            scampconffile=$(echo "$filename" |sed -e 's/\.fit$/\_full_scamp.conf/' |sed 's/_tnx//g')
        else
            xylsfile=${DASCH_ASTROMETRY}/${rootname}_full.xyls
            xylsfile1=${DASCH_ASTROMETRY}/${rootname}_none.db
            xylsfile2=${xylsfile}second

            scampobjectfile=${DASCH_ASTROMETRY}/${rootname}_full_scamp.object
            scampreffile=${DASCH_ASTROMETRY}/${rootname}_full_scamp.ref
            scampconffile=${DASCH_ASTROMETRY}/${rootname}_full_scamp.conf
            resultfile=${DASCH_ASTROMETRY}/${rootname}_full.solve
        fi

        if [ $debugMode -eq 1 ] ; then
            debugfile=$(echo "$filename" |sed -e 's/\.fit$/\_full.debug/' |sed 's/_tnx//g')
            debugflag="-d $debugfile"
        else
            debugflag=""
        fi

        mosaicsize=$(getlocation -a $plate)
        width=$(echo "$mosaicsize" |awk '{print $1}')
        height=$(echo "$mosaicsize" |awk '{print $2}')

        if [ $solutionNumber -eq 0 ] ; then
            if [ ! -e $sextractorfile ] ; then
                echo "ERROR: sextractor execution script is out-of-date for ${plate}${solutionString}  rerun run_sextractor.csh"
                continue
            fi
        fi

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
                maxdepth="--depth 1-120"
                maxstars=480
                cpulimit="--cpulimit 300"
            fi
            if [ $casenumber -eq 7 ] ; then
                filterFLAG="-e"
                ellipticityFLAG="-l 0.3"
                maxdepth="--depth 1-120"
                maxstars=480
                cpulimit="--cpulimit 300"
            fi
            if [ $casenumber -eq 8 ] ; then
                filterFLAG=""
                ellipticityFLAG=""
                maxdepth="--depth 1-120"
                maxstars=480
                cpulimit="--cpulimit 300"
            fi

            if [ $solutionNumber -eq 0 ] ; then
                echo "db2xyls -v $maskflag -i $sextractorfile -o $xylsfile2 -w $width -h $height -n $maxstars $fractionFLAG  $ellipticityFLAG $kronFLAG $filterFLAG $debugflag"
                db2xyls \
                    $maskflag \
                    -i $sextractorfile \
                    -o $xylsfile2 \
                    -w $width \
                    -h $height \
                    -n $maxstars \
                    $fractionFLAG $ellipticityFLAG $kronFLAG $filterFLAG $debugflag
            else
                # use reverse binning of 16 to undo binning by filter_multiple
                echo "db2xyls $maskflag -r 16 -i $xylsfile1 -o $xylsfile2 -w $width -h $height -n $maxstars $fractionFLAG  $ellipticityFLAG $kronFLAG $filterFLAG $debugflag"
                db2xyls \
                    $maskflag \
                    -r 16 \
                    -i $xylsfile1 \
                    -o $xylsfile2 \
                    -w $width \
                    -h $height \
                    -n $maxstars \
                    $fractionFLAG $ellipticityFLAG  $kronFLAG $filterFLAG $debugflag
            fi

            cp $xylsfile2 $xylsfile

            if [ $debugMode -eq 1 ] ; then
                column -a -b <$debugfile >${DASCH_ASTROMETRY}/${rootname}_full.debug
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
                --width $width \
                --height $height \
                $cpulimit \
                $maxdepth \
                --no-plots \
                --scamp $scampobjectfile \
                --scamp-ref $scampreffile \
                --scamp-conf $scampconffile \
                >$resultfile 2>&1
            set +x
            date

            centerra=$(cat $resultfile |grep "Field center: (RA,Dec)" |awk '{print $5}' |sed -e 's/[(),]//g')
            centerdec=$(cat $resultfile |grep "Field center: (RA,Dec)" |awk '{print $6}' |sed -e 's/[(),]//g')

            if [ -z "$centerdec" ] ; then
                cp $resultfile ${DASCH_ASTROMETRY}/fail_${plate}${solutionString}_full.txt
            else
                solvedFlag=1
                break
            fi
        done

        if [ $solvedFlag -eq 0 ] ; then
            echo "No result for ${plate}${solutionString}"
        else
            rm -f ${DASCH_ASTROMETRY}/fail_${plate}${solutionString}_full.txt
            solutiontext=$(grep "Solution found"  $resultfile)
            echo "$solutiontext for case $casenumber"

            listhead ${directory}/${rootname}_full.wcs >${rootname}_full.gethead
            naxis1=$(cat ${rootname}_full.gethead |grep 'IMAGEW' |awk '{print $3}')
            naxis2=$(cat ${rootname}_full.gethead |grep 'IMAGEH' |awk '{print $3}')
            ctype1=$(cat ${rootname}_full.gethead |grep 'CTYPE1' |awk '{print $3}' |sed 's/\-SIP//g' |sed "s/'//g")
            ctype2=$(cat ${rootname}_full.gethead |grep 'CTYPE2' |awk '{print $3}' |sed 's/\-SIP//g' |sed "s/'//g")
            crval1=$(cat ${rootname}_full.gethead |grep 'CRVAL1' |awk '{print $3}')
            crval2=$(cat ${rootname}_full.gethead |grep 'CRVAL2' |awk '{print $3}')
            crpix1=$(cat ${rootname}_full.gethead |grep 'CRPIX1' |awk '{print $3}')
            crpix2=$(cat ${rootname}_full.gethead |grep 'CRPIX2' |awk '{print $3}')
            cd1_1=$(cat ${rootname}_full.gethead |grep 'CD1_1' |awk '{print $3}')
            cd1_2=$(cat ${rootname}_full.gethead |grep 'CD1_2' |awk '{print $3}')
            cd2_1=$(cat ${rootname}_full.gethead |grep 'CD2_1' |awk '{print $3}')
            cd2_2=$(cat ${rootname}_full.gethead |grep 'CD2_2' |awk '{print $3}')
            scalex=0 # newer Astrometry.Net does not print this info
            scaley=0

            echo "plate ra dec scalex scaley binning naxis1 naxis2 ctype1 ctype2 crval1 crval2 crpix1 crpix2 cd1_1 cd1_2 cd2_1 cd2_2" >${DASCH_ASTROMETRY}/${plate}${solutionString}_full.tmp
            echo "----- -- --- ------ ------ ------- ------ ------ ------ ------ ------ ------ ------ ------ ----- ----- ----- -----" >>${DASCH_ASTROMETRY}/${plate}${solutionString}_full.tmp
            echo "$plate $centerra $centerdec $scalex $scaley $binning $naxis1 $naxis2 $ctype1 $ctype2 $crval1 $crval2 $crpix1 $crpix2 $cd1_1 $cd1_2 $cd2_1 $cd2_2" >>${DASCH_ASTROMETRY}/${plate}${solutionString}_full.tmp
            cat ${DASCH_ASTROMETRY}/${plate}${solutionString}_full.tmp |awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' >${DASCH_ASTROMETRY}/${plate}${solutionString}_full.db
            rm ${DASCH_ASTROMETRY}/${plate}${solutionString}_full.tmp
            rm ${rootname}_full.gethead
            echo "Solved: ra $centerra dec $centerdec for ${plate}${solutionString}"

            if [ $solutionNumber -eq 0 ] ; then
                mv $resultfile ${DASCH_ASTROMETRY}/${rootname}_full.solve
                echo $casenumber >${DASCH_ASTROMETRY}/${rootname}_full.casenumber.txt

                if [ -e $scampobjectfile ] ; then
                    mv $scampobjectfile ${DASCH_ASTROMETRY}/${rootname}_full_scamp.object
                else
                    echo "ERROR: scamp object file $scampobjectfile does not exist"
                fi
                if [ -e $scampreffile ] ; then
                    mv $scampreffile ${DASCH_ASTROMETRY}/${rootname}_full_scamp.ref
                else
                    echo "ERROR: scamp ref file $scampreffile does not exist"
                fi
                if [ -e $scampconffile ] ; then
                    mv $scampconffile ${DASCH_ASTROMETRY}/${rootname}_full_scamp.conf
                else
                    echo "ERROR: scamp conf file $scampconffile does not exist"
                fi
            fi

            # -a mode: "declare new algorithm astrometry.wcs success"
            setScampRMS -a -e $solutionNumber -p $plate
        fi

        if [ -e ${rootname}_full.wcs ] ; then
            mv ${rootname}_full.wcs ${DASCH_ASTROMETRY}/${rootname}_full.wcs
        fi

        if [ $solutionNumber -eq 0 ] ; then
            if [ -e ${directory}/${rootname}_full.wcs ] ; then
                mv ${directory}/${rootname}_full.wcs ${DASCH_ASTROMETRY}/${rootname}_full.wcs
            fi
        fi

        #echo "ERROR: early exit: directory $directory rootname is $rootname filename is $filename xylsfile $xylsfile xylsfile1 $xylsfile1 xylsfile2 $xylsfile2 width $width height $height"
        #echo "full.axy ${rootname}_full.axy and ${directory}/${rootname}_full.axy"
        #echo "scampobjectfile $scampobjectfile"
        #exit

        if [ $solutionNumber -gt 0 ] ; then
            rm -f $scampobjectfile
            rm -f $scampreffile
            rm -f $scampconffile
        fi

        rm -f ${rootname}_full.axy
        rm -f ${rootname}_full.match
        rm -f ${rootname}_full.rdls
        rm -f ${rootname}_full.solved
        rm -f ${directory}/${rootname}_full.axy
        rm -f ${directory}/${rootname}_full.match
        rm -f ${directory}/${rootname}_full.rdls
        rm -f ${directory}/${rootname}_full.solved
        rm -f ${directory}/${rootname}_full-indx.xyls
        rm -f ${directory}/${rootname}_full.corr
        rm -f $scampobjectfile
        rm -f $scampreffile
        rm -f $scampconffile
        rm -f ${rootname}_full-indx.xyls
        rm -f $resultfile
        rm -f $xylsfile
        rm -f $xylsfile2
    else
        echo "ERROR: $filename not found"
    fi
done

#  Nov 30, 2010 Edward J. Los - Adapt for bin 1 mosaics
#  Aug 29, 2011 Edward J. Los - Modify for astrometry.net-0.38 and extend cpulimit
#  Mar  2, 2012 Edward J. Los - Adopt for use with multiple exposure masks
#  Feb 20, 2013 Edward J. Los - Enable SCAMP output files
#  Mar  4, 2013 Edward J. Los - Set maxdepth to 120 for cases 6, 7, and 8
#  Mar 19, 2013 Edward J. Los - Support multipe exposures
#  Jul 23, 2013 Edward J. Los - Exit quietly if not a patrol plate
#  Apr 15, 2014 Edward J. Los - Correct the error message module name
