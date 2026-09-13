#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Astrometry correction of astrometry.net results for DASCH
#
# This script takes the results of find_astrometry01.sh, namely the star catalog,
# the source list, and the configuration file and provides them as inputs to scamp for
# a first pass solution based on UCAC4.
#
# TODO: annotate inputs/outputs; maybe also merge with `run_scamp.sh`.
#
# Invokes:
#
# - create_mosaic_scamp.py
# - getdirectory
# - getlocation
# - setScampRMS
# - sexcat_ndet.py
#
# External dependencies:
#
# - apps (FitsHeader, TranslateWCS)
# - SCAMP
#
# Variables used:
#
# - $DASCH_ASTROMETRY
# - $DASCH_MATCH
# - $DASCH_SCAMPIMAGES
# - $DASCH_SCRIPTS
# - $DASCH_SCRATCH

if [ $# -ne 2 ] ; then
    echo "usage: run_scamp_01.sh <solutionNumber> <list>"
    exit 1
fi

solutionNumber=$1
imagelist=$2

redo_filt=1
useTranslateWCS=1

if [ $solutionNumber -eq 0 ] ; then
    solutionString=""
else
    solutionString="_s$solutionNumber"
fi

echo "Processing $imagelist"

for mosaic in $(cat $imagelist) ; do
    #exit quietly if not a patrol plate
    getlocation -n $mosaic
    if [ $? -ne 0 ] ; then
        patrolPlate=0
    else
        patrolPlate=1
    fi

    if [ $patrolPlate -eq 0 ] ; then
        #echo "run_scamp_01.csh not patrol: $mosaic"
        continue;
    fi

    if [ $solutionNumber -eq 0 ] ; then
        mosaic_dir=$(getdirectory $mosaic)
        mosaic_fit=${mosaic_dir}/${mosaic}.fit
        mosaic_fit_tnx=${mosaic_dir}/${mosaic}_tnx.fit
        headeroutput=$DASCH_SCRATCH/${mosaic}_tnx.fit

        if [ ! -e $mosaic_fit -a ! -e $mosaic_fit_tnx ] ; then
            echo "ERROR run_scamp_01 can not find mosaic $mosaic_fit or $mosaic_fit_tnx"
            continue
        fi

        if [ -e $headeroutput ] ; then
            echo "ERROR output file $headeroutput already exists"
            continue
        fi

        if [ -e $mosaic_fit_tnx ] ; then
            imwcsinput=$mosaic_fit_tnx
        elif [ -e $mosaic_fit ] ; then
            imwcsinput=$mosaic_fit
        fi
    else
        echo "run_scamp_01.csh is not needed for solution $solutionNumber for $mosaic"
        continue
    fi

    echo "Run SCAMP 01 to correct distortions. Begin Process $(date) mosaic = ${mosaic}${solutionString}"

    #----------
    # Set parameters

    script_output=$(getlocation $mosaic -e $solutionNumber)
    scale=$(echo "$script_output" |awk '{print $4}')
    if [ $scale == "" ] ; then
        #echo "ERROR getlocation failed with args $mosaic -e $solutionNumber"
        continue
    fi

    aper_arcsec=10
    aper_pix=$(echo $aper_arcsec $scale |awk '{print $1/$2}')

    # Global parameters
    detthresh=5        #v1=5 v2=2
    maglim=8.0,16.0    #v1=v2=6,16
    snthresh=5.0,100.0 #v1=10,100 v2=5,100
    astrclip=5.0
    patternmatch=N     # done in the previous steps
    distortdegree=6   #v1=6 v2=5 v3=5
    #echo "ERROR: distortdegree is $distortdegree"
    ndet_max=10000     #v1=v2=10000
    #nref_max=10000    # will be determined after by:
    nref_max_factor=1  # nref_max = ndet_max / nref_max_factor
    ucac4_mag=MagM #v1=MagM v2=MagB
    ucac4_mag_imcat=m7
    match_radius_pix=10

    series=$(getlocation -s $mosaic |awk '{print $1}')
    if [ $series == "ac" -o $series == "am" -o $series == "ca" ] ; then
        match_radius_pix=40
    fi

    if [ $series == "ax" -o $series == "ay" ] ; then
        match_radius_pix=30
    fi

    match_radius=$(echo $match_radius_pix $scale |awk '{print $1*$2}')
    echo "Pixel Scale:  $scale Max distortion expected:  $match_radius N reference max:  $ndet_max for ${mosaic}${solutionString}"

    # take list from astrometry.net
    sex_cat=${DASCH_ASTROMETRY}/${mosaic}_full_scamp
    scampobjectfile=${DASCH_ASTROMETRY}/${mosaic}_full_scamp.object
    scampreffile=${DASCH_ASTROMETRY}/${mosaic}_full_scamp.ref
    scampconffile=${DASCH_ASTROMETRY}/${mosaic}_full_scamp.conf

    if [ ! -e ${DASCH_ASTROMETRY}/${mosaic}_full.db ] ; then
        echo "ERROR: run_scamp_01 ${DASCH_ASTROMETRY}/${mosaic}_full.db does not exist"
        continue
    fi

    if [ ! -e ${DASCH_ASTROMETRY}/${mosaic}_full_scamp.object ] ; then
        echo "ERROR: run_scamp_01 ${DASCH_ASTROMETRY}/${mosaic}_full_scamp.object does not exist"
        continue
    fi

    if [ ! -e ${DASCH_ASTROMETRY}/${mosaic}_full_scamp.ref ] ; then
        echo "ERROR: run_scamp_01 ${DASCH_ASTROMETRY}/${mosaic}_full_scamp.ref does not exist"
        continue
    fi

    #----------
    # Run Scamp

    checkplot_name=${DASCH_SCAMPIMAGES}/${mosaic}${solutionString}_astrfull_fgroups,${DASCH_SCAMPIMAGES}/${mosaic}${solutionString}_astrfull_distort,${DASCH_SCAMPIMAGES}/${mosaic}${solutionString}_astrfull_referr2d,${DASCH_SCAMPIMAGES}/${mosaic}${solutionString}_astrfull_referr1d

    if [ $solutionNumber -eq 0 -o $useTranslateWCS -eq 0 ] ; then
        while [ $distortdegree -gt 1 ] ; do
            set -x

            scamp \
                -c $scampconffile \
                -ASTREF_CATALOG FILE \
                -ASTREFCAT_NAME $scampreffile \
                -CROSSID_RADIUS $match_radius \
                -MATCH $patternmatch \
                -MATCH_RESOL $match_radius \
                -DISTORT_DEGREES $distortdegree \
                -ASTREFMAG_LIMITS $maglim \
                -SN_THRESHOLDS $snthresh \
                -ASTRCLIP_NSIGMA $astrclip \
                -CHECKPLOT_TYPE FGROUPS,DISTORTION,ASTR_REFERROR2D,ASTR_REFERROR1D \
                -CHECKPLOT_NAME $checkplot_name \
                -XML_NAME $DASCH_MATCH/${mosaic}${solutionString}_scamp_results.xml \
                -MERGEDOUTCAT_TYPE FITS_LDAC \
                -MERGEDOUTCAT_NAME $DASCH_MATCH/${mosaic}${solutionString}_scamp.refcat_used.fits \
                -HEADER_SUFFIX .head \
                -HEADER_TYPE NORMAL \
                -VERBOSE_TYPE QUIET \
                $scampobjectfile

            set +x

            if [ $? -eq 139 ] ; then
                echo "ERROR: run_scamp_01 segmentation fault for ${mosaic}${solutionString}"
            else
                if [ ! -e ${sex_cat}.head ] ; then
                    break
                fi

                FitsHeader -e 0 -i ${imwcsinput} -t ${sex_cat}.head -o ${sex_cat}.head2

                # -m mode: "check header output for scamp for inaccuracy"
                setScampRMS -p $mosaic -e $solutionNumber -m ${sex_cat}.head2

                if [ $? -eq 0 ] ; then
                    break
                fi
            fi

            ((distortdegree=$distortdegree - 1))
        done
    else
        set -x
        TranslateWCS -i $imwcsinput -o $headeroutputtmp -b $prevTnxFile
        FitsHeader -e $solutionNumber -i  $headeroutputtmp -o $headeroutput
        rm -f $headeroutputtmp
        set +x
    fi

    if [ -e $DASCH_MATCH/${mosaic}${solutionString}_scamp.refcat_used_1.fits ] ; then
        set -x
        nref_used=$($DASCH_SCRIPTS/sexcat_ndet.py $DASCH_MATCH/${mosaic}${solutionString}_scamp.refcat_used_1.fits)
        set +x

        if [ -e ${sex_cat}.head ] ; then
            rms_ra_deg=$(cat ${sex_cat}.head |grep ASTRRMS1= |awk '{print $2}')
            rms_ra=$(echo $rms_ra_deg 3600 $scale |awk '{print $1*$2/$3}')
            rms_dec_deg=$(cat ${sex_cat}.head |grep ASTRRMS2= |awk '{print $2}')
            rms_dec=$(echo $rms_dec_deg 3600 $scale |awk '{print $1*$2/$3}')

            echo "run_scamp_01 used $nref_used reference stars RMS RA = $rms_ra pixels RMS Dec = $rms_dec pixels distortdegree $distortdegree for ${mosaic}${solutionString}"

            # -b mode: "declare new algorithm SCAMP multiple exposure astrometry.wcs success"
            setScampRMS -b -e $solutionNumber -p $mosaic
            if [ $? -ne 0 ] ; then
                continue
            fi

            if [ $solutionNumber -eq 0 ] ; then
                # we are done.  The  save_header.csh file will write the "A" keywords
                set -x
                $DASCH_SCRIPTS/create_mosaic_scamp.py ${imwcsinput} ${sex_cat}.head ${headeroutput} $solutionNumber
                set +x
            else
                # pyfits can not handle headers without data properly
                # Use FitsHeader to write keywords "B" and up.
                rm -f ${headeroutput}
                set -x
                FitsHeader -p 16 -i $imwcsinput -t ${sex_cat}.head -o $headeroutput -e $solutionNumber
                set +x
            fi

            chmod g+w ${headeroutput}
        else
            echo "ERROR: File does not exist: ${sex_cat}.head"
            echo "ERROR EXIT"
            exit
        fi
    else
        if [ $solutionNumber -eq 0 -o $useTranslateWCS -eq 0 ] ; then
            echo "ERROR: File does not exist: $DASCH_MATCH/${mosaic}${solutionString}_scamp.refcat_used_1.fits"
        fi
    fi

    rm -f $DASCH_MATCH/${mosaic}${solutionString}_scamp_results.xml
    rm -f $DASCH_MATCH/${mosaic}${solutionString}_scamp.refcat_used_1.fits
    rm -f $DASCH_ASTROMETRY/${mosaic}${solutionString}_full_scamp.head
    rm -f $DASCH_ASTROMETRY/${mosaic}${solutionString}_full_scamp.head2
    rm -f $DASCH_ASTROMETRY/${mosaic}${solutionString}_full_scamp.head_cor
done

#  Feb 20, 2013 Edward J. Los     Adapted from Mathieu Servillat's run_scamp.csh
#  Mar  4, 2013 Edward J. Los     If an InaccuratePV error occurs within scamp, then reduce the polynomial order
#  Jul 23, 2013 Edward J. Los - Exit quietly if not a patrol plate
#  Apr 15, 2014 Edward J. Los - Correct module name in the error messages
