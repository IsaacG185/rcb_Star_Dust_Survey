#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Astrometry correction of DASCH mosaics with SCAMP.
#
# Inputs:
#
# - $RAID/ExposureData/Mosaics/${series}/${id}_${mosnum}/${series}${id}_${mosnum}_01${rot}ww.fit
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${sol_tag}.db
# - if solnum > 0:
#   - $DASCH_HEADERS/${series}${id}_${mosnum}_01${rot}ww_tnx.hdr
#   - $DASCH_HEADERS/${series}${id}_${mosnum}_01${rot}ww_s${solnum}.hdr
#
# Outputs:
#
# - if solnum = 0:
#   - $RAID/ExposureData/Mosaics/${series}/${id}_${mosnum}/${series}${id}_${mosnum}_01${rot}ww_tnx.fit
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_astr_fgroups_1.png
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_astr_referr1d_1.png
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_astr_referr2d_1.png
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_astr_distort_1.png
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_drad_map.png
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_drad_offset.png
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww.refcat.ucac5.fits.reg
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_scamp.sexcat_filtered.fits.reg
#   - $DASCH_SCAMPIMAGES/${series}${id}_${mosnum}_01${rot}ww_scamp.sexcat_filtered.fits.xy.reg
# - else:
#   - $DASCH_HEADERS/${series}${id}_${mosnum}_01${rot}ww_s${solnum}_tnx.hdr
#
# Database updates:
#
# - Set `PolyRefStars`, `PolyRARMS`, `PolyDecRMS` in `scanner.mosaics` row.
#
# Invokes:
#
# - comparefits
# - create_mosaic_scamp.py
# - FitsHeader
# - getdirectory
# - getlocation
# - match_cat_ref.py
# - plot_drad.py
# - refcat2fit.py
# - setScampRMS
# - sexcat_filter.py
# - sexcat_ndet.py
# - tnxdrad_addcols.py
# - update_sextractor
# - votable
#
# External dependencies:
#
# - SCAMP
#
# Variables used:
#
# - $DASCH_HEADERS
# - $DASCH_MATCH
# - $DASCH_SCAMPIMAGES
# - $DASCH_SCRIPTS
# - $DASCH_SCRATCH

secondpass=0

if [ $# -eq 3 ] ; then
    if [ "$3" = 2 ] ; then
        #echo "second pass active"
        secondpass=1
    else
        echo "usage: $0 <solutionNumber> <list> [pass]"
        exit
    fi
elif [ $# -ne 2 ] ; then
    echo "usage: $0 <solutionNumber> <list> [pass]"
    exit
fi

solutionNumber=$1
imagelist=$2

# if plot_drad is set to zero, then changes are necessary in the setSetScampRMS functionality
plot_drad=1
redo_filt=1
useTranslateWCS=1
# The "pass two" flag is set to one if find_astrometry_01.csh has produced an output
pass_two_flag=0

if [ $solutionNumber -eq 0 ] ; then
    solutionString=""
else
    solutionString="_s$solutionNumber"
fi

echo "Processing ${imagelist}"

for mosaic in $(cat $imagelist) ; do
    getlocation -n $mosaic
    if [ $? -ne 0 ] ; then
        patrolPlate=0
    else
        patrolPlate=1
    fi

    if [ $secondpass -eq 1 ] ; then
        if [ $patrolPlate -eq 0 ] ; then
            continue
        fi
    fi

    pass_two_error=0

    if [ $solutionNumber -eq 0 ] ; then
        mosaic_dir=$(getdirectory $mosaic)
        mosaic_fit=${mosaic_dir}/${mosaic}.fit
        mosaic_fit_tnx=${mosaic_dir}/${mosaic}_tnx.fit
        mosaic_pass_one=$DASCH_SCRATCH/${mosaic}_tnx.fit

        if [ -e $mosaic_pass_one ] ; then
            imwcsinput=$mosaic_pass_one
            pass_two_flag=1
            headeroutput=${mosaic_fit_tnx}.provisional
        else
            imwcsinput=$mosaic_fit
            headeroutput=${mosaic_fit_tnx}
        fi
    else
        imwcsinput=${DASCH_HEADERS}/${mosaic}${solutionString}.hdr
        headeroutput=${DASCH_HEADERS}/${mosaic}${solutionString}_tnx.hdr
        headeroutputtmp=${DASCH_HEADERS}/${mosaic}${solutionString}_tnx.hdrtmp
        mosaic_fit_tnx=${DASCH_MATCH}/${mosaic}${solutionString}_dummy_file_to_avoid_undefined_variable

        ((prevSolutionNumber=$solutionNumber - 1))
        if [ $prevSolutionNumber -eq 0 ] ; then
            prevSolutionString=""
        else
            prevSolutionString="_s$prevSolutionNumber"
        fi

        prevTnxFile=${DASCH_HEADERS}/${mosaic}${prevSolutionString}_tnx.hdr
        pass_two_result=${DASCH_ASTROMETRY}/${mosaic}${solutionString}_full.db

        if [ -e $pass_two_result ] ; then
            # -q mode: "declare new algorithm SCAMP USNO success"
            setScampRMS -q -e $solutionNumber -p $mosaic
            if [ $? -ne 0 ] ; then
                continue
            fi

            pass_two_flag=1
            imwcsinput=$pass_two_result
        fi
    fi

    if [ ! -e $imwcsinput ] ; then
        if [ $secondpass -eq 1 ] ; then
            echo "ERRORD: run_scamp cannot find mosaic $imwcsinput"
        else
            echo "ERROR: run_scamp cannot find mosaic $imwcsinput"
        fi

        continue
    fi

    echo "Run SCAMP to correct distortions. Begin Process $(date) mosaic = ${mosaic}${solutionString}"

    #----------
    # Set parameters

    script_output=$(getlocation  $mosaic -e $solutionNumber)
    scale=$(echo "$script_output" |awk '{print $4}')
    if [ -z "$scale" ] ; then
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
    distortdegree=6    #v1=6 v2=5 v3=5
    ndet_max=10000     #v1=v2=10000
    #nref_max=10000    # will be determined after by:
    nref_max_factor=1  # nref_max=ndet_max / nref_max_factor
    ucac5_mag=MagM
    ucac5_mag_imcat=m7
    testflag=""
    match_radius_pix=10

    series=$(getlocation -s $mosaic |awk '{print $1}')
    if [ $series = ac -o $series = am -o $series == ca ] ; then
        match_radius_pix=40
    elif [ $series = ax -o $series == ay ] ; then
        match_radius_pix=30
    fi

    #The following test performed on July 31, 2018 did not work.
    #if [ $patrolPlate -eq 1 ] ; then
    #    match_radius_pix=100
    #fi

    match_radius=$(echo $match_radius_pix $scale |awk '{print $1*$2}')

    echo "Pixel Scale:  $scale Max distortion expected:  $match_radius N reference max:  $ndet_max for ${mosaic}${solutionString}"

    #----------
    # Filter sextractor results for scamp input
    #
    # take list from run_sextractor.csh:
    sex_cat_init=$DASCH_MATCH/${mosaic}${solutionString}.db

    # Reduce file size and sort with starbase, which is faster then python
    if [ ! -e ${sex_cat_init} ] ; then
        if [ $pass_two_flag -eq 1 ] ; then
            sex_cat_init=$DASCH_MATCH/${mosaic}${solutionString}_tnx.db
            if [ ! -e ${sex_cat_init} ] ; then
                echo "ERROR: run_scamp cannot find sextractor file ${sex_cat_init}"
                continue
            fi
        else
            echo "ERROR: run_scamp cannot find sextractor file ${sex_cat_init}"
            continue
        fi
    fi

    sex_cat_init_tmp=$DASCH_MATCH/${mosaic}${solutionString}.tmp
    sex_cat=$DASCH_MATCH/${mosaic}${solutionString}.sexcat_filtered
    cmd="$DASCH_SCRIPTS/sexcat_filter.py ${sex_cat_init_tmp} $ndet_max ${sex_cat}.fits $imwcsinput "

    if [ $solutionNumber -eq 0 -o $useTranslateWCS -eq 0 ] ; then
        if [ ! -e ${sex_cat}.fits ] ; then
            set -x
            column \
                -i ${sex_cat_init} \
                X_IMAGE Y_IMAGE ERRA_IMAGE ERRB_IMAGE FLUX_ISO FLUXERR_ISO BFLAGS MAG_ISO ra dec \
                |sorttable -n MAG_ISO \
                |head -n 11000 \
                >$sex_cat_init_tmp
            $cmd
            set +x
        elif [ $redo_filt -eq 1 ] ; then
            set -x
            column \
                -i ${sex_cat_init} \
                X_IMAGE Y_IMAGE ERRA_IMAGE ERRB_IMAGE FLUX_ISO FLUXERR_ISO BFLAGS MAG_ISO ra dec \
                |sorttable -n MAG_ISO \
                |head -n 11000 \
                >$sex_cat_init_tmp
            $cmd
            set +x
        else
            echo "WARNING: Keep current filtered file ${sex_cat}.fits"
        fi
    fi

    #----------
    # Reference catalog for the field
    # Determine nref from ndet
    ref_cat=$DASCH_MATCH/${mosaic}${solutionString}.refcat.ucac5

    if [ $solutionNumber -eq 0 -o $useTranslateWCS -eq 0 ] ; then
        echo "$DASCH_SCRIPTS/sexcat_ndet.py ${sex_cat}.fits"
        ndet=$($DASCH_SCRIPTS/sexcat_ndet.py ${sex_cat}.fits)
        nref_max=$ndet
        if [ $nref_max_factor -ne 1 ] ; then
            nref_max=$(echo $ndet $nref_max_factor |awk '{print $1/$2}')
        fi

        nref=$(imcat -f -c ucac5 -m1 $maglim $imwcsinput |awk '{print $2}')
        nref_imcat=$(echo $nref_max 3 |awk '{print $1*$2}')
        echo "$nref reference stars available, keep at max $nref_max for ${mosaic}${solutionString}"

        cmd="imcat -dh -c ucac5 -mx $maglim -s $ucac5_mag_imcat -n $nref_imcat $imwcsinput"
        echo $cmd '>' ${ref_cat}'.db'
        $cmd > ${ref_cat}.db
        cmd="$DASCH_SCRIPTS/refcat2fit.py ${ref_cat}.db ${ref_cat}.fits $nref_max $ucac5_mag"
        echo $cmd
        $cmd
    fi

    #----------
    # Run Scamp

    checkplot_name=${DASCH_SCAMPIMAGES}/${mosaic}${solutionString}${testflag}_astr_fgroups,${DASCH_SCAMPIMAGES}/${mosaic}${solutionString}${testflag}_astr_distort,${DASCH_SCAMPIMAGES}/${mosaic}${solutionString}${testflag}_astr_referr2d,${DASCH_SCAMPIMAGES}/${mosaic}${solutionString}${testflag}_astr_referr1d
    set -x

    if [ $solutionNumber -eq 0 -o $useTranslateWCS -eq 0 ] ; then
        while [ $distortdegree -gt 1 ] ; do
            scamp_options="-c $DASCH_SCRIPTS/Sextractor/scamp.conf -ASTREF_CATALOG FILE -ASTREFCAT_NAME ${ref_cat}.fits -CROSSID_RADIUS $match_radius -MATCH $patternmatch -MATCH_RESOL $match_radius -DISTORT_DEGREES $distortdegree -ASTREFMAG_LIMITS $maglim -SN_THRESHOLDS $snthresh -ASTRCLIP_NSIGMA $astrclip -CHECKPLOT_NAME $checkplot_name -XML_NAME $DASCH_MATCH/${mosaic}${solutionString}_scamp_results.xml -MERGEDOUTCAT_NAME $DASCH_MATCH/${mosaic}${solutionString}_scamp.refcat_used.fits -VERBOSE_TYPE QUIET"
            cmd="scamp $scamp_options ${sex_cat}.fits"
            $cmd

            if [ $? -eq 139 ] ; then
                echo "ERROR: run_scamp segmentation fault for ${mosaic}${solutionString}"
                pass_two_error=1
            else
                if [ ! -e ${sex_cat}.head ] ; then
                    break
                fi

                #echo "FitsHeader -i ${imwcsinput} -t ${sex_cat}.head -o ${sex_cat}.head2"
                FitsHeader -e 0 -i ${imwcsinput} -t ${sex_cat}.head -o ${sex_cat}.head2
                # -m mode: "check header output for scamp for inaccuracy"
                setScampRMS -p $mosaic -e $solutionNumber -m ${sex_cat}.head2

                if [ $? -eq 0 ] ; then
                    break
                fi
            fi

            ((distortdegree-=1))
        done

        #echo "ERROR: early exit"
        #exit
        #echo Astrometry keywords:
        #cat ${sex_cat}.head
    else
        rm -f $headeroutputtmp

        TranslateWCS -i $imwcsinput -o $headeroutputtmp -b $prevTnxFile

        if [ $? -eq 0 ] ; then
            FitsHeader -e $solutionNumber -i  $headeroutputtmp -o $headeroutput

            if [ $? -eq 0 -a $pass_two_flag -ne 0 ] ; then
                # -l mode: "declare new algorithm SCAMP UCAC success for multiple exposures"
                setScampRMS -l -p $mosaic -e $solutionNumber

                if [ $? -eq 0 ] ; then
                    UpdateMosaicTable -e $solutionNumber -m $headeroutput
                fi
            fi
        fi

        rm -f $headeroutputtmp
    fi

    set +x

    if [ -e $DASCH_MATCH/${mosaic}${solutionString}_scamp.refcat_used_1.fits ] ; then
        echo "$DASCH_SCRIPTS/sexcat_ndet.py $DASCH_MATCH/${mosaic}${solutionString}_scamp.refcat_used_1.fits where distortdegree is $distortdegree"
        nref_used=$($DASCH_SCRIPTS/sexcat_ndet.py $DASCH_MATCH/${mosaic}${solutionString}_scamp.refcat_used_1.fits)

        if [ -e ${sex_cat}.head ] ; then
            rms_ra_deg=$(cat ${sex_cat}.head |grep ASTRRMS1= |awk '{print $2}')
            rms_ra=$(echo $rms_ra_deg 3600 $scale |awk '{print $1*$2/$3}')
            rms_dec_deg=$(cat ${sex_cat}.head |grep ASTRRMS2= |awk '{print $2}')
            rms_dec=$(echo $rms_dec_deg 3600 $scale |awk '{print $1*$2/$3}')
            #echo "Used $nref_used reference stars RMS RA = $rms_ra pixels RMS Dec = $rms_dec pixels for ${mosaic}${solutionString}"
            #Use a negative count for the first insertion
            if [ $pass_two_flag -eq 0 ] ; then
                # "original" mode: save SCAMP info to database
                setScampRMS -p $mosaic -e $solutionNumber -n  -${nref_used} -r $rms_ra -d $rms_dec
            fi

            #----------
            # Create new header and new mosaic
            #echo Create new header and mosaic
            if [ $solutionNumber -eq 0 ] ; then
                # we are done.  The  save_header.csh file will write the "A" keywords
                cmd="$DASCH_SCRIPTS/create_mosaic_scamp.py ${imwcsinput} ${sex_cat}.head ${headeroutput} $solutionNumber"
                echo $cmd
                $cmd
            else
                # pyfits can not handle headers without data properly
                # Use FitsHeader to write keywords "B" and up.
                rm -f ${headeroutput}
                cmd="FitsHeader -p 16 -i $imwcsinput -t ${sex_cat}.head -o  $headeroutput -e $solutionNumber"
                echo $cmd
                $cmd
            fi

            chmod g+w ${headeroutput}

            #cp -f $DASCH_SCRIPTS/${mosaic}${solutionString}_tnx.hdr $DASCH_HEADERS
            #echo Header saved in $DASCH_HEADERS

            #----------
            if [ $plot_drad -eq 1 ] ; then
                # Test: rerun sextractor on mosaic and match with reference catalog
                # would be better to get new ra,dec form x,y of the first sextractor run
                #echo
                #echo ----------
                #echo Run update_sextractor to get new RA and Dec
                sex_cat_init_tmp2=$DASCH_MATCH/${mosaic}${solutionString}.tmp2
                sex_cat_init_tmp3=$DASCH_MATCH/${mosaic}${solutionString}_scamp.db

                if [ $pass_two_flag -eq 1 ] ; then
                    # -m mode: "check header output for scamp for inaccuracy"
                    setScampRMS -p $mosaic -e $solutionNumber -m $headeroutput
                    if [ $? -ne 0 ] ; then
                        pass_two_error=1
                    fi
                fi

                if [ $pass_two_error -eq 0 ] ; then
                    cmd="update_sextractor -v -e $solutionNumber -i ${sex_cat_init} -o $sex_cat_init_tmp2 -m  $headeroutput "

                    echo $cmd
                    $cmd
                    echo "column -i $sex_cat_init_tmp2 X_IMAGE Y_IMAGE ERRA_IMAGE ERRB_IMAGE FLUX_ISO FLUXERR_ISO BFLAGS MAG_ISO ra dec |sorttable -n MAG_ISO |head -n 11000 > $sex_cat_init_tmp3"
                    column \
                        -i $sex_cat_init_tmp2 \
                        X_IMAGE Y_IMAGE ERRA_IMAGE ERRB_IMAGE FLUX_ISO FLUXERR_ISO BFLAGS MAG_ISO ra dec \
                        |sorttable -n MAG_ISO \
                        |head -n 11000 \
                        >$sex_cat_init_tmp3

                    #----------
                    # Filter sextractor results
                    cmd="$DASCH_SCRIPTS/sexcat_filter.py $sex_cat_init_tmp3 $ndet_max $DASCH_MATCH/${mosaic}${solutionString}_scamp.sexcat_filtered.fits $headeroutput"
                    echo $cmd
                    $cmd

                    #----------
                    # Match with reference catalog
                    cmd="$DASCH_SCRIPTS/match_cat_ref.py $DASCH_MATCH $DASCH_SCAMPIMAGES ${mosaic}${solutionString}.sexcat_filtered.fits ${mosaic}${solutionString}.refcat.ucac5.fits ${mosaic}${solutionString}_match.fits $match_radius $scale"
                    echo $cmd
                    $cmd

                    cmd="$DASCH_SCRIPTS/match_cat_ref.py $DASCH_MATCH $DASCH_SCAMPIMAGES ${mosaic}${solutionString}_scamp.sexcat_filtered.fits ${mosaic}${solutionString}.refcat.ucac5.fits ${mosaic}${solutionString}_match_scamp.fits $match_radius $scale"
                    echo $cmd
                    $cmd

                    rmsfile=$DASCH_SCAMPIMAGES/${mosaic}${solutionString}_scamp.sexcat_filtered.fits.rms
                    if [ -e $rmsfile ] ; then
                        rmstext=$(cat $rmsfile)
                        nref_used=$(echo "$rmstext" |awk '{print $2}')
                        rms_ra=$(echo "$rmstext" |awk '{print $5}')
                        rms_dec=$(echo "$rmstext" |awk '{print $8}')

                        if [ $pass_two_flag -eq 0 ] ; then
                            # "original" mode: save SCAMP info to database
                            setScampRMS -p $mosaic -e $solutionNumber -n $nref_used -r $rms_ra -d $rms_dec
                            pass_two_error=0
                        else
                            # -h mode: "declare new algorithm SCAMP UCAC success providing that the new header is acceptable"
                            setScampRMS -p $mosaic -e $solutionNumber -h ${headeroutput} -n $nref_used -r $rms_ra -d $rms_dec
                            #echo "ERROR: early exit"
                            #exit
                            if [ $? -eq 0 ] ; then
                                pass_two_error=0
                                if [ -e ${headeroutput} ] ; then
                                    mv ${headeroutput} ${mosaic_fit_tnx}
                                fi

                                echo "UpdateMosaicTable -e $solutionNumber -m ${mosaic_fit_tnx}"
                                UpdateMosaicTable -e $solutionNumber -m ${mosaic_fit_tnx}
                            else
                                pass_two_error=1
                            fi
                        fi
                    else
                        pass_two_error=1
                    fi

                    #----------
                    # Get previous results with IRAF + mean per tile
                    #tnx_drad=$DASCH_SCRIPTS/match/${mosaic}${solutionString}_drad.db
                    tnx_drad=${DASCH_MATCH}/${mosaic}${solutionString}_drad.db
                    if [ -e $tnx_drad ] ; then
                        #echo
                        #echo ----------
                        #echo Get drad from previous pipeline
                        cmd="votable -i $tnx_drad -o $DASCH_MATCH/${mosaic}${solutionString}_match_tnx.xml"
                        echo $cmd
                        $cmd

                        cmd="$DASCH_SCRIPTS/tnxdrad_addcols.py $DASCH_MATCH/${mosaic}${solutionString}_match_tnx.xml"
                        echo $cmd
                        $cmd
                    fi

                    #----------
                    # Plot e.g. distance vs. offset on the plate
                    #echo
                    #echo ----------
                    #echo Plot distance vs. offset
                    cmd="$DASCH_SCRIPTS/plot_drad.py $DASCH_SCAMPIMAGES/${mosaic}${solutionString} $DASCH_MATCH/${mosaic}${solutionString}_match.fits $DASCH_MATCH/${mosaic}${solutionString}_match_scamp.fits $DASCH_MATCH/${mosaic}${solutionString}_match_tnx.fits $scale ${testflag}"
                    echo $cmd
                    $cmd
                fi
            fi
        else
            pass_two_error=1
            echo "ERROR: File does not exist: ${sex_cat}.head"
        fi
    else
        pass_two_error=1
        if [ $solutionNumber -eq 0 -o $useTranslateWCS -eq 0 ] ; then
            echo "ERROR: File does not exist: $DASCH_MATCH/${mosaic}${solutionString}_scamp.refcat_used_1.fits"
        fi
    fi

    #echo "ERROR: skip cleanup for ${mosaic}${solutionString}"
    #date
    #continue

    if [ $pass_two_error -eq 0 ] ; then
        if [ -e $mosaic_fit_tnx ] ; then
            newsize=$(ls -s --block-size=1 $mosaic_fit_tnx |awk '{print $1}')
            #echo "$mosaic_fit_tnx exists size $newsize"

            if [ -e $mosaic_fit ] ; then
                oldsize=$(ls -s --block-size=1 $mosaic_fit |awk '{print $1}')
                echo "comparefits $mosaic_fit $mosaic_fit_tnx"
                comparefits $mosaic_fit $mosaic_fit_tnx
                if [ $? -eq 0 ] ; then
                    echo "comparefits succeeded for $mosaic_fit"
                    #chmod 664 $mosaic_fit
                    #chmod 444 $mosaic_fit_tnx
                    #rm $mosaic_fit  # change of 2019-05-15
                    rm -f $DASCH_SCRATCH/${mosaic}_tnx.fit
                else
                    echo "ERROR $newsize - $oldsize is too small for $mosaic_fit_tnx"
                fi
            elif [ -e $DASCH_SCRATCH/${mosaic}_tnx.fit ] ; then
                oldsize=$(ls -s --block-size=1 $DASCH_SCRATCH/${mosaic}_tnx.fit |awk '{print $1}')
                echo "comparefits $DASCH_SCRATCH/${mosaic}_tnx.fit $mosaic_fit_tnx"
                comparefits $DASCH_SCRATCH/${mosaic}_tnx.fit $mosaic_fit_tnx
                if [ $? -eq 0 ] ; then
                    echo "comparefits succeeded for $DASCH_SCRATCH/${mosaic}_tnx.fit"
                    rm $DASCH_SCRATCH/${mosaic}_tnx.fit

                    if [ $pass_two_flag -eq 1 ] ; then
                        # Force an update of the Sextractor file
                        if [ -e $DASCH_MATCH/${mosaic}${solutionString}_tnx.db ] ; then
                            echo "mv $DASCH_MATCH/${mosaic}${solutionString}_tnx.db $DASCH_MATCH/${mosaic}${solutionString}.db"
                            mv $DASCH_MATCH/${mosaic}${solutionString}_tnx.db $DASCH_MATCH/${mosaic}${solutionString}.db
                        fi
                    fi
                else
                    echo "ERROR $newsize - $oldsize is too small for $mosaic_fit_tnx"
                fi
            fi
        fi
    fi

    # Clean up stray files
    if [ $pass_two_flag -eq 1 ] ; then
        #echo "ERROR: early exit"
        #exit
        scampobjectfile=${DASCH_ASTROMETRY}/${mosaic}_full_scamp.object
        scampreffile=${DASCH_ASTROMETRY}/${mosaic}_full_scamp.ref
        scampconffile=${DASCH_ASTROMETRY}/${mosaic}_full_scamp.conf

        rm -f $DASCH_SCRATCH/${mosaic}_tnx.fit
        rm -f $scampobjectfile
        rm -f $scampreffile
        rm -f $scampconffile

        if [ $solutionNumber -eq 0 ] ; then
            rm -f $headeroutput
        fi
    fi

    rm -f $DASCH_MATCH/${mosaic}${solutionString}_scamp.hdr
    rm -f $DASCH_MATCH/${mosaic}${solutionString}.tmp
    rm -f $DASCH_MATCH/${mosaic}${solutionString}.sexcat_filtered.fits
    rm -f $DASCH_SCAMPIMAGES/${mosaic}${solutionString}.sexcat_filtered.fits.rms
    rm -f $DASCH_MATCH/${mosaic}${solutionString}.refcat.ucac5.db
    rm -f $DASCH_MATCH/${mosaic}${solutionString}.refcat.ucac5.fits
    rm -f $DASCH_MATCH/${mosaic}${solutionString}.sexcat_filtered.head
    rm -f $DASCH_MATCH/${mosaic}${solutionString}.sexcat_filtered.head2
    rm -f $DASCH_MATCH/${mosaic}${solutionString}_scamp_results.xml
    rm -f $DASCH_MATCH/${mosaic}${solutionString}_scamp.refcat_used_1.fits
    rm -f $DASCH_MATCH/${mosaic}${solutionString}.sexcat_filtered.head_cor
    rm -f $DASCH_MATCH/${mosaic}${solutionString}.tmp2
    rm -f $DASCH_MATCH/${mosaic}${solutionString}_scamp.db
    rm -f $DASCH_MATCH/${mosaic}${solutionString}_scamp.sexcat_filtered.fits
    rm -f $DASCH_SCAMPIMAGES/${mosaic}${solutionString}_scamp.sexcat_filtered.fits.rms
    rm -f $DASCH_MATCH/${mosaic}${solutionString}_match.fits
    rm -f $DASCH_MATCH/${mosaic}${solutionString}_match_tnx.xml
    rm -f $DASCH_MATCH/${mosaic}${solutionString}_match_tnx.fits
    rm -f $DASCH_MATCH/${mosaic}${solutionString}_match_scamp.fits
done

# 2010-06-01 Mathieu Servillat Initial Version
# 2010-09-20 Edward J. Los     Integrate into the DASCH photometry pipeline.
#                                Change result name from_scamp to _tnx
#                                Use $DASCH_SCRIPTS for python scripts
#                                Use $DASCH_MATCH   for intermediate files
# 2011-01-31 Edward J. Los     For solutionNumber > 0, replace scamp with TranslateWCS
# 2011-06-07 Edward J. Los     Add comparefits check because of an issue with Odyssey's filesystem overallocating bytes for mosaics
# 2012-09-17 Edward J. Los     Modify for 2GB mosaics
# 2013-02-01 Edward J. Los     Switch to UCAC4
# 2013-03-04 Edward J. Los     If an InaccuratePV error occurs within scamp, then reduce the polynomial order
# 2013-07-23 Edward J. Los     Exit quietly if second pass and the plate is not a patrol plate
# 2014-04-15 Edward J. Los     Mark miscellaneous errors deferred with ERRORD
# 2014-06-17 Edward J. Los     Correct "ERRORD: run_scamp can not find mosaic" for a second pass.
# 2016-04-08 Edward J. Los     Change the mosaic protection to read only
# 2016-09-12 Edward J. Los     Back out protection change
# 2018-05-28 Edward J. Los     Switch to UCAC5
# 2018-07-30 Edward J. Los     Increase the match radius for patrol plates to 100
# 2019-05-15 Edward J. Los     Defer deletion of comparefits source
