# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License
#
# 2010-06-01 Mathieu Servillat
#
# Astrometry correction of DASCH mosaics with Scamp
#

if ($#argv != 2) then
    echo "usage: run_scamp.csh  solutionNumber list"
    exit
endif

set solutionNumber = $1
set imagelist = $2
set scripts = $DASCH_SCRIPTS
set plot_drad = 1

set listfile = $imagelist:t
echo "Processing ${DASCH_SCRIPTS}/$listfile"

foreach mosaic ( `cat ${DASCH_SCRIPTS}/$listfile` )

    set mosaic_dir = `$DASCH_SCRIPTS/getdirectory $mosaic`
    set mosaic_fit = ${mosaic_dir}/${mosaic}.fit
    set mosaic_fit_tnx = ${mosaic_dir}/${mosaic}_tnx.fit
    set mosaic_fit_tnx_save = ${mosaic_dir}/${mosaic}_tnx_save.fit
    set mosaic_fit_scamp = ${mosaic_dir}/${mosaic}_scamp.fit

    set TSTART = `date`
    echo
    echo Begin Process $TSTART
    echo
    echo mosaic = $mosaic
    echo

    #----------
    # Set parameters

    set script_output = `$DASCH_SCRIPTS/getlocation  $mosaic -e $solutionNumber`
    set scale  = `echo "$script_output" | gawk '{print $4}'`
    if ($scale == "") then
	#echo "ERROR getlocation failed with args $plate -e $solutionNumber"
	continue
    endif
    set aper_arcsec = 10
    set aper_pix = `echo $aper_arcsec $scale | gawk '{print $1/$2}'`

    # Global parameters
    set detthresh = 5
    set maglim = 8.0,99.0
    set snthresh = 10.0,100.0
    set astrclip = 3.0
    set ndet_max = 10000
    set nref_max = 10000 # will be determined as ndet/2

    # Specific parameters
    #set match_radius = 10 # mc
    #set match_radius = 210 # ac
    echo "set: $DASCH_SCRIPTS/scamp/get_match_radius.py $mosaic"
    set match_radius = `$DASCH_SCRIPTS/scamp/get_match_radius.py $mosaic`

    echo "Pixel Scale: " $scale
    echo "Max distortion expected: " $match_radius

    #----------
    # Move to directory where all scamp images will be kept


    if ($solutionNumber == 0) then

	#----------
	#processing of the first exposure

	#----------
	# Recreate original mosaic if needed
	echo
	echo ----------
	if (! -e $mosaic_fit) then
	    echo "Recreate original mosaic"
	    if (-e $mosaic_fit_tnx) then
		if (! -e $mosaic_fit_tnx_save) then
		    echo "TNX mosaic found: $mosaic_fit_tnx"
		    echo "Backup TNX mosaic: $mosaic_fit_tnx_save"
		    cp $mosaic_fit_tnx $mosaic_fit_tnx_save
		endif
		echo "Prepare initial mosaic: $mosaic_fit"
		echo $mosaic >! temp.list
		$DASCH_SCRIPTS/prepare_mosaic.csh 0 temp.list
		rm temp.list
	    else
		echo "ERROR: No mosaic found in ${mosaic_dir}"
		stop
	    endif
	else
	    echo "Mosaic found: $mosaic_fit" 
	    if (-e $mosaic_fit_tnx_save) then
		echo "TNX mosaic backup found: $mosaic_fit_tnx_save"
	    else
		echo "WARNING: No TNX mosaic backup found"
	    endif
	endif

	#----------
	# Run Sextractor on mosaic
	echo
	echo ----------
	echo Run Sextractor
	set sex_options = "-c $DASCH_SCRIPTS/Sextractor/DASCH2_scamp.config -DETECT_THRESH $detthresh -ANALYSIS_THRESH $detthresh -PHOT_APERTURES $aper_pix -PARAMETERS_NAME $DASCH_SCRIPTS/Sextractor/DASCH_scamp.param -FILTER_NAME $DASCH_SCRIPTS/Sextractor/default.conv"
	echo sex $sex_options -CATALOG_NAME $DASCH_SCAMP/${mosaic}.cat.fits ${mosaic_fit}
	echo
	sex $sex_options -CATALOG_NAME $DASCH_SCAMP/${mosaic}.cat.fits ${mosaic_fit}

	#----------
	# Filter sextractor results
	echo "$DASCH_SCRIPTS/scamp/sexcat_filter.py $DASCH_SCAMP/${mosaic}.cat.fits $ndet_max"
	$DASCH_SCRIPTS/scamp/sexcat_filter.py $DASCH_SCAMP/${mosaic}.cat.fits $ndet_max

	#----------
	# Determine nref from ndet: nref = ndet/2
	echo "set: $DASCH_SCRIPTS/scamp/sexcat_ndet.py $DASCH_SCAMP/${mosaic}.cat.fits"
	set ndet = `$DASCH_SCRIPTS/scamp/sexcat_ndet.py $DASCH_SCAMP/${mosaic}.cat.fits`
	set nref_max = `echo $ndet 2 | gawk '{print $1/$2}'`

	#----------
	# Reference catalog for the field
	echo
	echo ----------
	set ref_cat = $DASCH_SCAMP/${mosaic}_ucac3
	echo Prepare reference catalog from UCAC-3
	echo ${ref_cat}.cat.fits
	echo
	echo "set: /home/scanner/junk/wcstools-3.8.1/bin/imcat -f -c ucac3 $mosaic_fit "
	set nref = `/home/scanner/junk/wcstools-3.8.1/bin/imcat -f -c ucac3 $mosaic_fit | gawk '{print $2}'`
	echo $nref reference stars available, keep at max $nref_max
	echo "/home/scanner/junk/wcstools-3.8.1/bin/imcat -dh -c ucac3 -n $nref_max $mosaic_fit > ${ref_cat}.db"
	/home/scanner/junk/wcstools-3.8.1/bin/imcat -dh -c ucac3 -n $nref_max $mosaic_fit > ${ref_cat}.db

	echo "$DASCH_SCRIPTS/scamp/refcat2fit.py ${ref_cat}.db ${ref_cat}.cat.fits"
	$DASCH_SCRIPTS/scamp/refcat2fit.py ${ref_cat}.db ${ref_cat}.cat.fits

	#----------
	# Run Scamp
	echo
	echo ----------
	set checkplot_type = FGROUPS,DISTORTION,ASTR_REFERROR2D,ASTR_REFERROR1D
	set checkplot_name = ${mosaic}_astr_fgroups,${mosaic}_astr_distort,${mosaic}_astr_referror2d,${mosaic}_astr_referror1d
	set scamp_options = "-c $DASCH_SCRIPTS/scamp/scamp.conf -ASTREF_CATALOG FILE -ASTREFCAT_NAME ${ref_cat}.cat.fits -CROSSID_RADIUS $match_radius -ASTREFMAG_LIMITS $maglim -SN_THRESHOLDS $snthresh -ASTRCLIP_NSIGMA $astrclip -CHECKPLOT_TYPE $checkplot_type -CHECKPLOT_NAME $checkplot_name -XML_NAME $DASCH_SCAMP/${mosaic}_scamp_results.xml"
	echo scamp $scamp_options $DASCH_SCAMP/${mosaic}.cat.fits
	echo
	scamp $scamp_options $DASCH_SCAMP/${mosaic}.cat.fits

	#----------
	# Create new header and new mosaic
	echo
	echo ----------
	echo Create new header and mosaic
	echo $DASCH_HEADERS/${mosaic}_scamp.hdr
	echo ${mosaic_fit_scamp}
	echo
	echo "$DASCH_SCRIPTS/scamp/create_mosaic_scamp.py ${mosaic_fit} $DASCH_SCAMP/${mosaic}.cat.head"
	$DASCH_SCRIPTS/scamp/create_mosaic_scamp.py ${mosaic_fit} $DASCH_SCAMP/${mosaic}.cat.head
	cp -f $DASCH_SCAMP/${mosaic}_scamp.hdr $DASCH_HEADERS

	#----------
	if ($plot_drad == 1) then
	    # Test: rerun sextractor on mosaic and match with reference catalog
	    # would be better to get new ra,dec form x,y of the first sextractor run
	    echo
	    echo ----------
	    echo sex $sex_options -CATALOG_NAME $DASCH_SCAMP/${mosaic}_scamp.cat.fits ${mosaic_fit_scamp}
	    echo
	    sex $sex_options -CATALOG_NAME $DASCH_SCAMP/${mosaic}_scamp.cat.fits ${mosaic_fit_scamp}

	    #----------
	    # Filter sextractor results
	    echo "$DASCH_SCRIPTS/scamp/sexcat_filter.py $DASCH_SCAMP/${mosaic}_scamp.cat.fits $ndet_max"
	    $DASCH_SCRIPTS/scamp/sexcat_filter.py $DASCH_SCAMP/${mosaic}_scamp.cat.fits $ndet_max

	    #----------
	    # Match with reference catalog
	    echo
	    echo ----------
	    echo Match detections with reference catalog
	    echo
	    #$DASCH_SCRIPTS/scamp/match_cat.py $DASCH_SCAMP/${mosaic}.cat.fits $DASCH_SCAMP/${mosaic}_scamp.cat.fits $DASCH_SCRIPTS/scamp/${mosaic}_match_shifts.fits $match_radius
	    #echo
	   
	    echo "$DASCH_SCRIPTS/scamp/match_cat_ref.py $DASCH_SCAMP/${mosaic}.cat.fits ${ref_cat}.cat.fits $DASCH_SCAMP/${mosaic}_match.fits $match_radius"
	    $DASCH_SCRIPTS/scamp/match_cat_ref.py $DASCH_SCAMP/${mosaic}.cat.fits ${ref_cat}.cat.fits $DASCH_SCAMP/${mosaic}_match.fits $match_radius
	    echo
	    echo "$DASCH_SCRIPTS/scamp/match_cat_ref.py $DASCH_SCAMP/${mosaic}_scamp.cat.fits ${ref_cat}.cat.fits $DASCH_SCAMP/${mosaic}_match_scamp.fits $match_radius"
	    $DASCH_SCRIPTS/scamp/match_cat_ref.py $DASCH_SCAMP/${mosaic}_scamp.cat.fits ${ref_cat}.cat.fits $DASCH_SCAMP/${mosaic}_match_scamp.fits $match_radius

	    #----------
	    # Get previous results with IRAF + mean per tile
	    #set tnx_drad = $DASCH_SCRIPTS/match/${mosaic}_drad.db
	    set tnx_drad = /home/scanner/Pipeline/match/${mosaic}_drad.db
	    if (-e $tnx_drad) then
		echo
		echo ----------
		echo Get drad from previous pipeline
		echo
		echo "$DASCH_SCRIPTS/votable -i $tnx_drad -o $DASCH_SCAMP/${mosaic}_match_tnx.xml"
		$DASCH_SCRIPTS/votable -i $tnx_drad -o $DASCH_SCAMP/${mosaic}_match_tnx.xml
		echo "$DASCH_SCRIPTS/scamp/tnxdrad_addcols.py $DASCH_SCAMP/${mosaic}_match_tnx.xml"
		$DASCH_SCRIPTS/scamp/tnxdrad_addcols.py $DASCH_SCAMP/${mosaic}_match_tnx.xml
	    endif

	    #----------
	    # Plot e.g. distance vs. offset on the plate
	    echo
	    echo ----------
	    echo Plot distance vs. offset
	    echo "$DASCH_SCRIPTS/scamp/plot_drad.py $DASCH_SCAMP/${mosaic} $DASCH_SCAMP/${mosaic}_match.fits $DASCH_SCAMP/${mosaic}_match_scamp.fits $DASCH_SCAMP/${mosaic}_match_tnx.fits"
	    $DASCH_SCRIPTS/scamp/plot_drad.py $DASCH_SCAMP/${mosaic} $DASCH_SCAMP/${mosaic}_match.fits $DASCH_SCAMP/${mosaic}_match_scamp.fits $DASCH_SCAMP/${mosaic}_match_tnx.fits
	endif

	#----------
	# Come back to TNX mosaic
	echo
	echo ----------
	echo "Delete scamp mosaic: $mosaic_fit_scamp"
	#rm $mosaic_fit_scamp
	if (! -e $mosaic_fit_tnx_save) then
	    echo "No TNX mosaic saved, keep original"
	else
	    echo "Recover TNX mosaic from backup"
	    echo "Delete original: $mosaic_fit"
	    #rm $mosaic_fit
	endif
	echo

    else

	#processing of additional exposures
	@ prevNumber = $solutionNumber - 1
	set solutionString = "_s$solutionNumber"
	if ($prevNumber == 0) then
	    set prevString = ""
	else
	    set prevString = "_s$prevNumber"
	endif

	set update_source =  ${DASCH_MATCH}/${plate}${prevString}.db
	set table         =  ${DASCH_MATCH}/${plate}${solutionString}.db
	set fits_image    =  ${DASCH_HEADERS}/${plate}${solutionString}.hdr
	set none_image    =  ${DASCH_ASTROMETRY}/${plate}${solutionString}_none.db

	echo "Additional exposures process not implemented yet"

    endif

end
