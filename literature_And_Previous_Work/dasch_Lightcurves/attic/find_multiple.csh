# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  find_multiple.csh - use David Hogg's astrometry.net to look for multiple exposures
#
#  NOTE: THIS FILE IS OBSOLETE
#
#  May 26, 2008 Edward J. Los - initial version
#  May 27, 2008 Edward J. Los - correct solvedFlag handling
#
#  NOTE: THIS FILE IS OBSOLETE
#
set imagelist = $1
set scripts = $DASCH_SCRIPTS
set binning = 1

set debugMode = 0
set fractionFLAG = "-f 0.01"
set kronFLAG = ""
set configSuffix = 1
date
cd ${DASCH_ASTROMETRY}

set listroot = $imagelist:t

echo "Processing ${DASCH_SCRIPTS}/$listroot"
foreach origplate ( `cat ${DASCH_SCRIPTS}/$listroot` ) 

  set solvedFlag = 0
  set plate = `echo  "$origplate" | sed 's/ww/wwe1/g'`

  if ($plate == $origplate) then
    echo "No WCS solution exists for $plate"
    continue

  endif

  set wcsfile = ${DASCH_ASTROMETRY}/${plate}.wcs
  if (-e $wcsfile) then
    echo "Solution already exists for $wcsfile"
    continue
  endif


  set allobjects = ${DASCH_INGEST}/${origplate}_allobjects.db

  if (!(-e $allobjects)) then
    echo "No allobjects file exists: $allobjects"
    continue

  endif

  set sextractorFile = ${DASCH_MATCH}/${origplate}_tnx.db
  if (!(-e $sextractorFile)) then
    echo "ERROR: No sextractor file found: $sextractorFile"
    continue
  endif

  set noRefFile = ${DASCH_ASTROMETRY}/${plate}_noref.db
  set sextractorCopy = ${DASCH_ASTROMETRY}/${plate}_tnx.db
  set joinFile = ${DASCH_ASTROMETRY}/${plate}_join.db


  row 'REF == "NONE"' < $allobjects | column NUMBER > $noRefFile
  cp $sextractorFile $sextractorCopy

  index -mb -n $sextractorCopy NUMBER
  index -mb -n $noRefFile  NUMBER
  jointable -j NUMBER -n $sextractorCopy $noRefFile > $joinFile 


  echo "Attempting solution for plate  $plate"


  set scaletext = `$DASCH_SCRIPTS/getseries $plate -f 1.1 -b $binning`
  set minscale = `echo $scaletext | gawk '{print $1}'`
  set maxscale = `echo $scaletext | gawk '{print $3}'`

    if (-e ${DASCH_ASTROMETRY}/fail_${plate}.txt) rm ${DASCH_ASTROMETRY}/fail_${plate}.txt
    if (-e ${DASCH_ASTROMETRY}/${plate}.db ) rm ${DASCH_ASTROMETRY}/${plate}.db 
  
    set xylsfile =  ${DASCH_ASTROMETRY}/${plate}.xyls
    set resultfile = ${DASCH_ASTROMETRY}/${plate}.solve
    if ($debugMode == 1) then
       set debugfile = ${DASCH_ASTROMETRY}/${plate}.debug2
       set debugflag = "-d $debugfile"
    else
       set debugflag = ""
    endif
    set mosaicsize = `$DASCH_SCRIPTS/getlocation -a $origplate`
    set width = `echo "$mosaicsize" | gawk '{print $1}'`
    set height  = `echo "$mosaicsize" | gawk '{print $2}'`


    foreach casenumber (1 2 3 4 5)
      if ($solvedFlag == 0) then
  
  
        if ($casenumber == 1) then
          set filterFLAG = "-e"
          set ellipticityFLAG = ""
          set maxdepth = 120
          set maxstars = 240
        endif
        if ($casenumber == 2) then
          set filterFLAG = "-e"
          set ellipticityFLAG = "-l 0.3"
          set maxdepth = 120
          set maxstars = 240
        endif
        if ($casenumber == 3) then
          set filterFLAG = ""
          set ellipticityFLAG = ""
          set maxdepth = 120
          set maxstars = 240
        endif
        if ($casenumber == 4) then
          set filterFLAG = "-e"
          set ellipticityFLAG = "-l 0.3"
          set maxdepth = 240
          set maxstars = 480
        endif
        if ($casenumber == 5) then
          set filterFLAG = "-e"
          set ellipticityFLAG = ""
          set maxdepth = 240
          set maxstars = 480
        endif
  
  
        echo "db2xyls -i $joinFile -o $xylsfile -w $width -h $height -n $maxstars $fractionFLAG $ellipticityFLAG  $kronFLAG $filterFLAG $debugflag"
        db2xyls -i $joinFile -o $xylsfile -w $width -h $height -n $maxstars $fractionFLAG $ellipticityFLAG  $kronFLAG $filterFLAG $debugflag 
        if ($debugMode == 1) then
           column -a -b < $debugfile  > ${DASCH_ASTROMETRY}/${plate}.debug
            rm $debugfile
        endif


        date
        echo "solve-field $xylsfile --no-tweak --scale-units arcsecperpix --scale-low $minscale --scale-high $maxscale --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_ISO --sort-ascending --overwrite --width $width --height $height  --no-fits2fits --depth 0-$maxdepth --no-plots >& $resultfile"
        solve-field $xylsfile --no-tweak --scale-units arcsecperpix --scale-low $minscale --scale-high $maxscale --x-column X_IMAGE --y-column Y_IMAGE --sort-column MAG_ISO --sort-ascending --overwrite --width $width --height $height  --no-fits2fits --depth 0-$maxdepth --no-plots  >& $resultfile
        date

        # FIXME: update these lines to adapt to different output format of newer Astrometry.Net code
        set centerra = `cat $resultfile | grep "Field center decimal" | awk '{print $6}'`
        set centerdec = `cat $resultfile | grep "Field center decimal" | awk '{print $7}'`
        set scalex = `cat $resultfile | grep "Plate scale" | awk '{print $3}'`
        set scaley = `cat $resultfile | grep "Plate scale" | awk '{print $5}'`
        if ($#centerdec == 0)  then
            cp $resultfile ${DASCH_ASTROMETRY}/fail_${plate}.txt
        else
            set solutiontext = `grep "Solution found"  $resultfile`
            echo "$solutiontext for case $casenumber"
            gethead ${DASCH_ASTROMETRY}/${plate}.wcs IMAGEW IMAGEH CTYPE1 CTYPE2 CRVAL1 CRVAL2 CRPIX1 CRPIX2 CD1_1 CD1_2 CD2_1 CD2_2 > ${DASCH_ASTROMETRY}/${plate}.gethead
            set naxis1 = `cat ${DASCH_ASTROMETRY}/${plate}.gethead | awk '{print $1}'`
            set naxis2 = `cat ${DASCH_ASTROMETRY}/${plate}.gethead | awk '{print $2}'`
            set ctype1 = `cat ${DASCH_ASTROMETRY}/${plate}.gethead | awk '{print $3}' | sed 's/\-SIP//g'`
            set ctype2 = `cat ${DASCH_ASTROMETRY}/${plate}.gethead | awk '{print $4}' | sed 's/\-SIP//g'`
            set crval1 = `cat ${DASCH_ASTROMETRY}/${plate}.gethead | awk '{print $5}'`
            set crval2 = `cat ${DASCH_ASTROMETRY}/${plate}.gethead | awk '{print $6}'`
            set crpix1 = `cat ${DASCH_ASTROMETRY}/${plate}.gethead | awk '{print $7}'`
            set crpix2 = `cat ${DASCH_ASTROMETRY}/${plate}.gethead | awk '{print $8}'`
            set cd1_1 = `cat ${DASCH_ASTROMETRY}/${plate}.gethead | awk '{print $9}'`
            set cd1_2 = `cat ${DASCH_ASTROMETRY}/${plate}.gethead | awk '{print $10}'`
            set cd2_1 = `cat ${DASCH_ASTROMETRY}/${plate}.gethead | awk '{print $11}'`
            set cd2_2 = `cat ${DASCH_ASTROMETRY}/${plate}.gethead | awk '{print $12}'`
            echo "plate ra dec scalex scaley binning naxis1 naxis2 ctype1 ctype2 crval1 crval2 crpix1 crpix2 cd1_1 cd1_2 cd2_1 cd2_2" >! ${DASCH_ASTROMETRY}/${plate}.tmp
            echo "----- -- --- ------ ------ ------- ------ ------ ------ ------ ------ ------ ------ ------ ----- ----- ----- -----" >> ${DASCH_ASTROMETRY}/${plate}.tmp
            echo "$plate $centerra $centerdec $scalex $scaley $binning $naxis1 $naxis2 $ctype1 $ctype2 $crval1 $crval2 $crpix1 $crpix2 $cd1_1 $cd1_2 $cd2_1 $cd2_2" >> ${DASCH_ASTROMETRY}/${plate}.tmp
            cat ${DASCH_ASTROMETRY}/${plate}.tmp | gawk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' >!  ${DASCH_ASTROMETRY}/${plate}.db  
            rm ${DASCH_ASTROMETRY}/${plate}.tmp
            rm ${DASCH_ASTROMETRY}/${plate}.gethead
            echo "Solved: ra $centerra dec $centerdec scale $scalex $scaley for $plate"

            set solvedFlag = 1
       endif


      # solvedFlag test
      endif
    # case loop
    end
  if ($solvedFlag == 0) then
     echo "No result for $plate"
  endif

  rm ${DASCH_ASTROMETRY}/${plate}.axy
  rm ${DASCH_ASTROMETRY}/${plate}.match
  rm ${DASCH_ASTROMETRY}/${plate}.rdls
  if (-e ${DASCH_ASTROMETRY}/${plate}.solved) rm ${DASCH_ASTROMETRY}/${plate}.solved 
  if (-e ${DASCH_ASTROMETRY}/${plate}-indx.xyls) rm ${DASCH_ASTROMETRY}/${plate}-indx.xyls
  rm $xylsfile


  rm $noRefFile
  rm ${noRefFile}.NUMBER-b
  rm ${sextractorCopy}.NUMBER-b
  rm $sextractorCopy
  rm $joinFile
end
