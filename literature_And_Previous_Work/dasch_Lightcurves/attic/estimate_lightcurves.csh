#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  estimate_lightcurves.csh
#
#   Prepare master file for use with estimate_lightcurves.c
#
#
#  Apr 14, 2012 Edward J. Los - Initial version
#  May 14, 2012 Edward J. Los - Correct for Julian Date prefixes
#                               Specify the working directory
#                               Correct a join table bug
#
#  Example: ./estimate_lightcurves.csh platelist.xxx /home/scanner/junk/epsilondasch
#
set imagelist = $1
set directory = $2
set local = $directory

if ($#argv == 2) then
  set qualifier = ""
  set suffix = "gsc"
  set qualcmd = ""
else 
  if ($#argv == 3) then
      set qualifier = `echo "_$3"`
      set suffix = "$3"
      set qualcmd =   `echo "-q $3"`
      
  else
     echo "usage: estimate_lightcurves.csh list directory [qualifier]"
     echo "       directory is the working directory"
     echo "       list      is the list of mosaics to process of the form mosaicname%rootname"
     echo "       qualifier is a string appended to output filenames"
      exit
  endif
endif

set aper_arcsec = 10
set scripts = $DASCH_SCRIPTS

date
set dbfile = ${local}/master.db
echo "New master file $dbfile"
if (-e $dbfile) then
  rm $dbfile
endif
set magisofile = ${local}/magiso.db
echo "New mag_iso file $magisofile"
if (-e $magisofile) then
  rm $magisofile
endif

set sequence = ${local}/sequence.db
echo "New sequence file: $sequence"
if (-e $sequence) then
  rm $sequence
endif

set listroot = $imagelist:t
set counter = 0
set counter2 = 0
cd $directory
pwd
echo "Processing ${local}/$listroot"
foreach plateentry ( `cat ${local}/$listroot` ) 
    set plate = `echo "$plateentry" | gawk -F% '{print $1}'`
    set root =  `echo "$plateentry" | gawk -F% '{print $2}'`
    set fitsfile = ${local}/${plate}.fit
    set catfile  = ${local}/${root}.cat
    #echo "plate $plate root $root fitsfile $fitsfile"
    set script_output = `$DASCH_SCRIPTS/getlocation -p ${root}_01`
    set scale  = `echo "$script_output" | gawk '{print $4}'` 
    set aper_pix = `echo $aper_arcsec $scale | gawk '{print $1/$2}'`
    set juliandate = `$DASCH_SCRIPTS/juliandate ${root}_01 | gawk '{print $1}'`
    set width = `listhead $fitsfile | grep "NAXIS1" | gawk '{print $3}'`
    set height = `listhead $fitsfile | grep "NAXIS2" | gawk '{print $3}'`
    set pixelRadius = "-r 9"
    set series = `$DASCH_SCRIPTS/getlocation -p -s ${root}_01 | gawk '{print $1}'`
    if (($series == "ac") || ($series == "am") || ($series == "ax") || ($series == "ay") || ($series == "ca") || ($series == "x")) then
      set set pixelRadius = "-r 21"
    endif
    set script_output = `$DASCH_SCRIPTS/getlocation -p  ${root}_01 $pixelRadius`
    set match_radius = `echo "$script_output" | gawk '{print $5}'`
    set plate_scale = `echo "$script_output" | gawk '{print $4}'`
    #echo "aper_arcsec $aper_arcsec scale $scale aper_pix $aper_pix script_output $script_output"

    echo "Executing: sex -c ${scripts}/Sextractor/DASCH.config $fitsfile -PARAMETERS_NAME ${scripts}/Sextractor/DASCH.param -FILTER_NAME ${scripts}/Sextractor/default.conv -PHOT_APERTURES $aper_pix  -CATALOG_NAME $catfile"
    sex -c ${scripts}/Sextractor/DASCH.config $fitsfile -PARAMETERS_NAME ${scripts}/Sextractor/DASCH.param -FILTER_NAME ${scripts}/Sextractor/default.conv -PHOT_APERTURES $aper_pix  -CATALOG_NAME $catfile

    sextotable < ${local}/${root}.cat | sed 's/FLAGS/BFLAGS/g' | sed 's/ALPHA_J2000/ra/g' | sed 's/DELTA_J2000/dec/g' > ${local}/${root}tmp.db
     column -i ${local}/${root}tmp.db -a plate_dra plate_ddec plate_dist | compute "plate_ddec=0; plate_dra=0; plate_dist=0" > ${local}/${root}.db
    echo "update_sextractor -v -p -u -e 0  -l 0 -r 0 -b 0 -t 0 -i ${local}/${root}.db -o ${local}/${root}ww.db  -f ${local}/${root}_magiso.db -m  ${local}/${plate}.fit"
    update_sextractor -v -p -u -e 0  -l 0 -r 0 -b 0 -t 0 -i ${local}/${root}.db -o ${local}/${root}ww.db -f ${local}/${root}_magiso.db -m  ${local}/${plate}.fit

    #echo "ERROR: early exit"
    #exit

		if (-e  ${local}/${root}_magiso.db) then
       @ counter = $counter + 1
       if ($counter == 1) then
           mv ${local}/${root}_magiso.db $magisofile
       else 
				   column -a -b -i ${local}/${root}_magiso.db > ${local}/${root}_magiso.tmp
           mv $magisofile ${magisofile}.tmp
           cat ${magisofile}.tmp ${local}/${root}_magiso.tmp > $magisofile
					 rm ${magisofile}.tmp
					 rm ${local}/${root}_magiso.tmp
			     rm  ${local}/${root}_magiso.db
       endif
    endif

    index -mb -n ${local}/${root}ww.db dec


    echo "matchstars $qualcmd -j $juliandate -w $width -h $height -r $match_radius -i ${local}/${root}ww.db -e ${local}/match_${root}_${suffix}_lim.db -o ${local}/match_${root}_${suffix}.db -s  $plate_scale"
    matchstars $qualcmd -j $juliandate -w $width -h $height -r $match_radius -i ${local}/${root}ww.db -e ${local}/match_${root}_${suffix}_lim.db -o ${local}/match_${root}_${suffix}.db -s $plate_scale 

    echo "sorttable -i ${local}/match_${root}_${suffix}.db -n Y_IMAGE > ${local}/match_${root}_${suffix}.dbtmp"
    sorttable -i ${local}/match_${root}_${suffix}.db -n Y_IMAGE > ${local}/match_${root}_${suffix}.dbtmp

    echo "$DASCH_SCRIPTS/filterblended  $qualcmd -p -t 7200 -s $plate_scale -w $width -h $height -e ${local}/match_${root}_${suffix}_lim.db  -r $match_radius ${local}/match_${root}_${suffix}.dbtmp"
    $DASCH_SCRIPTS/filterblended  $qualcmd -p -t 7200 -s $plate_scale -w $width -h $height -e ${local}/match_${root}_${suffix}_lim.db  -r $match_radius ${local}/match_${root}_${suffix}.dbtmp

    echo "column -i ${local}/match_${root}_${suffix}_u.db ra_2 dec_2 Stdmag REF REF | sorttable -n Stdmag | column -a src_name | compute 'src_name = REF' | sed 's/_2//g' >> ${sequence}"
    @ counter2 = $counter2 + 1


    if ($counter2 == 1) then
       column -i ${local}/match_${root}_${suffix}_u.db ra_2 dec_2 Stdmag REF REF | sorttable -n Stdmag | column -a src_name | compute 'src_name = REF' | sed 's/_2//g' >> ${sequence}
    else
       column -i ${local}/match_${root}_${suffix}_u.db ra_2 dec_2 Stdmag REF REF | sorttable -n Stdmag | column -a src_name | compute 'src_name = REF' | sed 's/_2//g' | column -a -b  >> ${sequence}

    endif

    sorttable -n NUMBER < ${root}ww.db > ${root}ww.dbtmp1
    sorttable -n NUMBER < match_${root}_${suffix}_u.db > ${root}ww.dbtmp2
    jointable -n -j NUMBER ${root}ww.dbtmp1 ${root}ww.dbtmp2 > ${root}ww.dbtmp3
    set ddd = Plate=\""${root}"\"
    set eee = Date=\"$juliandate\"

    if ($counter2 == 1) then
        column -i  ${root}ww.dbtmp3 FLUX_ISO AFLAGS_2 BFLAGS_2 REF gsc_bin_index ra_2 dec_2 Stdmag FLUX_MAX THRESHOLD | column -a AFLAGS BFLAGS Date spatial_bin Plate ra dec | compute "AFLAGS=AFLAGS_2;BFLAGS=BFLAGS_2;spatial_bin=0;ra=ra_2;dec=dec_2" | compute $ddd | compute $eee >> $dbfile
    else 
        column -i  ${root}ww.dbtmp3 FLUX_ISO AFLAGS_2 BFLAGS_2 REF gsc_bin_index ra_2 dec_2 Stdmag FLUX_MAX THRESHOLD | column -a AFLAGS BFLAGS Date spatial_bin Plate ra dec | compute "AFLAGS=AFLAGS_2;BFLAGS=BFLAGS_2;spatial_bin=0;ra=ra_2;dec=dec_2" | compute $ddd | compute $eee | column -a -b >> $dbfile

    endif

    #echo "ERROR: Early exit"
    #exit

    if (-e ${local}/match_${root}_${suffix}_b.db) then
      rm ${local}/match_${root}_${suffix}_b.db
    endif
    if (-e ${local}/match_${root}_${suffix}.db) then
      rm ${local}/match_${root}_${suffix}.db
    endif
    if (-e ${local}/match_${root}_${suffix}.dbtmp) then
      rm ${local}/match_${root}_${suffix}.dbtmp
    endif
    if (-e ${local}/match_${root}_${suffix}_lim1.db) then
      rm ${local}/match_${root}_${suffix}_lim1.db
    endif
    if (-e ${local}/match_${root}_${suffix}_lim2.db) then
      rm ${local}/match_${root}_${suffix}_lim2.db
    endif
    if (-e ${local}/match_${root}_${suffix}_lim3.db) then
      rm ${local}/match_${root}_${suffix}_lim3.db
    endif
    if (-e ${local}/match_${root}_${suffix}_lim.db) then
      rm ${local}/match_${root}_${suffix}_lim.db
    endif
    if (-e ${local}/match_${root}_${suffix}_u.db) then
      rm ${local}/match_${root}_${suffix}_u.db
    endif
    if (-e ${local}/${root}.cat) then
      rm ${local}/${root}.cat
    endif
    if (-e ${local}/${root}.db) then
      rm ${local}/${root}.db
    endif
    if (-e ${local}/${root}tmp.db) then
      rm ${local}/${root}tmp.db
    endif
    if (-e ${local}/${root}_${suffix}_drad.db) then
      rm ${local}/${root}_${suffix}_drad.db
    endif
    if (-e ${local}/${root}ww.db) then
      rm ${local}/${root}ww.db
    endif
    if (-e ${local}/${root}ww.db.dec.b) then
      rm ${local}/${root}ww.db.dec.b
    endif
    if (-e ${local}/${root}ww.dbtmp1) then
      rm ${local}/${root}ww.dbtmp1
    endif
    if (-e ${local}/${root}ww.dbtmp2) then
      rm ${local}/${root}ww.dbtmp2
    endif
    if (-e ${local}/${root}ww.dbtmp3) then
      rm ${local}/${root}ww.dbtmp3
    endif



end
date
