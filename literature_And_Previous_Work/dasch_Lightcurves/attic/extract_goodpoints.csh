# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Dec 13, 2007 Edward J. Los - Initial version
#
#  Given a list of plates, build a master database of all good datapoints
#  Then find out how many unique GSC2 objects we have and build a summary datbase
#  of these objects
#
#  Dec  5, 2007 Edward J. Los - Initial version
#  Dec 20, 2007 Edward J. Los - Reject stars if FLAGS >= 65536
#  Dec 29, 2007 Edward J. Los - Support recalibrated GSC magnitudes.
#                               Add counters.
#  Jan  6, 2008 Edward J. Los - Switch to using the output of the recover_points
#                               step.
#  Feb 19, 2008 Edward J. Los - Remove the julian date since recover_points.c now
#                               calculates it correctly
#  Feb 29, 2008 Edward J. Los - Add limiting magnitude, rms and bin to the master files
#  Mar  5, 2008 Edward J. Los - Add extinction to the master files
#  Mar  8, 2008 Edward J. Los - Add FLAGS and CAL_FLAG to the master file
#                               Include bin 9 entries in the master file.
#  Mar 10, 2008 Edward J. Los - Do not test for magcal_local_rms < 1
#                               Include npoints_local, magcal_iso_rms
#                               Output magcal_local_rms instead of magcal_local_error
#  Mar 24, 2008 Edward J. Los - Move counters into computmag
#  Apr  1, 2008 Edward J. Los - Include unmatched stars (computemag will not use them)
#                               Add Stdmag, dra, and ddec
#  Apr 22, 2008 Edward J. Los - add FLUX_ISO MAG_APER MAG_AUTO KRON_RADIUS BACKGROUND THRESHOLD 
#                               FLUX_MAX THETA_J2000 ELLIPTICITY ISOAREA_WORLD FWHM_IMAGE 
#                               FWHM_WORLD ISO0 ISO1 ISO2 ISO3 ISO4 ISO5 ISO6 ISO7 plate_dist
#  Apr 24, 2008 Edward J. Los - Now always use -t qualifier to throw away points within 0.5
#                               of limiting_mag_local. 
#                               Change median_local to clip_median_local
#                                      rms_local to clip_rms_local
#                                      ngood to clip_ngood
#  Jun 18, 2008 Edward J. Los - add magcal_iso
#
#  Jun 13, 2008 Edward J. Los - Add magcal_local_error
#  Jan 27, 2009 Edward J. Los - Split FLAGS into AFLAGS and BFLAGS, remove CAL_FLAG
#  Feb 24, 2009 Edward J. Los - Add local_bin_index to the master file
#  Mar 10, 2009 Edward J. Los - Add color
#  Jun 16, 2009 Edward J. Los - Add gsc_bin_index
set local = `pwd`
set imagelist = $1
set totalcounter = 0
set goodcounter = 0
set missingcounter = 0
# if mode is "_all", then we have V3.3 behavior where blended stars are included.
set mode = "_all"

if ($#argv == 1) then
  set qualifier = ""
else 
  if ($#argv == 2) then
      set qualifier = `echo "_$2"`
  else
     echo "usage: extract_goodpoints.csh list [qualifier]"
     echo "       list      is the list of mosaics to process"
     echo "       qualifier is a string appended to output filenames"
      exit
  endif
endif

set catalog = $DASCH_CATALOG
set scripts = $DASCH_SCRIPTS
set completedFlag = `echo "$DASCH_COMPLETED" | wc -c`
set dbfile = "${DASCH_INGEST}/master${mode}_${imagelist:t:r}${qualifier}.db"
set idfile = "${DASCH_INGEST}/id${mode}_${imagelist:t:r}${qualifier}.db"
set sumfile = "${DASCH_INGEST}/sum${mode}_${imagelist:t:r}${qualifier}.db"

date
echo "Writing $dbfile"
if (-e $dbfile) then
  rm $dbfile
endif

set headerflag = ""
set dateval = `date`
echo "$dateval Processing ${imagelist}${qualifier}"
foreach name ( `cat $imagelist` ) 

  set filename = "${DASCH_INGEST}/${name}${qualifier}_allobjects.db"
  @ totalcounter = $totalcounter + 1
  if (-e $filename) then
    echo "Reading $filename"
    @ goodcounter = $goodcounter + 1
    set plate = `echo "$name" | gawk -F"_" '{print $1}'`

    column -i $filename REF NUMBER ra dec magcal_local Date extinction AFLAGS BFLAGS limiting_mag_local magcal_local_error spatial_bin npoints_local magcal_iso magcal_iso_rms magcal_local_rms Stdmag color dra ddec FLUX_ISO MAG_APER MAG_AUTO KRON_RADIUS BACKGROUND THRESHOLD FLUX_MAX THETA_J2000 ELLIPTICITY ISOAREA_WORLD FWHM_IMAGE FWHM_WORLD ISO0 ISO1 ISO2 ISO3 ISO4 ISO5 ISO6 ISO7 plate_dist local_bin_index gsc_bin_index | setcolumn Plate ${plate} | column $headerflag >> $dbfile
  set headerflag = "-b"
  else
    echo "No file found: $filename"
    @ missingcounter = $missingcounter + 1
  endif
end
cd ${DASCH_INGEST}
pwd
index -mb $dbfile REF
#set dateval = `date`
#echo "$dateval Creating index"
#echo "column -b -i ${dbfile}  REF > ${idfile}.tmp"
#column -b -i ${dbfile}  REF > ${idfile}.tmp
#set dateval = `date`
#echo "$dateval Sorting index"
#sort -u ${idfile}.tmp > ${idfile}.tmp1
#rm ${idfile}.tmp
#wc -l $dbfile 
#wc -l ${idfile}.tmp1
echo "Total Plates: $totalcounter Found: $goodcounter Missing: $missingcounter"
set dateval = `date`
echo "$dateval Executing $DASCH_SCRIPTS/computemag -t -n 10 -c -b -r  ${dbfile}"

$DASCH_SCRIPTS/computemag -t -n 10  -c -b -r ${dbfile}

set dateval = `date`
echo "$dateval Creating ${sumfile}"
#rm ${idfile}.tmp1
#rm ${idfile}.tmp2
#rm ${idfile}.tmp3
#rm ${idfile}.tmp4
echo "minmag maxmag count median" >! ${sumfile}.tmp
echo "------ ------ ----- ---" >> ${sumfile}.tmp
foreach minmag (18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3) 
  @ maxmag = $minmag + 1
  set count  = `row -i $idfile "(clip_median_local >= $minmag) && (clip_median_local < $maxmag) && (clip_ngood >= 10)" | column -b clip_rms_local | wc -l`
  set median = `row -i $idfile "(clip_median_local >= $minmag) && (clip_median_local < $maxmag) && (clip_ngood >= 10)" | column clip_rms_local | statstable | grep Median | gawk '{print $2}'`


  echo "$minmag $maxmag $count $median"
  echo "$minmag $maxmag $count $median" >> ${sumfile}.tmp
 
end
set count  = `row -i $idfile "(clip_ngood >= 10)" | column -b clip_rms_local | wc -l`
set median = `row -i $idfile "(clip_ngood >= 10)" | column clip_rms_local | statstable | grep Median | gawk '{print $2}'`
echo "3 19 $count $median"
echo "3 19 $count $median" >> ${sumfile}.tmp
cat ${sumfile}.tmp | sed 's/ /\t/g' > ${sumfile}
rm  ${sumfile}.tmp
index -mb ${idfile} REF
set dateval = `date`
echo "$dateval Completed"
