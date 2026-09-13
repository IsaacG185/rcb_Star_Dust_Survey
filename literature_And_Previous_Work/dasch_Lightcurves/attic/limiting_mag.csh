# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Given a list of plates, prepare a table of limiting magnitudes
#
#  Jun 18, 2008 Edward J. Los - Initial version
#
#
set local = `pwd`
set imagelist = $1
set totalcounter = 0
set goodcounter = 0
set missingcounter = 0

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
set dbfile = "${DASCH_INGEST}/limiting_mag_${imagelist:t:r}${qualifier}.db"

date
echo "Writing $dbfile"
if (-e $dbfile) then
  rm $dbfile
endif

set headerflag = ""
set dateval = `date`
echo "$dateval Processing ${imagelist}${qualifier}"
foreach name ( `cat $imagelist` ) 

  set filename = "${DASCH_INGEST}/${name}${qualifier}.out.spatial_bins.db"
  @ totalcounter = $totalcounter + 1
  if (-e $filename) then
    echo "Reading $filename"
    @ goodcounter = $goodcounter + 1
    set plate = `echo "$name" | gawk -F"_" '{print $1}'`

    column -i $filename limiting_mag spatial_bin | row 'spatial_bin == 1' | setcolumn Plate ${plate} | column $headerflag >> $dbfile
  set headerflag = "-b"
  else
    echo "No file found: $filename"
    @ missingcounter = $missingcounter + 1
  endif
end
cd ${DASCH_INGEST}
pwd
set dateval = `date`
echo "$dateval Completed"
