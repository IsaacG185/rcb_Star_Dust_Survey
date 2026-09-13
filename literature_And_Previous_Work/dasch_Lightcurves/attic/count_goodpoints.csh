# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
# 
#  This script counts the number of good stars found in the list of images.  "good" is
#   magcal_iso < limiting_mag
#
# echo "count_goodpoints.csh good_scans.debug" | at now
#
#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead


set imagelist = $1

set ingestdir = $DASCH_INGEST

set totalplates = 0
set totalstars  = 0

cd $ingestdir
pwd
date
echo "Processing $imagelist"
foreach listentry ( `cat ${imagelist}` ) 
  set db = "${listentry}.out.local.db"

  if (-e $db) then
   date
   # First select the good stars only. These will be used to define the calibration map and will be corrected.
   # We choose only those stars brighter than their local limiting magnitude. Other selections could aso be used.
   set curcount = `row -i $db 'magcal_iso<limiting_mag' | wc -l`
   @ curcount = $curcount - 2
   @ totalplates = $totalplates + 1
   @ totalstars = $totalstars + $curcount
   echo "curcount: $curcount totalplates: $totalplates totalstars: $totalstars Database $db"
   else
    echo "Database $db does not exist"
  endif
end
echo "Total Plates $totalplates, total stars $totalstars"
date
