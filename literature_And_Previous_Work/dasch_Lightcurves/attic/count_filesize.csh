# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Feb 11, 2008 Edward J. Los - Initial version
#
#   This script counts miscellaneous file sizes for performance
#   characterization
#

set imagelist = $1
set catalog = $DASCH_CATALOG
set scripts = $DASCH_SCRIPTS
set completedFlag = `echo "$DASCH_COMPLETED" | wc -c`
set catdirectory = $catalog:h
echo "Processing $imagelist"
foreach plate ( `cat $imagelist` ) 
  echo $plate
#  set name = `echo $plate | sed 's/_tnx.fits//g'`
  set name = $plate
    set directory = `$DASCH_SCRIPTS/getdirectory $plate`

  date
    
# Sextractor file size
#    set sextr_filename = ${directory}/${plate}.db
#    set filesize = `wc -l $sextr_filename | gawk '{print $1}'`
#    echo "Sextractor $filesize for $plate"

# GSC2.3.2 catalog file size

  set fits_image = ${directory}/${plate}.fit
  set cal = ${catdirectory}/${plate}_cat.db
  set platelocation = `$DASCH_SCRIPTS/getlocation -g $fits_image`
#  set platelocation = "$platelocation -m2 17"
  set platelocation = "$platelocation -m2 18.5"
  echo "$platelocation > $cal"
  $platelocation > $cal
  set catsize = `wc -l $cal | gawk '{print $1}'`
  echo "Catalog size is $catsize for $plate"
  rm $cal
end
