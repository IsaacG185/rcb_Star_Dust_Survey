# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  display_thumbnail.csh  - display all thumbnails with ds9
#
#  May  2, 2008  Edward J. los
#
set imagelist = $1
set scripts = $DASCH_SCRIPTS
set binning = 16

echo "Processing $imagelist"
foreach plate ( `cat $imagelist` ) 
    set thumbnail = `echo "$plate.fit" | sed 's/_01\./_16\./g' | sed 's/_01w/_16w/g' | sed 's/_01r/_16r/g'`
    set directory = `$DASCH_SCRIPTS/getdirectory $plate`
    set filename = ${directory}/${thumbnail}

    if (-e $filename) then
    
    ds9 $filename

    else
        echo "ERROR: $filename not found"

    endif

end
