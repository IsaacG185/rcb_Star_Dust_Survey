# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Given a list of plates, a center position and a radius, extract images from the fits files.
# 
#  extractimage.csh /dasch/Pipeline/dell.list 08:40:24.000 +19:41:00.00 300 /dasch/junk/plates
#
# Apr 20, 2010 Edward J. Los - Initial implementation


if ($#argv != 5) then
    echo "usage: extractimage list ra dec radius directory"
    echo "       list      is the list of mosaics process"
    echo "       ra        is the right ascension"
    echo "       dec       is the declination"
    echo "       radius    is the radius in arcsec"
    echo "       directory is the target directory for the extracted images"
    exit
endif

set imagelist = $1
set rastring = $2
set decstring = $3
set arcsec = $4
set targetdir = $5

if (! -e $targetdir) then
    echo "ERROR: $targetdir does not exist"
    exit
endif






date
echo "Processing $imagelist"


foreach plate ( `cat $imagelist` ) 
    #echo "Extracting from plate $plate"
    set directory = `$DASCH_SCRIPTS/getdirectory $plate`
    set filename = ${directory}/${plate}_tnx.fit
    set targetname = ${targetdir}/${plate}_extract.fit
    set location =  `$DASCH_SCRIPTS/getlocation $plate -t $arcsec`
    set pixels = `echo "$location" | gawk '{print $1}'`
 
 

    if (-e $filename) then
       #echo "Found $filename"

 
       echo "getfits -i 0 -o $targetname -x $pixels $pixels $filename $rastring $decstring J2000"
             getfits -i 0 -o $targetname -x $pixels $pixels $filename $rastring $decstring J2000

    else
      set tmpflag = `echo  "$plate" | gawk '/_s/{print $1}'`
      #echo "TMPFLAG $tmpflag ${#tmpflag}"

      if (${#tmpflag} < 1) then
         echo "ERROR: $filename does not exist"
      endif
    endif

end
date 
