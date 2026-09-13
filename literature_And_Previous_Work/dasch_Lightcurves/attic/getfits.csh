# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  getfits_3c273.csh - Extract an object from mosaics
#
#  Oct 11, 2008 Edward J. Los - Initial version
#
set imagelist = $1
set scripts = $DASCH_SCRIPTS
date

set listroot = $imagelist:t

pwd
echo "Processing ${DASCH_SCRIPTS}/$listroot"
foreach plate ( `cat ${DASCH_SCRIPTS}/$listroot` ) 
    set directory = `$DASCH_SCRIPTS/getdirectory $plate`  
    set input_image = ${directory}/${plate}_tnx.fit
    set output_image = ${directory}/${plate}.fit
    set input_sextractor = ${DASCH_MATCH}/${plate}_tnx.db
    set output_sextractor = ${DASCH_MATCH}/${plate}.db

    if (-e $input_image) then
       echo "getfits -i 0 -o stdout -x 300 300   $input_image 5:34:31.950 +22:00:52.10 J2000 > ~/junk/plates/${plate}.fit"
             getfits -i 0 -o stdout -x 300 300   $input_image 5:34:31.950 +22:00:52.10 J2000 > ~/junk/plates/${plate}.fit
    endif
end
date
