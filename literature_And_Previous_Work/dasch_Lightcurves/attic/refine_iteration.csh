# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  refine_iteration.csh - Prepare for iterative WCS fitting with IRAF ccmap
#
#  NOTE: THIS FILE IS OBSOLETE
#
#  Sep 29, 2008 Edward J. Los - Initial version
#  Oct  7, 2008 Edward J. Los - remove unused plotting commands
#
#  NOTE: THIS FILE IS OBSOLETE
#
set imagelist = $1
set scripts = $DASCH_SCRIPTS
set binning = 16
set useRaDec = 0
date
cd ${DASCH_MATCH}
set listroot = $imagelist:t

pwd
echo "Processing ${DASCH_SCRIPTS}/$listroot"
foreach plate ( `cat ${DASCH_SCRIPTS}/$listroot` ) 
    set directory = `$DASCH_SCRIPTS/getdirectory $plate`  
    set input_image = ${directory}/${plate}_tnx.fit
    set output_image = ${directory}/${plate}.fit
    set input_sextractor = ${DASCH_MATCH}/${plate}_tnx.db
    set output_sextractor = ${DASCH_MATCH}/${plate}.db

    if (-z $input_image) then
      echo "ERROR: refine_iteration $input_image has zero size"
    else
      if (-e $input_image) then
        if (-e $input_sextractor) then
          echo "mv $input_image $output_image"
          mv $input_image $output_image
          echo "${DASCH_SCRIPTS}/update_sextractor -e 0 -v -i $input_sextractor -o $output_sextractor -m $output_image"
          ${DASCH_SCRIPTS}/update_sextractor -e 0 -v -i $input_sextractor -o $output_sextractor -m $output_image
        else
          echo "ERROR: refine_iteration can not find $input_sextractor"
        endif
      else
        echo "ERROR: refine_iteration can not find $input_image"
      endif
    endif
end
date
