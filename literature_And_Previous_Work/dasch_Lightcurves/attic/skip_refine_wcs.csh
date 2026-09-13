# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
# This is a placeholder for refine_wcs.csh.  It simply
# renames the TAN file to TNX
#

set scripts = $DASCH_SCRIPTS

set imagelist = $1
set local = `pwd`
set completedFlag = `echo "$DASCH_COMPLETED" | wc -c`
set run_filter_ccmap = 0
set plot_filter_ccmap = 0

if (-e ${local}/${imagelist:t:r}_refine.tmp) then
   rm  ${local}/${imagelist:t:r}_refine.tmp
endif

echo "Processing $imagelist"
date
foreach plate ( `cat $imagelist` )
  set inputfile = $plate.fit
  if ($completedFlag > 2) then
      set directory = $DASCH_COMPLETED
      set subdir = ${directory}/${plate}
  else
      set directory = `$DASCH_SCRIPTS/getdirectory $inputfile`
  endif

  cd $directory
  pwd

#
#  Clean up old fits files if present
#
    if (-e ${plate:r}.fit) then
#      Here our base file exists
       if (-e ${plate:r}_tnx.fits) then
          rm ${plate:r}_tnx.fits
       endif
       if (-e ${plate:r}_tnx.fit) then
          rm ${plate:r}_tnx.fit
       endif
       mv ${plate:r}.fit ${plate:r}_tnx.fit
    else
#      Here our base file does not exist.  
       if (-e ${plate:r}_tnx.fit) then
          echo "skip_refine_wcs: ${plate:r}_tnx.fit already exists"
          if (-e ${plate:r}_tnx.fits) then
             rm ${plate:r}_tnx.fits
          endif
       else 
          if (-e ${plate:r}_tnx.fits) then
             echo "skip_refine_wcs: ${plate:r}_tnx.fit already exists"
             mv ${plate:r}_tnx.fits  ${plate:r}_tnx.fit
          endif
       endif
    endif
end


date

