# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
# Test routine for SCAMP mosaic cleanup
#
# 2011-06-07 Edward J. Los Initial version adopted from run_scamp.
# 

if ($#argv != 2) then
    echo "usage: run_scamp.csh  solutionNumber list"
    exit
endif

set solutionNumber = $1
set imagelist = $2
set scripts = $DASCH_SCRIPTS

set plot_drad = 1
set redo_filt = 1
set useTranslateWCS = 1

if ($solutionNumber == 0) then
    set solutionString = ""
else
    set solutionString = "_s$solutionNumber"
endif

set listfile = $imagelist:t
echo "Processing ${DASCH_SCRIPTS}/$listfile"

foreach mosaic ( `cat ${DASCH_SCRIPTS}/$listfile` )

    if ($solutionNumber == 0) then
       set mosaic_dir = `$DASCH_SCRIPTS/getdirectory $mosaic`
       set mosaic_fit = ${mosaic_dir}/${mosaic}.fit
       set mosaic_fit_tnx = ${mosaic_dir}/${mosaic}_tnx.fit
       set imwcsinput = $mosaic_fit
       set headeroutput = $mosaic_fit_tnx
    else 
       set imwcsinput = ${DASCH_HEADERS}/${mosaic}${solutionString}.hdr
       set headeroutput =  ${DASCH_HEADERS}/${mosaic}${solutionString}_tnx.hdr
       set headeroutputtmp =  ${DASCH_HEADERS}/${mosaic}${solutionString}_tnx.hdrtmp
       set mosaic_fit_tnx = ${DASCH_MATCH}/${mosaic}${solutionString}_dummy_file_to_avoid_undefined_variable
       @ prevSolutionNumber = $solutionNumber - 1;
       if ($prevSolutionNumber == 0) then
         set prevSolutionString = ""
       else
         set prevSolutionString = "_s$prevSolutionNumber"
       endif
       set prevTnxFile = ${DASCH_HEADERS}/${mosaic}${prevSolutionString}_tnx.hdr
    endif

    if (-e $mosaic_fit_tnx) then
      set newsize = `ls -s --block-size=1 $mosaic_fit_tnx | awk '{print $1}'`
      echo "$mosaic_fit_tnx exists size $newsize"
        if (-e $mosaic_fit) then
          set oldsize = `ls -s --block-size=1 $mosaic_fit | awk '{print $1}'`
          echo "$mosaic_fit exists $oldsize"
            echo "comparefits $mosaic_fit $mosaic_fit_tnx"
            comparefits $mosaic_fit $mosaic_fit_tnx
            if ($status == 0) then
              echo "comparefits succeeded for $mosaic_fit"
              rm $mosaic_fit
            else
                echo "ERROR $newsize - $oldsize is too small for $mosaic_fit_tnx"

            endif

        endif
    endif
    #echo
    #echo


end
