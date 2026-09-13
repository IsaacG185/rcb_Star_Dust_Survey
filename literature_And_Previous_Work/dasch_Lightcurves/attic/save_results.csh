# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  save_results.csh.  
#
#
#
#    if DASCH_COMPLETED is defined as a valid directory, copy the mosaic
#    to that directory.
#    
#  Nov  10, 2007 Edward J. Los - First version to replace copyback.perl
#  Jan   8, 2008 Edward J. Los - Disable this script until the pipeline is stable
#
echo "DO NOT USE save_results.csh TO AVOID CONFUSION OF PIPELINE VERSIONS"
exit
set imagelist = $1
set scripts = $DASCH_SCRIPTS
set data = $DASCH_COMPLETED
set ingestdir = $DASCH_INGEST
set counter = 1
set completedFlag = `echo "$DASCH_COMPLETED" | wc -c`
echo "Processing $imagelist Copy location $DASCH_COMPLETED"
foreach plate ( `cat $imagelist` ) 
  set datestr = `date`
  echo "$datestr $counter $plate"
  @ counter = $counter + 1

  set directory = `$DASCH_SCRIPTS/getdirectory $plate`
  set origfilename = $directory/${plate}.fit 
  set tnxfilename = $directory/${plate}_tnx.fit

#
# Begin by copying the mosaic from the completed directory
# assuming that we are using the completed directory
#  

  if ($completedFlag > 2) then
    set origcomplname = $DASCH_COMPLETED/${plate}.fit 
    set tnxcomplname = $DASCH_COMPLETED/${plate}_tnx.fit
    set matchname = $DASCH_COMPLETED/${plate}_tnx/match_${plate}_tnx_u.db
    if (-e $tnxcomplname) then
        echo "$tnxcomplname found"
        mv $tnxcomplname $tnxfilename
        if (-e $origcomplname) then
            rm $origcomplname
        endif
    else
        echo "$tnxcomplname does not exist"
    endif
    if (-e $matchname) then
        echo "cp $matchname $directory"
        cp $matchname $directory
    else 
       echo "$matchname does not exist"
    endif
  endif
#
# Now remove the old mosaic if we have successfully copied
#  the new one
#
  if (-e $tnxfilename) then
     echo "$tnxfilename found"
     if (-e $origfilename) then
         rm $origfilename
     endif
  else
     echo "$tnxfilename does not exist"
  endif
#
#   Now copy over the octave results
#
    set gridfile = `ls ${DASCH_BINOUTPUT}/${plate}_a?.grid`
    if ($#gridfile > 0) then
        foreach gridfile (${DASCH_BINOUTPUT}/${plate}_a?.grid)
           set root = ${gridfile:r}

           cp $gridfile $directory
           cp ${root}.grid $directory
           cp ${root}_a.eps $directory
           cp ${root}_b.eps $directory
        end

    else

        echo "${DASCH_BINOUTPUT}/${plate}_an.grid not found"

    endif
#
#  Copy the ingested data
#
#
    if (-e ${DASCH_INGEST}/${plate}.out.spatial_bins.db) then
        cp ${DASCH_INGEST}/${plate}.out.spatial_bins.db $directory
        cp ${DASCH_INGEST}/${plate}.out.local.db $directory
    else
        echo "${DASCH_INGEST}/${plate}.out.spatial_bins.db does not exist"

    endif
end
