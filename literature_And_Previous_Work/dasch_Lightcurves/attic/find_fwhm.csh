# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  find_fwhm.csh.  
#
#   Search for match_<plate name> files
#
#    
#  Dec  5, 2007 Edward J. Los - First version
#  Apr 16, 2011 Edward J. Los - Change /media to /dasch
#
#  TODO: Will not work for 14x17 plates
#  TODO: Reject galaxies
#
setenv DASCH_COMPLETED "/dasch/mosaic/Pipeline/completed"

set imagelist = $1
set scripts = $DASCH_SCRIPTS
set data = $DASCH_COMPLETED
set counter = 1
set completedFlag = `echo "$DASCH_COMPLETED" | wc -c`
echo "Processing $imagelist Copy location $DASCH_COMPLETED"
foreach plate ( `cat $imagelist` ) 
  set datestr = `date`
  echo "$datestr $counter $plate"
  @ counter = $counter + 1
  set foundflag = 0;  
  set directory1 = `$DASCH_SCRIPTS/getdirectory $plate`
  set directorylist = $directory1
  if ($completedFlag > 2) then
    set directory2 = $DASCH_COMPLETED/${plate}_tnx
    set directory3 = $DASCH_COMPLETED/${plate}
    set directorylist = "$directory1 $directory2 $directory3"
  endif
  set filelist = "match_${plate}_tnx.db match_${plate}.db match_${plate}_tnx_u.db match_${plate}_u.db"
  foreach filetest ($filelist)
     foreach directorytest ($directorylist)
        if ($foundflag == 0) then
            set testname = "${directorytest}/${filetest}"
            if (-e $testname) then
#               echo "Found  $testname"
                if ($foundflag == 0) then
                    set foundflag = 1
                    set fileroot = $filetest
                    set filename =  $testname
                    set directory = ${directorytest}
                endif
            else
#               echo "Not Found  $testname"
            endif
        endif
     end
  end
 if ($foundflag > 0) then
    echo "Using file: $filename"
 else
    echo "No match file for ${plate}"
    continue
 endif
 set horizontalFlag = `echo "$filename" | egrep "r90|r270"`
 set horizontalFlag = $#horizontalFlag
# 1: extract fields of interest
# 2: exclude 10% margins
# 3: sort reverse by Stdmag
# 4: take the highest Stdmag

 if ($horizontalFlag == 0) then
  set magline = `column -i $filename X_IMAGE Y_IMAGE Stdmag FWHM_IMAGE | row  '(Y_IMAGE > 2204) && (Y_IMAGE < 19841) && (X_IMAGE > 1741) && (X_IMAGE < 15670)' | sorttable -n -r Stdmag | tail -n 1`
 else
  set magline = `column -i $filename X_IMAGE Y_IMAGE Stdmag FWHM_IMAGE | row  '(X_IMAGE > 2204) && (X_IMAGE < 19841) && (Y_IMAGE > 1741) && (Y_IMAGE < 15670)' | sorttable -n -r Stdmag | tail -n 1`
 endif

  set xposition = `echo $magline | awk '{print $1}'`
  set yposition = `echo $magline | awk '{print $2}'`
  set maxmag    = `echo $magline | awk '{print $3}'`
  set fwhm      = `echo $magline | awk '{print $4}'`
  echo "fwhm $fwhm for Stdmag $maxmag at $xposition $yposition for ${fileroot}"
end
date
