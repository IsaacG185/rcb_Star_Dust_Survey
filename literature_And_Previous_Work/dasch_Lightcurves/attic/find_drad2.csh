# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  find_drad2.csh   This script extracts X_IMAGE Y_IMAGE dra and ddec from
#  a list of plates presumably all in the same series and with the same 
#  scan pattern and creates a table for plotting these deviations as a 
#  function of position on the plate.
#
#
#  Feb  8, 2008 Edward J. Los - Initial version
#  May 31, 2008 Edward J. Los - Correct for new mosaic location
#  Jan 28, 2009 Edward J. Los - Split FLAGS into AFLAGS and BFLAGS
# 


set imagelist = $1
set catalog = $DASCH_CATALOG
set scripts = $DASCH_SCRIPTS
if (-e find_drad.db) then
  rm find_drad.db
endif
   echo "X_IMAGE Y_IMAGE dra ddec" >! find_drad2.tmp
   echo "------- ------- --- ----" >> find_drad2.tmp


echo "Processing $imagelist"
foreach plate ( `cat $imagelist` ) 
#  set name = `echo $plate | sed 's/_tnx.fits//g'`
  set name = $plate
  echo "Processing $name"
  set filename = ${DASCH_MATCH}/match_${plate}_tnx_u.db
  set platelocation = `$DASCH_SCRIPTS/getlocation -s $plate`
  set series = `echo "$platelocation" | gawk '{print $1}'`
  set scale  = `echo "$platelocation" | gawk '{print $2}'`
  set radius = `echo "$platelocation" | gawk '{print $3}'`
  if (-e $filename) then
    row  -i $filename 'AFLAGS < 1' | column -b X_IMAGE Y_IMAGE dra ddec >> find_drad2.tmp
  endif

end
cat find_drad2.tmp | gawk '{OFS="\t"}{print $1,$2,$3,$4}' >! find_drad2.db
rm find_drad2.tmp
