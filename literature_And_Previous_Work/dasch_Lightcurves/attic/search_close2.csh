# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  search_close2.csh - Look for close multiple exposures
#
#  Dec 20, 2011 Edward J. Los - Adapted from filter_wedge.csh
#  Jun  1, 2012 Edward J. Los - Add background search of calibration squares
#  Jul  3, 2012 Edward J. Los - Integrate high background object support with rest of pipeline
#  Jul 23, 2012 Edward J. Los - Skip this for solutionNumber != 0;
#  Oct 24, 2017 Edward J. Los - set the analysis threshold to zero
#
if ($#argv != 2) then
    echo "usage: search_close2.csh <solutionNumber> <list>"
    exit
endif
#if ($#argv != 2) then
#    echo "usage: search_close2.csh <list> <analysis_threshold>"
#    exit
#endif
set solutionNumber = $1
set imagelist = $2
set scripts = $DASCH_SCRIPTS
set binning = 16
set doplots = 1
#set dobackground = 1
set dobackground = 0 
set plotstr = ""
set calibrationSquare = " "
set analysis_threshold = 2.2
#set analysis_threshold = 0.

if ($doplots == 1) then
    set plotstr = " -q -p "
endif
cd ${DASCH_MATCH}
set listroot = $imagelist:t

if ($solutionNumber != 0) then
#   We are done.  this code is independent of the solution number
    exit
endif


pwd
echo "Processing ${DASCH_SCRIPTS}/$listroot"
foreach plate ( `cat ${DASCH_SCRIPTS}/$listroot` ) 

  set backgroundname = ${DASCH_MATCH}/${plate}_background.db
  if (-e $backgroundname) rm $backgroundname
  set backgroundtmp = ${DASCH_MATCH}/${plate}_background.tmp
  if (-e $backgroundtmp) rm $backgroundtmp

  if ($dobackground == 1) then
     set mosaicDirectory = `getdirectory $plate`
     set calibrationSquare = " -b 1 -d $backgroundtmp ${mosaicDirectory}/background.reg"
  endif

  set series = `$DASCH_SCRIPTS/getlocation -s $plate | gawk '{print $1}'`

  set outputname = ${DASCH_MATCH}/${plate}_close.db
  if (-e $outputname) rm $outputname
  set gmt0name = ${DASCH_MATCH}/${plate}_close_gmt0.txt
  set gmt1name = ${DASCH_MATCH}/${plate}_close_gmt1.txt
  set gmt2name = ${DASCH_MATCH}/${plate}_close_gmt2.txt
  set gmt3name = ${DASCH_MATCH}/${plate}_close_gmt3.txt
  set gmt4name = ${DASCH_MATCH}/${plate}_close_gmt4.txt
  set emulsionname = ${DASCH_MATCH}/${plate}_close_plot.db

  if (-e $gmt0name) rm $gmt0name
  if (-e $gmt1name) rm $gmt1name
  if (-e $gmt2name) rm $gmt2name
  if (-e $gmt3name) rm $gmt3name
  if (-e $gmt4name) rm $gmt4name
  if (-e $emulsionname) rm $emulsionname


  if (-e ${DASCH_MATCH}/${plate}_tnx.db) then
     set inputname = ${DASCH_MATCH}/${plate}_tnx.db
  else 
    if (-e ${DASCH_MATCH}/${plate}.db) then
       set inputname = ${DASCH_MATCH}/${plate}.db
    else
       echo "ERROR: search_close2 found no sextractor file for $plate"
       continue
    endif
  endif

  

  set mosaicsize = `$DASCH_SCRIPTS/getlocation -a $plate `
  set width = `echo "$mosaicsize" | gawk '{print $1}'`
  set height  = `echo "$mosaicsize" | gawk '{print $2}'`
  set platescale =    `echo "$mosaicsize" | gawk '{print $3}'`
  if ($platescale == "") then
#       echo "ERROR getlocation failed with args -a $image"
    continue
  endif


  if (-e $inputname) then
    echo "$DASCH_SCRIPTS/search_close2  -v $plotstr -t $analysis_threshold  $calibrationSquare  -s $platescale -w $width -h $height -r ${plate} -i $inputname -o $outputname"
    #echo "ERROR: Early exit"
    #exit
    $DASCH_SCRIPTS/search_close2  $plotstr  -t $analysis_threshold  $calibrationSquare -w $width -s $platescale -h $height -r ${plate} -i $inputname -o $outputname

    if (-e $backgroundtmp) then
        #echo "sorttable -n NUMBER < $backgroundtmp > $backgroundname"
        sorttable -n NUMBER < $backgroundtmp > $backgroundname
        rm $backgroundtmp

    endif


  else 
      echo "ERROR: Can not find input file $inputname"
  endif
  
end
