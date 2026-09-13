# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#  Nov  7, 2007 Edward J. Los - Operate on the mosaics in their working directory
#  Nov 23, 2007 Edward J. Los - Support filtering of blended stars
#  Dec  5, 2007 Edward J. Los - Create a new catalog for each plate
#  Sep 29, 2009 Edward J. Los - add multiple exposure support
#  Oct 12, 2009 Edward J. Los - restructure by removing dasch_match
#  Oct 15, 2009 Edward J. Los - replace with scamp

echo "ERROR: run_match.csh is deprecated"
exit

if ($#argv != 2) then
    echo "usage: run_match.csh  solutionNumber list"
    exit
endif

set solutionNumber = $1
set imagelist = $2

set scripts = $DASCH_SCRIPTS

if ($solutionNumber == 0) then
    set solutionString = ""
else
    set solutionString = "_s$solutionNumber"
endif


echo "Processing $imagelist"
foreach plate ( `cat $imagelist` ) 
  echo $plate
#  set name = `echo $plate | sed 's/_tnx.fits//g'`
  set name = $plate
  set directory = `$DASCH_SCRIPTS/getdirectory $plate`
  set filename = ${directory}/${plate}.fit

#     ${scripts}/dasch_match.csh $filename $catalog 10 20 1 1 0 0

  set qualifier = ""
  set qualcmd = ""

  set scripts = $DASCH_SCRIPTS

  set fits_image = $filename
  set image = $filename:t:r

  set maglimit = " -m 17.0"
  set plate_scale = ""
  set pixelRadius = "-r 6"
# NOTE: any change in the search radius should be reflected in MULTIPLE_EXPOSURE_TOLERANCE of 
#       AstrometryWCS.cpp for multiple exposures
 
  set script_output = `$DASCH_SCRIPTS/getlocation  $fits_image $pixelRadius`
  set match_radius = `echo "$script_output" | gawk '{print $5}'`


  set mosaicsize   = `$DASCH_SCRIPTS/getlocation -e $solutionNumber -a $image`
  set width = `echo "$mosaicsize" | gawk '{print $1}'`
  set height  = `echo "$mosaicsize" | gawk '{print $2}'`
  set plate_scale_tycho = `echo "$mosaicsize" | gawk '{print $3}'`
  set leftMargin   = `echo "$mosaicsize" | gawk '{print $5}'`
  set rightMargin  = `echo "$mosaicsize" | gawk '{print $6}'`
  set bottomMargin = `echo "$mosaicsize" | gawk '{print $7}'`
  set topMargin    = `echo "$mosaicsize" | gawk '{print $8}'`

  if ($plate_scale_tycho == "") then
#    echo "ERROR getlocation failed with args -e $solutionNumber -a $image"
    continue
  endif

  set tychoRadius = 20
  set series = `$DASCH_SCRIPTS/getlocation -s $image | gawk '{print $1}'` 
  if (($series == "ac") || ($series == "am") || ($series == "ax") || ($series == "ca")) then 
    set tychoRadius = 60
  endif
  set directory = `$DASCH_SCRIPTS/getdirectory $image`




  set TSTART = `date`
  echo Begin Process $TSTART
  echo fits_image $fits_image
  echo image $image 
  echo match_radius $match_radius $pixelRadius
  echo qualifier $qualifier
  echo qualcmd $qualcmd


#cd ${fits_image:h}
cd $DASCH_MATCH
pwd

  set table = ${image:r}.db 
  set matchtable = match_${table:r}${qualifier}.db


# 3. Match Detected Stars with Reference Stars
  if (-e match_${table}${solutionString}${qualifier}_u.db) rm match_${table}${solutionString}${qualifier}_u.db
  
 



  date
  if ($solutionNumber == 0) then
    set imty2input = ${directory}/${image}.fit
    set imty2output = ${DASCH_MATCH}/${image}.fit.tycho2  
  else
    set imty2input = ${DASCH_HEADERS}/${plate}${solutionString}.hdr
    set imty2output =  ${DASCH_MATCH}/${image}${solutionString}.hdr.tycho2
  endif
  if (!(-e $imty2input)) then
    echo "ERROR: run_match did not find $imty2input"
    continue
  endif
  

  echo "imty2  -w -d -n 5000000 $imty2input"
  imty2  -w -d -n 5000000 $imty2input 
  #mv ${image}.fit.tycho2 ${DASCH_MATCH}
  echo "matchtycho2 -v -w $width -h $height -c $imty2output -i ${DASCH_MATCH}/${image}${solutionString}.db -o ${DASCH_MATCH}/match_${image}${solutionString}${qualifier}_u.db -p $tychoRadius -l $leftMargin -r $rightMargin -b $bottomMargin -t $topMargin -s $plate_scale_tycho -n 5000"
  matchtycho2 -v -w $width -h $height -c $imty2output -i ${DASCH_MATCH}/${image}${solutionString}.db -o ${DASCH_MATCH}/match_${image}${solutionString}${qualifier}_u.db -p $tychoRadius -l $leftMargin -r $rightMargin -b $bottomMargin -t $topMargin -s $plate_scale_tycho -n 5000
  rm  $imty2output


  date



end
