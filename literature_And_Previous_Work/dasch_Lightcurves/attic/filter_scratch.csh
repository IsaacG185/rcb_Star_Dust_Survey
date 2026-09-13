# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  filter_scratch.csh
#
#  From correspondance of Mon 3/02/15 10:47 AM:
#
# Borrowing an idea from Sumin Tang, if there are N objects evenly distributed in an area A of the plate, then the mean 
# separation R between these objects is approximately  PI * (R**2) * N ~ A, or  R ~ sqrt(A/(PI*N)).  In the attached figures, 
# I set R = K*sqrt(A/(PI*N)) where experimentation shows that K = 0.3 appears to be the minimum necessary to form reasonably 
# long chains from the scratches.  To find N without including lots of grain noise while recognizing the presence of multiple 
# exposures, I selected a cutoff MAG_ISO such that N = 4 * (number of images matched to the GSC2.3.2 catalog).    I found that 
# for deeper scratches, SExtractor already recognizes components of the scratch through the ELLIPTICITY, THETA_IMAGE, aLength, 
# and bLength parameters.  I therefore plotted rectangles using these SExtractor parameters, but expanding aLength and bLength 
# by R.  Results are given for the six examples listed below where the defect filter did not work.  Two of the images (i49523 
# and ac17627) were too small for this algorithm to be effective, and might be addressed by enhancements in Sumin's defect algorithm.  
# The other four images do show good chains that do not require straight line, vertical defects.  Needless to say, there are a lot of 
# false positives in the images.  Additional experimental parameters will be needed to reduce the number of false positives.  These 
# additional parameters could include the length of the chain and the percentage of area of the bounding rectangle occupied by the 
# chain.  The additional parameters, however, could also raise the false negative rate for groups of parallel scratches.
#
#  Inputs:
#    ${DASCH_MATCH}/match_${name}_tnx_u.db - All objects matched with the GSC catalog
#    ${DASCH_MATCH}/${name}_tnx.db - All objects found by sextractor.
#
#
#  Mar  9, 2015 Edward J.Los - Initial version, adopted from filter_defect
#
set local = `pwd`
set totalcounter = 0
set doregion = 1
set plotstr = ""
set solutionNumber = $1
set imagelist = $2
set magisodiff = ""
set magisodiff = "-b -0.625" 
if ($#argv == 2) then
  set qualifier = ""
  set qualcmd = ""
else 
  if ($#argv == 3) then
      set qualifier = `echo "_$3"`
      set qualcmd = "-q $qualifier" 
  else
     echo "usage: filter_scratch.csh solutionNumber list [qualifier]"
     echo "       solutionNumber is the multiple exposure iteration"
     echo "       list      is the list of mosaics to process"
     echo "       qualifier is a string appended to output filenames"
      exit
  endif
endif

if ($solutionNumber != 0) then
#   We are done.  multiple solutions work only with headers, not the full fits file.
    exit
endif

    set solutionString = ""



set catalog = $DASCH_CATALOG
set scripts = $DASCH_SCRIPTS
set completedFlag = `echo "$DASCH_COMPLETED" | wc -c`
set dateval = `date`
echo "$dateval Processing ${imagelist}${qualifier}"
foreach name ( `cat $imagelist` ) 

  set mosaicsize = `$DASCH_SCRIPTS/getlocation -a $name -e $solutionNumber`
  set width = `echo "$mosaicsize" | gawk '{print $1}'`
  set height  = `echo "$mosaicsize" | gawk '{print $2}'`
  set platescale =    `echo "$mosaicsize" | gawk '{print $3}'`
  if ($height == "") then
    echo "ERROR filter_scratch getlocation failed with args -e solutionNumber -a $name"
    continue
  endif
  if ($platescale == "") then
    echo "ERROR filter_scratch getlocation failed with args -a $image"
    continue
  endif
 

  set matchfilename = "${DASCH_MATCH}/match_${name}${solutionString}_tnx${qualifier}_u.db"
  set sextractorfilename = "${DASCH_MATCH}/${name}${solutionString}_tnx.db"
  set outputfilename = "${DASCH_MATCH}/${name}${solutionString}${qualifier}_scratch.db"
  set summaryfilename = "${DASCH_MATCH}/${name}${solutionString}${qualifier}_summary_scratch.db"

 
  if (-e $outputfilename) then
    rm $outputfilename
  endif
  if (-e $summaryfilename) then
    rm $summaryfilename
  endif

  if ($doregion == 1) then
    set directory = `$DASCH_SCRIPTS/getdirectory $name`
    set plotstr = " -g $directory/${name}${solutionString}${qualifier}_scratch.reg"
  endif



  @ totalcounter = $totalcounter + 1
  set dateval = `date`

  echo "Processing ${name}${solutionString} at $dateval"

  cd $DASCH_MATCH
  pwd

  if (-e $matchfilename) then 


      if (-e $sextractorfilename) then 

          echo "$DASCH_SCRIPTS/filter_scratch -c 0.30 $magisodiff -v $plotstr -f $qualcmd -p ${name}${solutionString} -t $platescale  -w $width -h $height -m $matchfilename -s $sextractorfilename -o $outputfilename"
          $DASCH_SCRIPTS/filter_scratch       -c 0.30 $magisodiff -v $plotstr -f $qualcmd -p ${name}${solutionString} -t $platescale  -w $width -h $height -m $matchfilename -s $sextractorfilename -o $outputfilename

 
      else
  
         echo "ERROR: No filter_scratch sextractor file found: $sextractorfilename"
  
      endif


  else

     echo "ERROR: No filter_scratch match file found: $matchfilename"

  endif


end
set dateval = `date`
echo "$dateval Completed ${imagelist}${qualifier}"

