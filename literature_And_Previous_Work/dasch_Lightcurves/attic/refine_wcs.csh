# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
# Script to refine the astrometric WCS for images
# Begin with a fits image that has an approximate WCS
# (1) Run Sextractor
# (2) Match catalog against GSC2.2
# (3) compute astrometric errors
# (4) Run IRAF CCMAP to fit a new WCS (TNX, order 6,6,full)
# (5) Repeat steps 2 & 3 to assess improvement
#
#
#  NOTE: AS OF Jun 6, 2009, THIS FILE USES TRICKS INTO FOOLING IRAF TO HANDLE A 2GB FILE
#        1.  Only the header is presented to IRAF
#        2.  BITPIX is changed to 8 before running IRAF and back to 16 after running IRAF
#        3.  The FitsHeader program recombines the new header with the old body
#
#
#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#  Nov  7, 2007 Edward J. Los - Operate on the mosaics in their working directory
#  Mar 25, 2008 Edward J. Los - Use ra_2 and dec_2 instead of ra_1 and dec_1
#                               Support match files in $DASCH_MATCH
#  Aug 26, 2008 Edward J. Los - Add filter_ccmap to pick good images for refinement
#  Oct 10, 2008 Edward J. Los - Improve error handling for zero size *.coo files
#  Oct 14, 2008 Edward J. Los - correct zero size error
#                               back out filter_ccmap
#  Jun  6, 2009 Edward J. Los - Avoiding copying or writing the entire file because of an issue with
#                               2GB+ mosaics
#  Jun 15, 2009 Edward J. Los - Change .hdr to .hdr8 to avoid confusion with BITPIX = 16 saved headers
#  Sep 29, 2009 Edward J. Los - add multiple exposure support
#  Oct 15, 2009 Edward J. Los - replace with scamp

echo "ERROR: refine_wcs.csh is deprecated"
exit

set scripts = $DASCH_SCRIPTS

if ($#argv != 2) then
    echo "usage: refine_wcs.csh  solutionNumber list"
    exit
endif

set solutionNumber = $1
set imagelist = $2
set local = `pwd`
set run_filter_ccmap = 1
set plot_filter_ccmap = 0

if ($solutionNumber == 0) then
    set solutionString = ""
else
    set solutionString = "_s$solutionNumber"
endif

if (-e ${local}/${imagelist:t:r}${solutionString}_refine.tmp) then
   rm  ${local}/${imagelist:t:r}${solutionString}_refine.tmp
endif

echo "Processing $imagelist"
date
foreach plate ( `cat $imagelist` )
  set inputfile = $plate.fit
  set directory = `$DASCH_SCRIPTS/getdirectory $inputfile`
  

  if ($solutionNumber == 0) then
    cd $directory
    pwd

#
#  Clean up old fits files if present
#
    if (-e ${plate:r}${solutionString}.fit) then
#      Here our base file exists
       if (-e ${plate:r}${solutionString}_tnx.out) then
          rm ${plate:r}${solutionString}_tnx.out
       endif
       if (-e ${plate:r}${solutionString}_tnx.fit) then
          rm ${plate:r}${solutionString}_tnx.fit
       endif
    else
#      Here our base file does not exist.  
       if (-e ${plate:r}${solutionString}_tnx.fit) then
          mv ${plate:r}${solutionString}_tnx.fit ${plate:r}${solutionString}.fit
          if (-e ${plate:r}${solutionString}_tnx.out) then
             rm ${plate:r}${solutionString}_tnx.out
          endif
       else 
          if (-e ${plate:r}${solutionString}_tnx.out) then
             mv ${plate:r}${solutionString}_tnx.out  ${plate:r}${solutionString}.fit
          endif
       endif
    endif

  else
    cd $DASCH_HEADERS
    pwd
    if (-e ${plate:r}${solutionString}_tnx.out) then
      rm ${plate:r}${solutionString}_tnx.out
    endif
    if (-e ${plate:r}${solutionString}_tnx.hdr) then
      rm ${plate:r}${solutionString}_tnx.hdr
    endif
  endif

  if (-e ${plate:r}${solutionString}_tnx.fits) rm ${plate:r}${solutionString}_tnx.fits 


# Step (2) Match catalog against GSC2.2 and extract co-ordinate list.
# (for now, just use existing match-lists from previous run)
 
  set mosaicsize = `$DASCH_SCRIPTS/getlocation -a $plate -e $solutionNumber`
  set width = `echo "$mosaicsize" | gawk '{print $1}'`
  set height  = `echo "$mosaicsize" | gawk '{print $2}'`
  if ($height == "") then
#    echo "ERROR getlocation failed with args -a $plate -e $solutionNumber"
    continue
  endif

  if (-e ${plate:r}${solutionString}.coo) rm ${plate:r}${solutionString}.coo
  if (!(-e ${DASCH_MATCH}/match_${plate:r}${solutionString}_u.db)) then
    echo "ERROR: ${DASCH_MATCH}/match_${plate:r}${solutionString}_u.db does not exist"
    continue
  endif


  # Match everything
    if ($run_filter_ccmap == 0) then
      column -i ${DASCH_MATCH}/match_${plate:r}${solutionString}_u.db -b X_IMAGE Y_IMAGE ra_2 dec_2 >! ${plate:r}${solutionString}.coo  
    else
  # Selective Match Pipeline V3.4.11
      if (-e ${DASCH_MATCH}/${plate:r}${solutionString}_ccmap.db) rm ${DASCH_MATCH}/${plate:r}${solutionString}_ccmap.db 
      echo "$DASCH_SCRIPTS/filter_ccmap -v -i ${DASCH_MATCH}/match_${plate:r}${solutionString}_u.db -o ${plate:r}${solutionString}.coo -l 2 -u 9 -t -w $width -h $height -s ${DASCH_MATCH}/${plate:r}${solutionString}_ccmap.db -g ${DASCH_MATCH}/${plate:r}${solutionString}_ccmapgmt.txt"
      $DASCH_SCRIPTS/filter_ccmap -i ${DASCH_MATCH}/match_${plate:r}${solutionString}_u.db -o ${plate:r}${solutionString}.coo -l 2 -u 9 -t -w $width -h $height -s ${DASCH_MATCH}/${plate:r}${solutionString}_ccmap.db 

  # Create plots
      if ($plot_filter_ccmap == 1) then
        source ${DASCH_SCRIPTS}/initgmt
        minmax ${DASCH_MATCH}/${plate:r}${solutionString}_ccmapgmt.txt
#         psxy   ${DASCH_MATCH}/${plate:r}${solutionString}_ccmapgmt.txt -R0/23000/0/23000 -Y4i -JX6i/5.5i -P -B20000:"X_IMAGE":/20000:"Y_IMAGE":WSne -Sy0.05 -K >!  ${DASCH_MATCH}/${plate:r}${solutionString}.ps

        psxy   ${DASCH_MATCH}/${plate:r}${solutionString}_ccmapgmt.txt -R0/6000/0/100 -Y4i -JX6i/5.5i -P -B1000:"FLUX_MAX (ADU)":/10:"FWHM_IMAGE (sq pixels)":WSne -Sy0.05 -K >!  ${DASCH_MATCH}/${plate:r}${solutionString}.ps
        pwd
        ps2pdf ${DASCH_MATCH}/${plate:r}${solutionString}.ps 
        mv ${plate:r}${solutionString}.pdf ${DASCH_MATCH}/${plate:r}${solutionString}.pdf
      endif

      if (-z ${plate:r}${solutionString}.coo) then 
         echo "WARNING: Reverting to old refine_wcs algorithm ${plate:r}${solutionString}.coo "

         echo "column -i ${DASCH_MATCH}/match_${plate:r}${solutionString}_u.db -b X_IMAGE Y_IMAGE ra_2 dec_2 >! ${plate:r}${solutionString}.coo"  
         column -i ${DASCH_MATCH}/match_${plate:r}${solutionString}_u.db -b X_IMAGE Y_IMAGE ra_2 dec_2 >! ${plate:r}${solutionString}.coo  
      endif


    endif
    #echo "ERROR: early exit 1"
    #exit

  if (-z ${plate:r}${solutionString}.coo) then 
    echo "ERROR: ${plate:r}${solutionString}.coo has zero size"
    continue
  endif
 
  if ($solutionNumber == 0)  then
    set sourcefile = ${directory}/${plate}${solutionString}.fit
    set coofile = ${directory}/${plate}${solutionString}.coo
    set outfile = ${directory}/${plate}${solutionString}_tnx
  else
    set coofile = ${DASCH_HEADERS}/${plate}${solutionString}.coo
    set sourcefile = ${DASCH_HEADERS}/${plate}${solutionString}.hdr
    set outfile = ${DASCH_HEADERS}/${plate}${solutionString}_tnx
  endif  
  if (!(-e $sourcefile)) then
    echo "ERROR: $sourcefile does not exist"
    continue
  endif


  echo "$sourcefile $coofile  $outfile $outfile" >>!  ${local}/${imagelist:t:r}${solutionString}_refine.tmp

  #cp ${plate:r}${solutionString}.fit ${plate:r}${solutionString}_tnx.fits  

  echo "FitsHeader -i ${sourcefile} -p 8 -o  ${plate:r}${solutionString}_tnx.fits"
  FitsHeader -i ${sourcefile} -p 8 -o  ${plate:r}${solutionString}_tnx.fits 

# Step (3) compute astrometric errors and plot.
# (already done from previous run)    
  
end

###echo "ERROR: early exit 2"
###exit


date
# Step (4) Run IRAF CCMAP to fit a new WCS (TNX, order 6,6,full)
# Invoke IRAF 

if !(-e ${local}/${imagelist:t:r}${solutionString}_refine.tmp) then
    echo "ERROR: List file does not exist: ${local}/${imagelist:t:r}${solutionString}_refine.tmp"
    exit
endif

limit stacksize unlimited
cd ~/iraf
echo "Invoking iraf"
cl << ENDCL > ${local}/${imagelist:t:r}${solutionString}.log
  cd $local

  set uparm = $scripts/
  string image, coords, out, wcsdir, s4
  struct *list1
  
  list1 = "${imagelist:t:r}${solutionString}_refine.tmp"
  time
  while ( fscan(list1,s1,s2,s3,s4) != EOF ) {
    pwd
    image=(s1)
    coords=(s2)
    out=(s3)
    wcsdir=(s4)

    time
    print ("in: "//image//"  out: "//out//" coords: "//coords//" wcsdir: "//wcsdir)

    
#    imdelete (out//".fits")
#    imcopy (image,out//".fits")
#    time
 
#    ccmap(xxorder=6,xyorder=6,yxorder=6,yyorder=3,xxterms="full",yxterms="full",input=coords,database=wcsdir//".wcs",images=out//".fits",inter-,update+)
#    ccmap(input=coords,database=wcsdir//".wcs",images=out//".fits",update-)
     ccmap(input=coords,database=wcsdir//".wcs",images=out//".fits",update+)
    pwd
    time

  }

  logout

ENDCL


date

cd ${scripts}
pwd



foreach plate ( `cat ${imagelist}` )
  set inputfile = $plate.fit
  set directory = `$DASCH_SCRIPTS/getdirectory $inputfile`
  if ($solutionNumber == 0) then
    set bodyString = "-b ${plate:r}${solutionString}.fit"
    set outputFile = ${plate:r}${solutionString}_tnx.fit
    cd $directory
    pwd
  else 
    set bodyString = ""
    set outputFile = ${plate:r}${solutionString}_tnx.hdr
    cd $DASCH_HEADERS
    pwd
  endif
  
    if (!(-e ${plate:r}${solutionString}_tnx.fits)) then
     continue
    endif

    echo "mv ${plate:r}${solutionString}_tnx.fits ${plate:r}${solutionString}_tnx.hdr8"
    mv ${plate:r}${solutionString}_tnx.fits ${plate:r}${solutionString}_tnx.hdr8


    echo "FitsHeader -p 16 ${bodyString} -i ${plate:r}${solutionString}_tnx.hdr8 -o  ${plate:r}${solutionString}_tnx.out -e $solutionNumber"

    FitsHeader -p 16  ${bodyString}  -i ${plate:r}${solutionString}_tnx.hdr8 -o  ${plate:r}${solutionString}_tnx.out -e $solutionNumber 
    rm ${plate:r}${solutionString}_tnx.hdr8

#
#  Clean up fits files
#
    if (-e ${plate:r}${solutionString}_tnx.out) then
#      Here our base file exists
       mv ${plate:r}${solutionString}_tnx.out ${outputFile}
       if (-e ${plate:r}${solutionString}.fit) then
          rm ${plate:r}${solutionString}.fit
       endif
    endif
    if (-e ${plate:r}${solutionString}_tnx.wcs) then
        mv ${plate:r}${solutionString}_tnx.wcs $DASCH_MATCH
    else
      echo "WARNING ${plate:r}${solutionString}_tnx.wcs does not exist"
      pwd
    endif
    if (-e ${plate:r}${solutionString}.coo) then
        mv ${plate:r}${solutionString}.coo $DASCH_MATCH
    endif
    if (-e ${DASCH_SCRIPTS}/${plate:r}${solutionString}.log) then
        mv ${DASCH_SCRIPTS}/${plate:r}${solutionString}.log $DASCH_MATCH
    endif
   
end

if (-e ${local}/${imagelist:t:r}${solutionString}.log) then
    mv ${local}/${imagelist:t:r}${solutionString}.log $DASCH_MATCH
endif

rm  ${local}/${imagelist:t:r}${solutionString}_refine.tmp
