# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License
#
#  clean_all.csh - Clean up pipeline files after entry into the photometry database
#
#  Mar 26, 2012 Edward J. Los - initial version, adapted from clean_intermediate.ch
#  Jul  3, 2012 Edward J. Los - remove *_background.db
#  Jan 21, 2014 Edward J. Los - Clean filter_defect plots
#  Jun 24, 2018 Edward J. Los - Remover the allobjects.db file
#
#
set imagelist = $1
set scripts = $DASCH_SCRIPTS

if ($#argv == 1) then
  set qualifier = ""
  set qualifier2 = ""
else 
  if ($#argv == 2) then
      set qualifier = `echo "_$2"`
      set qualifier2 = $2
  else
     echo "usage: clean_imtermediate.csh list [qualifier]"
     echo "       list      is the list of mosaics to process"
     echo "       qualifier is a string appended to output filenames"
      exit
  endif
endif


date
set listroot = $imagelist:t
pwd
echo "Processing ${DASCH_SCRIPTS}/$listroot"
set platecount = 0
set totalfiles = 0
foreach plate ( `cat ${DASCH_SCRIPTS}/$listroot` ) 
   
  set solutionNumber = 0
  while (1) 
    set filecount = 0
    if ($solutionNumber == 0) then
        set solutionString = ""
    else
        set solutionString = "_s$solutionNumber"
    endif
    echo "$plate" > clean_all.xxx
    $DASCH_SCRIPTS/clean_intermediate.csh  $solutionNumber clean_all.xxx $qualifier2

    foreach bin (1 2 3 4 5 6 7 8 9)
      set filename = ${DASCH_BINOUTPUT}/${plate}${solutionString}${qualifier}_a${bin}_a.eps
      if (-e $filename) then 
        @ filecount = $filecount + 1
        rm $filename
      endif

      set filename = ${DASCH_BINOUTPUT}/${plate}${solutionString}${qualifier}_a${bin}_b.eps
      if (-e $filename) then 
        @ filecount = $filecount + 1
        rm $filename
      endif

      set filename = ${DASCH_BINOUTPUT}/${plate}${solutionString}${qualifier}_a${bin}.grid
      if (-e $filename) then 
        @ filecount = $filecount + 1
        rm $filename
      endif

      set filename = ${DASCH_BINOUTPUT}/${plate}${solutionString}${qualifier}_a${bin}.para
      if (-e $filename) then 
        @ filecount = $filecount + 1
        rm $filename
      endif

      set filename = ${DASCH_BINOUTPUT}/${plate}${solutionString}${qualifier}_a${bin}.out
      if (-e $filename) then 
        @ filecount = $filecount + 1
        rm $filename
      endif

    #END OF SPATIAL BIN BIN SPECIFIC PROCESSING
    end
 



    set filename = ${DASCH_ASTROMETRY}/${plate}${solutionString}_none.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename = ${DASCH_CATALOGBIN}/${plate}${solutionString}${qualifier}_a${bin}.colorterm.txt
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename = ${DASCH_CATALOGBIN}/${plate}${solutionString}${qualifier}_a${bin}.colorterm.eps
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}_astr_distort_1.png
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}_astr_fgroups_1.png
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}_astr_referr1d_1.png
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}_astr_referr2d_1.png
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}_astrfull_distort_1.png
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}_astrfull_fgroups_1.png
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}_astrfull_referr1d_1.png
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}_astrfull_referr2d_1.png
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}_drad_map.png
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}_drad_offset.png
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}.refcat.ucac4.fits.reg
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}_scamp.sexcat_filtered.fits.reg
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}_scamp.sexcat_filtered.fits.xy.reg
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}.sexcat_filtered.fits.reg
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}.sexcat_filtered.fits.xy.reg
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}_bright1.ps
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}_bright2.ps
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}_brightcalibrate.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}_dmagcor.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}_dmagcor_copy.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.db.dec.i
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.db.ra.b
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.db.REF.i
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.local.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.local.db.dec-i
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.local.db.ra-b
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.local.db.REF-i
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.spatial_bins.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}_dmag.ps
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}_magdep1.ps
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}_magdep2.ps
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_INGEST}/${plate}${solutionString}${qualifier}_magdepcalibrate.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_CATALOGALL}/${plate}${solutionString}${qualifier}_extinction.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/match_${plate}${solutionString}${qualifier}_tnx_b.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/match_${plate}${solutionString}_tnx${qualifier}_b.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/match_${plate}${solutionString}${qualifier}_tnx_lim1.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/match_${plate}${solutionString}${qualifier}_tnx_lim2.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/match_${plate}${solutionString}${qualifier}_tnx_lim3.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/match_${plate}${solutionString}${qualifier}_tnx_lim.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/match_${plate}${solutionString}${qualifier}_tnx_u.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/match_${plate}${solutionString}_tnx${qualifier}_u.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_background.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_close.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_close_plot.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_close_gmt0.txt
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_close_gmt1.txt
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_close_gmt2.txt
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_close_gmt3.txt
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_close_gmt4.txt
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_defect.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_drad.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_summary_defect.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    set filename =  ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_wedge.db
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

		set filelist = `ls ${DASCH_MATCH}/${plate}${solutionString}${qualifier}*_plot.eps`
		echo "filelist is $filelist size $#filelist"
		if ($#filelist > 0) then
		   foreach filename ($filelist)
		     echo $filename
         #if (-e $filename) then 
         #  @ filecount = $filecount + 1
         #  rm $filename
         #endif
       end
    endif


    #END OF SOLUTION NUMBER SPECIFIC PROCESSING   
    @ solutionNumber = $solutionNumber + 1

    #special case - we do not delete the original sextractor file
    set filename =  ${DASCH_MATCH}/${plate}_s${solutionNumber}${qualifier}_tnx.db
    #echo $filename
    if (-e $filename) then 
      @ filecount = $filecount + 1
      rm $filename
    endif

    @ totalfiles = $totalfiles + $filecount
    #set datestr = `date`
    #echo "totalfiles $totalfiles filecount $filecount for ${plate}_s${solutionNumber} at $datestr"
    if ($filecount == 0) then
      break
    endif
    set filecount = 0   

  end

  set filename =  ${DASCH_INGEST}/${plate}_s0${qualifier}_allobjects.db
  if (-e $filename) then 
    @ filecount = $filecount + 1
    rm $filename
  endif

  # added Jun 24, 2018
  set filename =  ${DASCH_INGEST}/${plate}${qualifier}_allobjects.db
  if (-e $filename) then 
    @ filecount = $filecount + 1
    rm $filename
  endif




  #END OF PLATE SPECIFIC PROCESSING
  @ platecount = $platecount + 1
  @ platemod   = $platecount % 100
  if ($platemod == 0) then
    set datestr = `date`
    echo "plates $platecount totalfiles $totalfiles at $datestr"
  endif

end
set datestr = `date`
echo "plates $platecount totalfiles $totalfiles at $datestr"
if (-e clean_all.xxx) rm clean_all.xxx

