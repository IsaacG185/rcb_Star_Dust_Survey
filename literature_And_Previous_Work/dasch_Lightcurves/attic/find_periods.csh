# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
# This script controls the complete process of generating a power density spectrum 
# (lomb-scargle periodogram) detecting significant periods and making folded lightcurves
# The input file is the lightcurve (Time, Magnitude, Error)
# Outputs are:  power-spectrum (frequency, Power)
#		peaks list (Period,Amplitude,period_errors(1-5), significance)		
#		folded lightcurves, zero or more per object, folded at the detected periods.
#
# The fifth argument, if 1, enables the scargle version
# The sixth argument, if 1, enables the pdf version
#
#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#  Mar 24, 2008 Edward J. Los - Move magnitude from third column to fifth column

set scripts = $DASCH_SCRIPTS
#setenv LD_LIBRARY_PATH /usr/local/pgplot


set lc = $1
set p1 = $2
set p2 = $3
set intF = $4
set scargle = $5
set pdm = $6
set threshold = 90
set bins = 10

if ( $#argv != 6 ) then
	echo "usage: lc p1 p2 intF scargle pdm"
	exit
endif

if ($scargle == 1) then
  set powspec_scargle = ${lc}.pds
  set peaks_scargle = ${lc}.peaks

  touch scarglepars.tmp
  touch peakpars.tmp
endif

if ($pdm == 1) then

  set powspec_pdm = ${lc}.pdm
  set peaks_pdm = ${lc}.best

  touch pdmpars.tmp

endif

# sort lightcurve in TIME order
sort -n -k 1,1 $lc >! sorted.tmp
mv sorted.tmp $lc

if ($scargle == 1) then
  echo $lc >> scarglepars.tmp
  echo $p1 $p2 $intF >> scarglepars.tmp

  if (-e scarglefast_logfile) then
    rm scarglefast_logfile
  endif

  echo "executing scarglefast_SFS2"
  ${scripts}/scarglefast_SFS2 < scarglepars.tmp


  set variance = `grep variance VARMNT | gawk '{print $2}'`
  set M = `grep "independent_freq" VARMNT | gawk '{print $2}'`
  set N = `grep datapoints VARMNT | gawk '{print $2}'`
  set exptime = `grep Timespan VARMNT  | gawk '{print $2}'`
  
  mv VARMNT ${lc}.VARMNT
  mv scarglefast_logfile ${lc}.logfile

  if (-e scarglefast_pds.txt)  then 
    mv scarglefast_pds.txt $powspec_scargle
  
    echo $powspec_scargle >> peakpars.tmp
    echo $peaks_scargle >> peakpars.tmp
    echo $M $variance $N $exptime $threshold $p1 $p2 >> peakpars.tmp
  
    echo "executing find_periods"
    ${scripts}/find_periods < peakpars.tmp
  
  # If any significant peaks were found, fold the lightcurve at the best period.
    if ( `cat $peaks_scargle | wc -l` > 1 ) then
      set period = `grep -v "Period" $peaks_scargle | sort -nr -k 2,2 | head -1 | gawk '{print $1}'`
      set folded = ${lc}.p=${period}.fold
      cat $lc | gawk '{print $2, $5, 0.0}' | sort -n -k 1,1 >! lc_sorted.tmp
      echo lc_sorted.tmp >! foldpars.tmp
      echo $period 0.0 >> foldpars.tmp
      echo "Executing fold"
      ${scripts}/fold < foldpars.tmp >! $folded	
      cat ${scripts}/plot_fold_template.wip | sed "s/FILE/${folded}/g" | sed "s/PERIOD/${period}/g" >! ${folded}.wip
      echo "creating ${folded}.wip"
      rm lc_sorted.tmp
      rm foldpars.tmp
     #wip ${folded}.wip
     endif
  
  
  
    echo "creating ${powspec_scargle}.wip"
    cat ${scripts}/plot_pds.wip | sed "s/DATA/${powspec_scargle}/g" >! ${powspec_scargle}.wip
    set L99 = `grep "99%_level" ${lc}.VARMNT | gawk '{print $2}'`
    set L95 = `grep "95%_level" ${lc}.VARMNT | gawk '{print $2}'`
    set L90 = `grep "90%_level" ${lc}.VARMNT | gawk '{print $2}'`
    echo 0.000001 $L99 $L95 $L90 >! ${powspec_scargle}.limits
    echo 10000000 $L99 $L95 $L90 >> ${powspec_scargle}.limits
    #wip ${powspec_scargle}.wip
  endif
     rm scarglepars.tmp 
     rm peakpars.tmp


endif

if ($pdm == 1) then
  if ( -e pdm_result.tmp) then 
    rm pdm_result.tmp
  endif
  echo $lc >> pdmpars.tmp
  echo $bins >> pdmpars.tmp
  echo $p1 $p2 $intF >> pdmpars.tmp
  echo 'pdm_result.tmp' >> pdmpars.tmp
  ${scripts}/pdm_freq < pdmpars.tmp > ${peaks_pdm}
  if ( -e pdm_result.tmp) then

    mv pdm_result.tmp ${lc}.pdm
  
    # If any significant peaks were found, fold the lightcurve at the best period.
    if ( `cat $peaks_pdm | wc -l` > 1 ) then
      set period = `grep "Period" $peaks_pdm | gawk '{print $3}'`
      set folded = ${lc}.p=${period}.pdm.fold
      cat $lc | gawk '{print $2, $5, 0.0}' | sort -n -k 1,1 >! lc_sorted.tmp
      echo lc_sorted.tmp >! foldpars.tmp
      echo $period 0.0 >> foldpars.tmp
      ${scripts}/fold < foldpars.tmp >! $folded	
      cat ${scripts}/plot_fold_template.wip | sed "s/FILE/${folded}/g" | sed "s/PERIOD/${period}/g" >! ${folded}.wip
      #wip ${folded}.wip
         rm foldpars.tmp
         rm lc_sorted.tmp
    endif
  endif
  
  cat ${scripts}/plot_pdm.wip | sed "s/DATA/${powspec_pdm}/g" >! ${powspec_pdm}.wip
  
  #wip ${powspec_pdm}.wip
  
     rm pdmpars.tmp
endif
