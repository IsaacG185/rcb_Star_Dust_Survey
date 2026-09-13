# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#  Nov 14, 2007 Edward J. Los - Eliminate JulianDates.txt and access the MySQL database directly
#  Mar 11, 2007 Edward J. Los - Replace extract_lightcurves_by_id5_local.csh with extract_lightcurves.c
#  Sep 15, 2007 Edward J. Los - Add 3c273 blend
#  Oct 14, 2008 Edward J. Los - Ignore the defect flag when plotting
#  May  4, 2009 Edward J. Los - Change lc_scatter_flags_local.csh to plot_lightcurves.csh

# run lightcurve script on a bunch of source lists

if ($#argv == 2) then
  set qualifier = ""
  set qualcmd = ""
else 
  if ($#argv == 3) then
      set qualifier = `echo "_$3"`
      set qualcmd = "-q $qualifier"
  else
     echo "usage: run_lc.csh <target catalog> <plate list> [<qualifier>]"
     echo "       target catalog is a list of catalogs, ending in _REF.db"
     echo "       plate list      is the list of mosaics to process"
     echo "       qualifier is a string appended to output filenames"
     echo "For example ./run_lc M44 good_plates_local.list"
      exit
  endif
endif



set imagelist = $2
set lightcurvedir = $DASCH_LIGHTCURVES
set scripts = $DASCH_SCRIPTS

echo "Number of arguments $#argv"
echo "cataloglist is $1"
echo "imagelist is $imagelist"
echo "lightcurvedir is $lightcurvedir"

foreach set ( $1 )
  date
  cd $scripts
  pwd
  echo "index -mb ${set}_REF.db REF"
  index -mb ${set}_REF.db REF 

  cd $lightcurvedir
  if ( ! -e $set ) mkdir $set 
  cd ${lightcurvedir}/${set}
  pwd
echo "${scripts}/extract_lightcurves  -v $qualcmd -f -l ${scripts}/${imagelist} -t ${scripts}/${set}_REF.db  -d ${lightcurvedir}/${set} -b N120013341 -b N120013339 "
  ${scripts}/extract_lightcurves  $qualcmd -f -l ${scripts}/${imagelist} -t ${scripts}/${set}_REF.db -d ${lightcurvedir}/${set}  -b N120013341 -b N120013339
echo "${scripts}/extract_lightcurves -g -v -r $qualcmd -f -l ${scripts}/${imagelist} -t ${scripts}/${set}_REF.db  -d ${lightcurvedir}/${set} -b N120013341 -b N120013339 "
  ${scripts}/extract_lightcurves -g -r $qualcmd -f -l ${scripts}/${imagelist} -t ${scripts}/${set}_REF.db -d ${lightcurvedir}/${set}  -b N120013341 -b N120013339
  ${scripts}/plot_lightcurves.csh

  cd ..
end
date


