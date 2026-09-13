#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#  noise_repair.csh
#
#  run the noise repair algorithm on the file list
#
#  May   2, 2018 Edward J. Los - perform image repair processing for a row noise problem first manifested in Januray, 2018
#  May   4, 2018 Edward J. Los - Add the noise database file and optional deletion
#  May  11, 2018 Edward J. Los - Add the -BACKGROUND (nobackgroundfile) check file
#
echo "ERROR: noise_repair.csh is deprecated: Use ProcNoiseRepair"
exit

if ($#argv == 2) then
    set deleteflag = 0
else 
    if ($#argv == 3) then
       if ($3 != "delete") then
            echo "usage: third parameter [delete]"
            exit
         else
       set deleteflag = 1
       endif
    else 
       echo "usage: noise_repair.csh  <list>  <db file> [delete]"
       exit
    endif
endif
set imagelist = $1
set dbsummary = $2
set scripts = $DASCH_SCRIPTS
set counter = 1

echo "Processing $imagelist Copy location $DASCH_COMPLETED"
foreach inputfile ( `cat $imagelist` ) 

  set directory = $inputfile:h
  set inputname = $inputfile:t
  set inputroot = $inputname:r
  set inputext  = $inputfile:e

  set datestr = `date`
  echo "noise_repair begin $datestr $counter $directory/$inputname"
  @ counter = $counter + 1
  
  set invertfile = $directory/${inputroot}_inv.{$inputext}

  set pixelnoisefile = $directory/${inputroot}_pixelnoise.{$inputext}
  set nobackgroundfile = $directory/${inputroot}_nobackground.{$inputext}
  set maskfile = $directory/${inputroot}_mask.{$inputext}
  set outfile = $directory/${inputroot}_out.{$inputext}

  set catfile    = $directory/${inputroot}.cat
  set databasefile    = $directory/${inputroot}.db


  echo "FitsMath -i -o $invertfile $inputfile"
  FitsMath -i -o $invertfile $inputfile
  echo "sex -c  ${scripts}/Sextractor/DASCH.config $invertfile -PARAMETERS_NAME ${scripts}/Sextractor/DASCH.param -FILTER_NAME  ${scripts}/Sextractor/default.conv -PHOT_APERTURES 1.67092  -CATALOG_NAME $catfile -CHECKIMAGE_NAME $pixelnoisefile,$maskfile,$nobackgroundfile -CHECKIMAGE_TYPE -OBJECTS,-MASK,-BACKGROUND"
  sex -c  ${scripts}/Sextractor/DASCH.config $invertfile -PARAMETERS_NAME ${scripts}/Sextractor/DASCH.param -FILTER_NAME  ${scripts}/Sextractor/default.conv -PHOT_APERTURES 1.67092  -CATALOG_NAME $catfile -CHECKIMAGE_NAME $pixelnoisefile,$maskfile,$nobackgroundfile -CHECKIMAGE_TYPE -OBJECTS,-MASK,-BACKGROUND

  echo "NoiseRepair -r  $invertfile -i -p $pixelnoisefile -m $maskfile -b $nobackgroundfile -o $outfile -s $databasefile -n $dbsummary"
  NoiseRepair -r  $invertfile -i -p $pixelnoisefile -m $maskfile -b $nobackgroundfile -o $outfile -s $databasefile -n $dbsummary 

  if ($deleteflag == 1) then
    rm $catfile
    rm $invertfile
    rm $outfile
    rm $pixelnoisefile 
    rm $maskfile
    rm $nobackgroundfile
    rm $databasefile
  endif

  set datestr = `date`
  echo "noise_repair done $datestr $counter $outfile"

end
