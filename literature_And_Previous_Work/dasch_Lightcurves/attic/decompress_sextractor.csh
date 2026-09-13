#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Dec 28, 2015 Edward J. Los - Compress Sextractor output files
#  
#  In Dec, 2015 a test was conducted on 154 Sextractor files showing that bzip2 provided superior
#  compression over gzip with only 4% more CPU time.
#
# gzip    29071 MB -> 11277 MB 2.57:1 compression  0.7039 hours gzip1.log    (Limited by CPU)
# gunzip  11277 MB -> 29071 MB        expansion    0.1500 hours gunzip2.log  (Limited by disk)
# bzip2   29071 MB ->  8910 MB 3.26:1 compression  0.7319 hours bzip23.log   (Limited by CPU)
# bunzip2  8910 MB -> 29071 MB        expansion    0.3402 hours bunzip24.log (Limited by CPU)
# gzip    29071 MB -> 11277 MB        compression  0.7022 hours gzip5.log    (Limited by CPU)
#

if ($#argv != 2) then
    echo "usage: compress_sextractor.csh  solutionNumber list"
    exit
endif
set solutionNumber = $1
set imagelist = $2
#
# IMPORTANT NOTE: force_update needs to be set to 1 if sextractor parameters change
#
set force_update = 0
set scripts = $DASCH_SCRIPTS
set aper_arcsec = 10
set r_match = 20
set blendUpdateFlag = ""

set listfile = $imagelist:t
echo "Processing ${DASCH_SCRIPTS}/$listfile"
foreach plate ( `cat ${DASCH_SCRIPTS}/$listfile` ) 
  if ($solutionNumber == 0) then 
    #processing of the first exposure

    set tnximage = ${plate}_tnx
    
    if (-e ${DASCH_MATCH}/${tnximage}.db.bz2) then
      if (-e ${DASCH_MATCH}/${tnximage}.db) then
         echo "ERROR: both ${DASCH_MATCH}/${tnximage}.db and ${DASCH_MATCH}/${tnximage}.db.bz2 exist"
      else 
        echo "Decompressing $tnximage"
        date
        bunzip2 ${DASCH_MATCH}/${tnximage}.db.bz2
        date
      endif
    else
      if (-e ${DASCH_MATCH}/${plate}.db.bz2) then
        if (-e ${DASCH_MATCH}/${plate}.db) then
          echo "ERROR: both ${DASCH_MATCH}/${plate}.db and ${DASCH_MATCH}/${plate}.db.bz2 exist"
        else 
          echo "Decompressing $plate"
          date
          bunzip2 ${DASCH_MATCH}/${plate}.db.bz2
          date
        endif
      else 
        echo "ERROR: File ${DASCH_MATCH}/${plate}.db.bz2 not found"
      endif
    endif
  else
    echo "ERROR: only solution 0 is valid for compress_sextractor using $plate";
  endif
end
