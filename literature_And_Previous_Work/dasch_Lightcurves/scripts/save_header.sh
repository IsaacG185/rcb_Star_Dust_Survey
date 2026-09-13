#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Save the headers of the mosaic in $DASCH_HEADERS. No-op for solnum > 0.
#
# Inputs:
#
# - $RAID/ExposureData/Mosaics/${series}/${id}_${mosnum}/${series}${id}_${mosnum}_16${rot}ww.fit
# - $RAID/ExposureData/Mosaics/${series}/${id}_${mosnum}/${series}${id}_${mosnum}_01${rot}ww_tnx.fit
#
# Outputs:
#
# - $DASCH_HEADERS/${series}${id}_${mosnum}_16${rot}ww.hdr
# - $DASCH_HEADERS/${series}${id}_${mosnum}_01${rot}ww_tnx.hdr
#
# No database updates.
#
# Invokes:
#
# - getdirectory
#
# External dependencies:
#
# - apps (FitsHeader)
#
# Variables used:
#
# - $DASCH_HEADERS

if [ $# -ne 2 ] ; then
    echo "usage: $0 <solutionNumber> <list>"
    exit 1
fi

solutionNumber=$1
imagelist=$2

binning=16

if [ $solutionNumber -ne 0 ] ; then
    exit
fi

counter=0
echo "Processing $imagelist"
date

for plate in $(cat $imagelist) ; do
  ((counter=$counter + 1))
  solvedFlag=0
  directory=$(getdirectory $plate)

  # save the thumbnail header

  thumbnail=$(echo "$plate.fit" |sed 's/_01w/_16w/g' |sed 's/_01r/_16r/g')
  filename=${directory}/${thumbnail}
  header=$(echo "$plate.hdr" |sed 's/_01w/_16w/g' |sed 's/_01r/_16r/g')
  headername=${DASCH_HEADERS}/${header}

  if [ -e $filename ] ; then
    if [ ! -e $headername ] ; then
       FitsHeader -i $filename -o $headername -e 0
    fi
  else
    echo "ERROR: $filename not found"
  fi

  # now save the full mosaic header

  filename=${directory}/${plate}_tnx.fit
  headername=${DASCH_HEADERS}/${plate}_tnx.hdr

  if [ -e $filename ] ; then
    FitsHeader -i $filename -o $headername -e 0
  else
    filename=${directory}/${plate}.fit
    headername=${DASCH_HEADERS}/${plate}.hdr

    if [ -e $filename ] ; then
      echo "Found $filename"
      if [ -e $headername ] ; then
         echo "Found $headername"
      else
         FitsHeader -i $filename -o $headername -e 0
      fi
    else
      echo "ERROR: $filename not found"
    fi
  fi
done

echo "$(date) save_header for $counter mosaics"

#  Sep 22, 2009 Edward J. Los - initial version
#  Sep 29, 2009 Edward J. Los - add multiple exposure support
#  Sep 29, 2010 Edward J. Los - unconditionally save tnx headers because of new PV coefficients
