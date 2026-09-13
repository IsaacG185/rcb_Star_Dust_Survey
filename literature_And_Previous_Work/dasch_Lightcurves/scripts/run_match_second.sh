#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Inputs:
#
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx.db
#
# Outputs:
#
# - $DASCH_MATCH/match_${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx${qual}_u.db
# - $DASCH_MATCH/match_${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx${qual}_b.db
# - $DASCH_MATCH/match_${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx${qual}_lim.db
# - $DASCH_MATCH/match_${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx${qual}_lim1.db
# - $DASCH_MATCH/match_${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx${qual}_lim2.db
# - $DASCH_MATCH/match_${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx${qual}_lim3.db
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_drad.db
#
# No database updates.
#
# Invokes:
#
# - filterblended
# - getdirectory
# - getlocation
# - juliandate
# - matchstars
#
# External dependencies:
#
# - Starbase (index, column, matchstars, sorttable)
#
# Variables used:
#
# - $DASCH_CATALOG (implicitly by matchstars)
# - $DASCH_HEADERS
# - $DASCH_MATCH
# - $DASCH_SCRIPTS

solutionNumber=$1
imagelist=$2

if [ $# == 2 ] ; then
  qualifier=""
  qualcmd=""
elif [ $# == 3 ] ; then
  qualifier="_$3"
  qualcmd="-q $3"
else
  echo "usage: $0 <solutionNumber> <list> [qualifier]"
  exit 1
fi

if [ $solutionNumber -eq 0 ] ; then
  solutionString=""
else
  solutionString="_s$solutionNumber"
fi

echo "Processing $imagelist"

for plate in $(cat $imagelist) ; do
  echo $plate
  directory=$(getdirectory $plate)

  if [ $solutionNumber == 0 ] ; then
    filename=${directory}/${plate}${solutionString}_tnx.fit
  else
    filename=${DASCH_HEADERS}/${plate}${solutionString}_tnx.hdr
  fi

  fits_image=$filename
  image=${plate}${solutionString}_tnx

  # NOTE: any change in the search radius should be reflected in
  #       MULTIPLE_EXPOSURE_TOLERANCE of AstrometryWCS.cpp for multiple
  #       exposures
  r_match=20

  series=$(getlocation -s $fits_image |awk '{print $1}')
  script_output=$(getlocation -e $solutionNumber $fits_image)
  plate_scale=$(echo "$script_output" |awk '{print $4}')

  if [ -z "$plate_scale" ] ; then
    # echo "ERROR getlocation failed with args -e $solutionNumber $fits_image"
    continue
  fi

  mosaicsize=$(getlocation -e $solutionNumber -a $fits_image)
  width=$(echo "$mosaicsize" |awk '{print $1}')
  height=$(echo "$mosaicsize" |awk '{print $2}')

  echo Begin Process $(date)
  echo fits_image $fits_image
  echo image $image
  echo qualifier $qualifier
  echo qualcmd $qualcmd

  cd $DASCH_MATCH

  estimate=match_${image}${qualifier}_lim.db
  table=${image}.db
  tmptable=${image}${qualifier}_tmp.db

  # 3. Match Detected Stars with Reference Stars
  rm -f match_${table}${qualifier}_u.db
  rm -f match_${table}${qualifier}_b.db
  rm -f match_${table}${qualifier}_lim.db
  rm -f match_${table}${qualifier}_lim1.db
  rm -f match_${table}${qualifier}_lim2.db
  rm -f match_${table}${qualifier}_lim3.db
  rm -f $estimate

  dradfile=$(echo ${image}${qualifier} |sed 's/_tnx//g')

  #NOTE: for the purposes of calculating proper motion, the solutionNumber == 0 julian date is good enough
  juliandate=$(juliandate $dradfile |awk '{print $1}')
  #echo "Julian Date: $juliandate for $dradfile"

  dradfile=${dradfile}_drad.db
  rm -f $dradfile

  set -ex
  column \
    -i ${table} \
    NUMBER MAG_ISO X_IMAGE Y_IMAGE ra dec FWHM_IMAGE FWHM_WORLD AFLAGS BFLAGS THRESHOLD FLUX_MAX aLength bLength THETA_J2000 \
    >$tmptable
  index -mb -n $tmptable dec

  # You can tell we're close to wrapping up because I'm just hardcoding the path
  # of the global big scratch filesystem here.
  matchtable=$(mktemp -p /n/netscratch/dasch_project/Lab match_${image}${qualifier}.XXXXXX.db)

  date
  matchstars \
    $qualcmd \
    -j $juliandate \
    -w $width \
    -h $height \
    -i $tmptable \
    -e $estimate \
    -o $matchtable \
    -s $plate_scale
  date

  ls -lh $matchtable

  # If matchstars emits no rows, subsequent processing will fail. In that case
  # the file has two lines: the column names and the delimiter hyphens. Note that the file
  # might also be 60 GB large, so we don't want to `wc -l` it willy-nilly

  if [ "$(head "$matchtable" |wc -l |awk '{print $1}')" -lt 3 ] ; then
    echo >&2 "error: no data in matchstars output file $matchtable"
    rm -f \
      $(echo ${matchtable} |sed -e 's/\.db$/_u.db/') \
      $(echo ${matchtable} |sed -e 's/\.db$/_b.db/') \
      $(dirname "$matchtable")/$(basename "$dradfile")
    exit 33
  fi

  # This program creates a temporary file in $TMPDIR which might also be too large for
  # /tmp or /scratch
  TMPDIR=$(dirname $matchtable) sorttable -i $matchtable -n Y_IMAGE >${matchtable}tmp

  filterblended \
    -j $juliandate \
    -f $solutionNumber \
    $qualcmd \
    -t 7200 \
    -s $plate_scale \
    -w $width \
    -h $height \
    -e $estimate \
    ${matchtable}tmp

  # filterblended's output filenames are derived from the input filename. Since
  # the input is a tempfile in our large scratch area, we need to put outputs
  # back into $DASCH_MATCH under the expected name.

  udb_before=$(echo ${matchtable} |sed -e 's/\.db$/_u.db/')
  udb_after="$DASCH_MATCH/match_${image}${qualifier}_u.db"
  mv "$udb_before" "$udb_after"

  bdb_before=$(echo ${matchtable} |sed -e 's/\.db$/_b.db/')
  bdb_after="$DASCH_MATCH/match_${image}${qualifier}_b.db"
  mv "$bdb_before" "$bdb_after"

  drad_before="$(dirname "$matchtable")/$(basename "$dradfile")"
  mv "$drad_before" "$dradfile"

  rm ${matchtable}
  rm ${matchtable}tmp
  rm -f ${tmptable}.dec.b
  rm -f ${tmptable}
  set +ex
done

#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#  Nov  7, 2007 Edward J. Los - Operate on the mosaics in their working directory
#  Nov 23, 2007 Edward J. Los - Support filtering of blended stars
#  Dec  5, 2007 Edward J. Los - Create a new catalog for each plate
#  Aug 18, 2009 Edward J. Los - Add kepler input catalog support
#  Sep 29, 2009 Edward J. Los - add multiple exposure support
#  Oct 13, 2009 Edward J. Los - Merge code from dasch_match.csh
#  Feb 24, 2010 Edward J. Los - Double the maximum timeout
#  May 28, 2010 Edward J. Los - Add an optional region file for visualization
#  Nov 29, 2010 Edward J. Los - Add new SCAMP search radii
#  Jun 27, 2010 Edward J. Los - Prevent use of a stale multiple exposure source list.
#  Aug 27, 2013 Edward J. Los - Remove the match radius, now set by object dimensions
