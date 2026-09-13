#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Convert the match_* files from the catalog match to files to be
# used by the octave scripts.
#
# Inputs:
#
# - ${DASCH_MATCH}/match_${name}${sol_tag}_tnx${qual}_u.db
#
# Outputs:
#
# - $DASCH_CATALOGALL/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}.db
# - $DASCH_CATALOGALL/${series}${id}_${mosnum}_01${rot}ww${sol_tag}${qual}_extinction.db
#
# Database updates: none.
#
# Invokes:
#
# - getdirectory
# - prepare_octave
#
# Variables used:
#
# - $DASCH_CATALOGALL
# - $DASCH_INGEST, optionally
# - $DASCH_MATCH


#   Fields from match_${table}
#   1: REF 2: ra_1   3: dec_1  4: Stdmag  5: color  7: MAG_ISO  9: MAG_APER 17: ISOAREA_IMAGE  39: plate_dist
#   19: X_IMAGE   20: Y_IMAGE
#   41: ddec 42: drad  6: NUMBER  36: FLAGS
#
#   Output Fields: X_WIDTH and Y_WIDTH are the width and height in pixels of the mosaic.
#   1: REF 2: ra_1   3: dec_1  4: Stdmag  5: color  6: MAG_ISO  7: MAG_APER 8: ISOAREA_IMAGE  8: plate_dist
#   9: ddec 10: drad 11: NUMBER 12: FLAGS 13: X_IMAGE 14: Y_IMAGE 15: X_WIDTH  16: Y_WIDTH 16: PLATE_SCALE

solutionNumber=$1
imagelist=$2

if [ $# -eq 2 ] ; then
  qualifier=""
  magnitudefile=""
  qualcmd=""
elif [ $# -eq 3 ] ; then
  qualifier="_$3"
  qualcmd="-q $qualifier"
elif [ $# -eq 4 ] ; then
  qualifier="_$3"
  if [ $qualifier == "_magtest" ] ; then
    magnitudefile=""
    qualcmd=""
  else
    magnitudefile=$4
    qualcmd="-q $qualifier -c ${DASCH_INGEST}/$magnitudefile"
  fi
else
  echo "usage: $0 <solutionNumber> <list> [qualifier] [database]"
  echo "  solutionNumber is the multiple exposure iteration"
  echo "  list      is the list of mosaics to process"
  echo "  qualifier is a string append to output filenames"
  echo "  database  is the recalibration datbase"
  exit 1
fi

if [ $solutionNumber -eq 0 ] ; then
  solutionString=""
else
  solutionString="_s$solutionNumber"
fi


for plate in $(cat $imagelist) ; do
  echo "Plate: ${plate}${solutionString}"

  datapath=$(getdirectory ${plate}${solutionString})
  fits_image=${datapath}/${plate}${solutionString}_tnx.fit
  image=${plate}${solutionString}_tnx
  input_file=match_${image}${qualifier}_u.db

  output_file=$(echo ${DASCH_CATALOGALL}/${image}${qualifier}.db | sed 's/_tnx//g')
  extinction_file=$(echo ${DASCH_CATALOGALL}/${plate}${solutionString}${qualifier}_extinction.db)

  cd $DASCH_MATCH

  set -x
  prepare_octave \
    $qualcmd \
    -s $solutionNumber \
    -m $DASCH_MATCH/${input_file} \
    -o $output_file \
    -e $extinction_file \
    -p ${plate}${solutionString}
  set +x
done

#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#  Nov  7, 2007 Edward J. Los - Operate on the mosaics in their working directory
#  Dec 11, 2007 Edward J. Los - Add FLAGS and NUMBER columns
#  Dec 29, 2007 Edward J. Los - Support recalibrated GSC magnitudes.
#  Jan 18, 2008 Edward J. Los - Include X_IMAGE, Y_IMAGE and the maximum values of these two parameters
#                               Filter out FLAGS >= 65536 on first pass
#  Feb 15, 2008 Edward J. Los - Add the plate scale to the output
#  Mar  5, 2008 Edward J. Los - Add extincation.  Move to a C-language program
#  Mar 15, 2008 Edward J. Los - Use "magtest" as a special exception from recalibration
#  Mar 28, 2008 Edward J. Los - Use the new match directory
#  Aug 20, 2009 Edward J. Los - Add kepler catalog support
#  Sep 29, 2009 Edward J. Los - add multiple exposure support
