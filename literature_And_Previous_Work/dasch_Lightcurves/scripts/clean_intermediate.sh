#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Clean up intermediate files after processing is completed
#
# What REMAINS:
#
# - ${DASCH_BINOUTPUT}/${series}${id}_${scannum}_01${rot}ww${qual}_a{1..9}.grid
# - ${DASCH_BINOUTPUT}/${series}${id}_${scannum}_01${rot}ww${qual}_a{1..9}.para
# - ${DASCH_CATALOGALL}/${series}${id}_${scannum}_01${rot}ww${qual}_extinction.db
# - ${DASCH_HEADERS}/${series}${id}_${scannum}_01${rot}ww${qual}_tnx.hdr
# - ${DASCH_HEADERS}/${series}${id}_${scannum}_16${rot}ww${qual}.hdr
# - ${DASCH_INGEST}/${series}${id}_${scannum}_01${rot}ww${qual}_allobjects.db
# - ${DASCH_INGEST}/${series}${id}_${scannum}_01${rot}ww${qual}_dmagcor.db
# - ${DASCH_INGEST}/${series}${id}_${scannum}_01${rot}ww${qual}_dmagcor_copy.db
# - ${DASCH_INGEST}/${series}${id}_${scannum}_01${rot}ww${qual}_magdepcalibrate.db
# - ${DASCH_INGEST}/${series}${id}_${scannum}_01${rot}ww${qual}.out.local.db
# - ${DASCH_INGEST}/${series}${id}_${scannum}_01${rot}ww${qual}.out.spatial_bins.db
# - ${DASCH_MATCH}/${series}${id}_${scannum}_01${rot}ww${qual}_background.db
# - ${DASCH_MATCH}/${series}${id}_${scannum}_01${rot}ww${qual}_defect.db
# - ${DASCH_MATCH}/${series}${id}_${scannum}_01${rot}ww${qual}_drad.db
# - ${DASCH_MATCH}/${series}${id}_${scannum}_01${rot}ww${qual}_summary_defect.db
# - ${DASCH_MATCH}/${series}${id}_${scannum}_01${rot}ww${qual}_tnx.db
# - ${DASCH_MATCH}/${series}${id}_${scannum}_01${rot}ww${qual}_tnx_b.db
# - ${DASCH_MATCH}/${series}${id}_${scannum}_01${rot}ww${qual}_tnx_u.db
#
# Invokes:
#
# - checksize
#
# Variables used:
#
# - $DASCH_BINOUTPUT
# - $DASCH_CATALOGALL
# - $DASCH_CATALOGBIN
# - $DASCH_INGEST
# - $DASCH_MATCH
# - $DASCH_SCAMPIMAGES
# - $DASCH_SCRATCH

solutionNumber=$1
imagelist=$2
scripts=$DASCH_SCRIPTS

if [ $# -eq 2 ] ; then
  qualifier=""
elif [ $# -eq 3 ] ; then
  qualifier="_$3"
else
  echo "usage: $0 <solutionNumber> <list> [qualifier]"
  echo "  solutionNumber is the multiple exposure iteration"
  echo "  list      is the list of mosaics to process"
  echo "  qualifier is a string appended to output filenames"
  exit 1
fi

if [ $solutionNumber -eq 0 ] ; then
  solutionString=""
else
  solutionString="_s$solutionNumber"
fi

for plate in $(cat $imagelist) ; do
  rm -f ${DASCH_BINOUTPUT}/${plate}${solutionString}${qualifier}.annular9_marker.txt

  for bin in 1 2 3 4 5 6 7 8 9 ; do
    rm -f ${DASCH_BINOUTPUT}/${plate}${solutionString}${qualifier}_a${bin}.out
    rm -f ${DASCH_CATALOGBIN}/${plate}${solutionString}${qualifier}_a${bin}.db
    rm -f ${DASCH_CATALOGBIN}/${plate}${solutionString}${qualifier}_a${bin}.colorterm.txt
  done

  rm -f ${DASCH_CATALOGBIN}/${plate}${solutionString}${qualifier}.colorterm_marker.txt
  rm -f ${DASCH_CATALOGBIN}/${plate}${solutionString}${qualifier}_9bins_a.eps

  rm -f ${DASCH_SCRATCH}/${plate}${solutionString}${qualifier}.divide_annul9.txt

  rm -f ${DASCH_CATALOGALL}/${plate}${solutionString}${qualifier}.db
  rm -f ${DASCH_CATALOGALL}/${plate}${solutionString}${qualifier}_a0.db

  rm -f ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.db
  rm -f ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.db.dec.i
  rm -f ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.db.REF.i
  rm -f ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.db.ra.b
  rm -f ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.local.db.dec-i
  rm -f ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.local.db.REF-i
  rm -f ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out.local.db.ra-b
  rm -f ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.para
  rm -f ${DASCH_INGEST}/${plate}${solutionString}${qualifier}.out

  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}.coo
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}.db
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}.log
  # SCAMP intermediate files
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_scamp.hdr
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}.tmp
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}.sexcat_filtered.fits
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}.refcat.ucac4.db
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}.refcat.ucac4.fits
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}.refcat.ucac5.db
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}.refcat.ucac5.fits
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}.sexcat_filtered.head
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_scamp_results.xml
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_scamp.refcat_used_1.fits
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}.sexcat_filtered.head_cor
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}.tmp2
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_scamp.db
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_scamp.sexcat_filtered.fits
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_match.fits
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_match_tnx.xml
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_match_tnx.fits
  rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_match_scamp.fits

  rm -f ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}.sexcat_filtered.fits.rms
  rm -f ${DASCH_SCAMPIMAGES}/${plate}${solutionString}${qualifier}_scamp.sexcat_filtered.fits.rms

  for bin in 0 1 2 3 ; do
    rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_gmt${bin}.pdf
    rm -f ${DASCH_MATCH}/${plate}${solutionString}${qualifier}_gmt${bin}.ps
  done

  rm -f ${DASCH_MATCH}/match_${plate}${solutionString}${qualifier}_u.db

  rm -f ${DASCH_MATCH}/match_${plate}${solutionString}_tnx${qualifier}_lim.db
  rm -f ${DASCH_MATCH}/match_${plate}${solutionString}_tnx${qualifier}_lim1.db
  rm -f ${DASCH_MATCH}/match_${plate}${solutionString}_tnx${qualifier}_lim2.db
  rm -f ${DASCH_MATCH}/match_${plate}${solutionString}_tnx${qualifier}_lim3.db

  # Get rid of zero length files or stub headers

  filename=${DASCH_CATALOGALL}/${plate}${solutionString}${qualifier}.db
  checksize -f $filename
  if [ $? -gt 1 ] ; then
    rm $filename
  fi

  filename=${DASCH_CATALOGALL}/${plate}${solutionString}${qualifier}_extinction.db
  checksize -f $filename
  if [ $? -gt 1 ] ; then
    rm $filename
  fi

  filename=${DASCH_MATCH}/match_${plate}${solutionString}_tnx${qualifier}_b.db
  checksize -f $filename
  if [ $? -gt 1 ] ; then
    rm $filename
  fi

  filename=${DASCH_MATCH}/match_${plate}${solutionString}_tnx${qualifier}_u.db
  checksize -f $filename
  if [ $? -gt 1 ] ; then
    rm $filename
  fi

  # Get rid of octave stubs (8 lines long)

  filename=${DASCH_CATALOGALL}/${plate}${solutionString}${qualifier}.db
  checksize -l 9 -f $filename
  if [ $? -gt 1 ] ; then
    rm $filename
  fi
done

#  Jul 28, 2008 Edward J. Los - initial version
#  Jan 12, 2009 Edward J. Los - Remove colorterm files
#  Jul 28, 2009 Edward J. Los - Do not delete graphs
#  Sep  7, 2009 Edward J. Los - Add pass2 support
#  Sep 29, 2009 Edward J. Los - add multiple exposure support
#  Mar 15, 2010 Edward J. Los - Remove marker files
#  Sep 29, 2010 Edward J. Los - Add scamp intermediate files
#  Apr 22, 2010 Edward J. Los - Remove .db file in DASCH_CATALOGALL
