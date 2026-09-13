#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Extract sources from a SCAMPed mosaic using SExtractor.
#
# Inputs:
#
# - if solnum = 0:
#   - $RAID/ExposureData/Mosaics/${series}/${id}_${mosnum}/${series}${id}_${mosnum}_01${rot}ww_tnx.fit
# - else:
#   - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${prev_sol_tag}_tnx.db
#   - $DASCH_HEADERS/${series}${id}_${mosnum}_01${rot}ww_s${solnum}_tnx.hdr
#
# Outputs:
#
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${sol_tag}_tnx.db
#
# No database updates.
#
# Invokes:
#
# - getdirectory
# - getlocation
# - update_sextractor
#
# External dependencies:
#
# - Funtools (funhead)
# - SExtractor (sex, sextotable)
# - Starbase (column, etc)
# - Wcstools (xy2sky)
#
# Variables used:
#
# - $DASCH_MATCH
# - $DASCH_SCRIPTS
#
# TODO: redundancy with `run_sextractor.sh`

secondpass=0

if [ $# -eq 3 ] ; then
  if [ $3 -eq 2 ] ; then
    #echo "second pass active"
    secondpass=1
  else
    echo "usage: run_sextractor_second.csh  solutionNumber list pass"
    exit 1
  fi
elif [ $# -ne 2 ] ; then
  echo "usage: run_sextractor_second.csh  solutionNumber list pass"
  exit 1
fi

solutionNumber=$1
imagelist=$2

force_update=0 # set to 1 if sextractor parameters change
aper_arcsec=10
r_match=20
blendUpdateFlag="-u"

echo "Processing ${imagelist}"

for plate in $(cat ${imagelist}) ; do
  if [ $secondpass -eq 1 ] ; then
    #exit quietly if not a patrol plate
    getlocation -n $plate
    if [ $? -ne 0 ] ; then
      #echo "run_scamp.csh not patrol: $plate"
      continue
    fi
  fi

  echo $plate
  directory=$(getdirectory $plate)
  fits_image=${directory}/${plate}_tnx.fit
  echo Begin Process $(date)
  echo fits_image $fits_image
  echo image ${plate}_tnx
  echo phot_aperture $aper_arcsec
  echo match_radius $r_match

  cd $DASCH_MATCH

  script_output=$(getlocation $plate -e $solutionNumber)
  scale=$(echo "$script_output" |awk '{print $4}')
  if [ -z "$scale" ] ; then
    # echo "ERROR getlocation failed with args $plate -e $solutionNumber"
    continue
  fi

  mosaicsize=$(getlocation -a ${plate}_tnx -e $solutionNumber)
  leftMargin=$(echo "$mosaicsize" |awk '{print $5}')
  rightMargin=$(echo "$mosaicsize" |awk '{print $6}')
  bottomMargin=$(echo "$mosaicsize" |awk '{print $7}')
  topMargin=$(echo "$mosaicsize" |awk '{print $8}')

  echo "Pixel Scale: " $scale
  # TODO(?): This value should increase with arcsec/mm.
  aper_pix=$(echo $aper_arcsec $scale |awk '{print $1/$2}')
  echo aper_pix is $aper_pix

  if [ $solutionNumber -eq 0 ] ; then
    #processing of the first exposure

    baseimage=$(echo ${plate}_tnx |sed 's/_tnx//g')
    tnximage=${baseimage}_tnx

    # This script used to have logic to check if a SExtractor database from "WW"
    # processing was present; if so, the update_sextractor tool was used rather
    # than full reprocessing. However, this approach yields different results
    # than what you get if you actually re-extract the TNX file, and while the
    # differences are small, they do affect pipeline results. So to obtain
    # consistency between first-time processing and reprocessing, we now always
    # re-SExtract.

    echo "soln 0 margins $leftMargin $rightMargin $bottomMargin $topMargin for ${plate}_tnx"
    table=${plate}_tnx.db

    set -x
    sex \
      -c $DASCH_SCRIPTS/Sextractor/DASCH.config \
      $fits_image \
      -PARAMETERS_NAME $DASCH_SCRIPTS/Sextractor/DASCH.param \
      -FILTER_NAME $DASCH_SCRIPTS/Sextractor/default.conv \
      -PHOT_APERTURES $aper_pix \
      -CATALOG_NAME ${DASCH_MATCH}/${plate}_tnx.cat

    sextotable \
      <${DASCH_MATCH}/${plate}_tnx.cat \
      |sed -e 's/FLAGS/BFLAGS/g' -e 's/ALPHA_J2000/ra/g' -e 's/DELTA_J2000/dec/g' \
      >${plate}_tnxcatalog.tmp

    rm ${DASCH_MATCH}/${plate}_tnx.cat
    set +x

    # determine the plate center from header keywords
    xc=$(funhead $fits_image |grep NAXIS1 |awk '{print $3/2}')
    yc=$(funhead $fits_image |grep NAXIS2 |awk '{print $3/2}')
    rac=$(xy2sky -jd $fits_image $xc $yc |awk '{print $1}')
    decc=$(xy2sky -jd $fits_image $xc $yc |awk '{print $2}')

    echo Plate Center $rac $decc
    set -x

    column \
      -i ${plate}_tnxcatalog.tmp \
      -a plate_dra plate_ddec plate_dist \
      |compute  \
      "plate_ddec=(dec - $decc); plate_dra=((ra - $rac)*(cos(((dec + $decc)/2)/57.29577951))); plate_dist=sqrt(plate_dra^2+plate_ddec^2)" \
      >${plate}_tnxcatalog2.tmp

    update_sextractor \
      -e $solutionNumber \
      -s \
      $blendUpdateFlag \
      -l $leftMargin \
      -r $rightMargin \
      -b $bottomMargin \
      -t $topMargin \
      -i ${plate}_tnxcatalog2.tmp \
      -o $table \
      -m $fits_image

    rm ${plate}_tnxcatalog.tmp
    rm ${plate}_tnxcatalog2.tmp
    set +x
  else
    #processing of additional exposures
    ((prevNumber=$solutionNumber - 1))
    solutionString="_s$solutionNumber"

    if [ $prevNumber -eq 0 ] ; then
      prevString=""
    else
      prevString="_s$prevNumber"
    fi

    update_source=${DASCH_MATCH}/${plate}${prevString}_tnx.db
    table=${plate}${solutionString}_tnx.db
    fits_image=${DASCH_HEADERS}/${plate}${solutionString}_tnx.hdr

    if [ ! -e $update_source ] ; then
      update_source=${DASCH_MATCH}/${plate}${prevString}_tnx.db
    fi

    testfile=$(head -n 1 $update_source |grep ERRA_IMAGE)
    if [ -z "$testfile" ] ; then
      echo "ERROR: ERRA_IMAGE (3) not found in $update_source"
      continue
    fi

    rm -f $table

    if [ -e $update_source -a -e $fits_image ] ; then
      set -x
      update_sextractor \
        $blendUpdateFlag \
        -e $solutionNumber \
        -l $leftMargin \
        -r $rightMargin \
        -b $bottomMargin \
        -t $topMargin \
        -i $update_source \
        -o $table \
        -m $fits_image
      set +x
    else
      if [ ! -e $update_source ] ; then
        echo "ERROR: run_sextractor $update_source is missing"
      fi
      if [ ! -e $fits_image ] ; then
        echo "ERROR: run_sextractor $fits_image is missing"
      fi
    fi
  fi

  if [ ! -s ${DASCH_MATCH}/${table} ] ; then
    echo "ERROR File ${DASCH_MATCH}/${table} has zero size, deleting result"
    rm ${DASCH_MATCH}/${table}
  fi
done

#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#  Nov  7, 2007 Edward J. Los - Operate on the mosaics in their working directory
#  Sep 29, 2009 Edward J. Los - add multiple exposure support
#                 remove reference to dasch_sextractor.csh
#  Sep 20, 2010 Edward J. Los - Add new parameters for scamp support
#  Jul  3, 2012 Edward J. Los - Add high background object support
#  Jul 23, 2012 Edward J. Los - Move background support to filter_wedge.csh
#  Jul 23, 2013 Edward J. Los   Exit quietly if second pass and the plate is not a patrol plate
#  Mar 30, 2020 Edward J. Los - Delete the result if it has zero bytes
#  Apr 16, 2020 Edward J. Los - Correct previous change
