#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Inputs:
#
# - if solnum = 0:
#   - $RAID/ExposureData/Mosaics/${series}/${id}_${mosnum}/${series}${id}_${mosnum}_01${rot}ww.fit
# - else:
#   - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${prev_sol_tag}_tnx.db
#   - $DASCH_HEADERS/${series}${id}_${mosnum}_01${rot}ww_s${solnum}.hdr
#   - $DASCH_ASTROMETRY/${series}${id}_${mosnum}_01${rot}ww_s${solnum}_none.db
#
# Outputs:
#
# - $DASCH_MATCH/${series}${id}_${mosnum}_01${rot}ww${sol_tag}.db
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
# - $DASCH_ASTROMETRY
# - $DASCH_MATCH
# - $DASCH_SCRIPTS

if [ $# -ne 2 ] ; then
    echo "usage: $0 <solutionNumber> <list>"
    exit 1
fi

solutionNumber=$1
imagelist=$2

force_update=0 # set to 1 if sextractor parameters change
aper_arcsec=10
blendUpdateFlag=""

echo "Processing ${imagelist}"

for plate in $(cat $imagelist) ; do
  echo $plate
  directory=$(getdirectory $plate)
  fits_image=${directory}/${plate}.fit
  echo Begin Process $(date)
  echo fits_image $fits_image
  echo image $plate
  echo phot_aperture $aper_arcsec

  cd $DASCH_MATCH

  script_output=$(getlocation $plate -e $solutionNumber)
  scale=$(echo "$script_output" |awk '{print $4}')
  if [ -z "$scale" ] ; then
    # echo "ERROR getlocation failed with args $plate -e $solutionNumber"
    continue
  fi

  mosaicsize=$(getlocation -a $plate -e $solutionNumber)
  leftMargin=$(echo "$mosaicsize" |awk '{print $5}')
  rightMargin=$(echo "$mosaicsize" |awk '{print $6}')
  bottomMargin=$(echo "$mosaicsize" |awk '{print $7}')
  topMargin=$(echo "$mosaicsize" |awk '{print $8}')

  echo "Pixel Scale: " $scale
  aper_pix=$(echo $aper_arcsec $scale |awk '{print $1/$2}')
  echo aper_pix is $aper_pix

  if [ $solutionNumber -eq 0 ] ; then
    #processing of the first exposure

    baseimage=$(echo $plate |sed 's/_tnx//g')
    tnximage=${baseimage}_tnx

    if [ $baseimage = $plate ] ; then
       if [ $force_update -ne 0 ] ; then
         updatedb=0
       else
         if [ -e ${DASCH_MATCH}/${tnximage}.db ] ; then
            if [ ! -s ${DASCH_MATCH}/${tnximage}.db ] ; then
               echo "File ${DASCH_MATCH}/${tnximage}.db has zero size, repeating sextractor"
               updatedb=0
            else
               testfile=$(head -n 1 ${DASCH_MATCH}/${tnximage}.db |grep ERRA_IMAGE)
               if [ -z "$testfile" ] ; then
                 echo "WARNING: ERRA_IMAGE (1) not found, repeating sextractor"
                 updatedb=0
               else
                 echo "WARNING: force_update is set to zero assuming that sextractor parameters have not changed"
                 updatedb=1
                 update_source=${DASCH_MATCH}/${tnximage}.db
               fi
            fi
         else
            echo "File ${DASCH_MATCH}/${tnximage}.db not found, repeating sextractor"
            updatedb=0
         fi
       fi
    else
      echo "ERROR: base image $baseimage does not equal image $plate"
      continue
    fi

    echo "Update Database $updatedb margins $leftMargin $rightMargin $bottomMargin $topMargin for $plate"
    table=${plate}.db

    if [ $updatedb -eq 0 ] ; then
      set -x

      sex \
        -c ${DASCH_SCRIPTS}/Sextractor/DASCH.config \
        $fits_image \
        -PARAMETERS_NAME ${DASCH_SCRIPTS}/Sextractor/DASCH.param \
        -FILTER_NAME ${DASCH_SCRIPTS}/Sextractor/default.conv \
        -PHOT_APERTURES $aper_pix \
        -CATALOG_NAME ${DASCH_MATCH}/${plate}.cat

      sextotable \
        <${DASCH_MATCH}/${plate}.cat \
        |sed -e 's/FLAGS/BFLAGS/g' -e 's/ALPHA_J2000/ra/g' -e 's/DELTA_J2000/dec/g' \
        >${plate}catalog.tmp

      rm ${DASCH_MATCH}/${plate}.cat
      set +x

      # determine the plate center from header keywords
      xc=$(funhead $fits_image |grep NAXIS1 |awk '{print $3/2}')
      yc=$(funhead $fits_image |grep NAXIS2 |awk '{print $3/2}')
      rac=$(xy2sky -jd $fits_image $xc $yc |awk '{print $1}')
      decc=$(xy2sky -jd $fits_image $xc $yc |awk '{print $2}')

      echo Plate Center $rac $decc
      set -x

      column \
        -i ${plate}catalog.tmp \
        -a plate_dra plate_ddec plate_dist \
        |compute \
          "plate_ddec=(dec - $decc); plate_dra=((ra - $rac)*(cos(((dec + $decc)/2)/57.29577951))); plate_dist=sqrt(plate_dra^2+plate_ddec^2)" \
        >${plate}catalog2.tmp

      update_sextractor \
        -e $solutionNumber \
        -s \
        $blendUpdateFlag \
        -l $leftMargin \
        -r $rightMargin \
        -b $bottomMargin \
        -t $topMargin \
        -i ${plate}catalog2.tmp \
        -o $table \
        -m $fits_image

      rm ${plate}catalog.tmp
      rm ${plate}catalog2.tmp
      set +x
    else
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
    fi
  else
    #processing of additional exposures
    ((prevNumber=$solutionNumber - 1))
    solutionString="_s$solutionNumber"

    if [ $prevNumber -eq 0 ] ; then
      prevString=""
    else
      prevString="_s$prevNumber"
    fi

    update_source=${DASCH_MATCH}/${plate}${prevString}.db
    table=${plate}${solutionString}.db
    fits_image=${DASCH_HEADERS}/${plate}${solutionString}.hdr
    none_image=${DASCH_ASTROMETRY}/${plate}${solutionString}_none.db

    if [ ! -e $update_source ] ; then
      update_source=${DASCH_MATCH}/${plate}${prevString}_tnx.db
    fi

    testfile=$(head -n 1 $update_source |grep ERRA_IMAGE)
    if [ -z "$testfile" ] ; then
      echo "ERROR: ERRA_IMAGE (2) not found in $update_source"
      continue
    fi

    rm -f $table

    if [ -e $update_source -a -e $fits_image -a -e $none_image ] ; then
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
        -m $fits_image \
        -n $none_image
      set +x
    else
      if [ ! -e $none_image ] ; then
        echo "ERROR: run_sextractor $none_image is missing"
      elif [ ! -e $update_source ] ; then
        echo "ERROR: run_sextractor $update_source is missing"
      elif [ ! -e $fits_image ] ; then
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
#                               Remove reference to dasch_sextractor.csh
#  Sep 20, 2010 Edward J. Los - Add new parameters for scamp support
#  Apr 16, 2020 Edward J. Los - Delete the result if it has zero bytes
