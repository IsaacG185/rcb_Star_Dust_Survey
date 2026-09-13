#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Look for Pickering Wedge Images
#
# TODO: annotate inputs/outputs
#
# Invokes:
#
# - filter_wedge
# - getlocation
# - initgmt
#
# External dependencies:
#
# - Starbase (column, index, jointable)
# - psxy, ps2pdf
#
# Variables used:
#
# - $DASCH_MATCH
# - $DASCH_SCRIPTS

if [ $# -ne 2 ] ; then
    echo "usage: $0 <solutionNumber> <list>"
    exit 1
fi

solutionNumber=$1
imagelist=$2

binning=16
doplots=0
plotstr=""
useRaDec=0

if [ $solutionNumber -eq 0 ] ; then
    solutionString=""
else
    solutionString="_s$solutionNumber"
fi

if [ $doplots -eq 1 ] ; then
    # TODO probably won't work in modernized pipeline
    echo "source ${DASCH_SCRIPTS}/initgmt"
    source ${DASCH_SCRIPTS}/initgmt
    plotstr=" -p "
fi

date

echo "Processing $imagelist"

cd "$DASCH_MATCH"

for plate in $(cat $imagelist) ; do
    series=$(getlocation -s $plate |awk '{print $1}')

    if [ $series == "i" -o $series == "mc" -o $series == "b" -o $series == "dsy" ] ; then
        inputname=${plate}${solutionString}_tnx.db
        outputname=${plate}${solutionString}_wedge.db
        rm -f $outputname

        backgroundname=${plate}_background.db
        if [ -e $backgroundname ] ; then
            backgroundparam="-c $backgroundname"
        else
            backgroundparam=""
        fi

        mosaicsize=$(getlocation -a $plate -e $solutionNumber)
        width=$(echo "$mosaicsize" |awk '{print $1}')
        height=$(echo "$mosaicsize" |awk '{print $2}')
        platescale=$(echo "$mosaicsize" |awk '{print $3}')
        if [ -z "$platescale" ] ; then
            # echo "ERROR getlocation failed with args -e $solutionNumber -a $image"
            continue
        fi

        if [ -e $inputname ] ; then
            tempname=${plate}${solutionString}_tnx.tmp
            tempname1=${plate}${solutionString}_tnx.tmp1
            tempname2=${plate}${solutionString}_tnx.tmp2
            tempname3=${plate}${solutionString}_tnx.tmp3
            tempname4=${plate}${solutionString}_tnx.tmp4
            gmt0name=${plate}${solutionString}_gmt0.txt
            gmt0ps=${plate}${solutionString}_gmt0.ps
            gmt1name=${plate}${solutionString}_gmt1.txt
            gmt1ps=${plate}${solutionString}_gmt1.ps
            gmt2name=${plate}${solutionString}_gmt2.txt
            gmt2ps=${plate}${solutionString}_gmt2.ps
            gmt3name=${plate}${solutionString}_gmt3.txt
            gmt3ps=${plate}${solutionString}_gmt3.ps
            gmt4name=${plate}${solutionString}_gmt4.txt
            gmt4ps=${plate}${solutionString}_gmt4.ps
            gmt5name=${plate}${solutionString}_gmt5.txt
            gmt5ps=${plate}${solutionString}_gmt5.ps
            gmt6name=${plate}${solutionString}_gmt6.txt
            gmt6ps=${plate}${solutionString}_gmt6.ps

            if [ $useRaDec -eq 1 ] ; then
                column -i ${inputname} NUMBER MAG_ISO ra dec FWHM_WORLD AFLAGS X_IMAGE Y_IMAGE >$tempname
                index -mb -n $tempname dec
            else
                column -i ${inputname} NUMBER MAG_ISO FWHM_WORLD AFLAGS X_IMAGE Y_IMAGE >$tempname
                index -mb -n $tempname Y_IMAGE
            fi

            set -x
            filter_wedge \
                $backgroundparam \
                $plotstr \
                -e $solutionNumber \
                -w $width \
                -s $platescale \
                -h $height \
                -r ${plate}${solutionString} \
                -i $tempname \
                -o $outputname
            set +x

            if [ $doplots -eq 1 ] ; then
                psxy $gmt0name -R-350/350/-350/350 -Y4i -JX6i/5.5i -P -B100:"RA error (arcsec)":/100:"DEC error (arcsec)":WSne -Sc0.03  >$gmt0ps
                ps2pdf $gmt0ps
                psxy $gmt1name -R0/360/0/400 -Y4i -JX6i/5.5i -P -B100:"Bin Angle (degrees)":/100:"Bin Count":WSne -W0.5p  >$gmt1ps
                ps2pdf $gmt1ps

                xlimtext=$(minmax $gmt4name |awk '{print $5}' |sed 's/<//g' |sed 's/>//g')
                ylimtext=$(minmax $gmt4name |awk '{print $6}' |sed 's/<//g' |sed 's/>//g')

                psxy $gmt4name -R${xlimtext}/${ylimtext} -Y4i -JX6i/5.5i -P -K -B5:"RA error (arcsec)":/5:"DEC error (arcsec)":WSne -Sx0.2  >$gmt2ps
                ps2pdf $gmt2ps

                psxy $gmt6name -R-20/-6/0/10 -Y4i -JX6i/5.5i -P -K -B1:"MAG_ISO (pri)":/1:"MAG_ISO (sec) - MAG_ISO (pri) ":WSne -W2.0p  >$gmt3ps
                psxy $gmt5name -R -JX  -P -B  -K -O -Sx0.2 >>$gmt3ps
                psxy $gmt3name -R -JX  -P -B  -K -O -Sp >>$gmt3ps

                ps2pdf $gmt3ps
            fi

            if [ $useRaDec -eq 1 ] ; then
                rm ${tempname}.dec.b
            else
                rm ${tempname}.Y_IMAGE.b
            fi

            # if we have an output file, replace the AFLAGS column of the original sextractor file with it
            if [ -e $outputname ] ; then
                tnxheader=$(column -h -i $inputname |head -n 1 |sed 's/\t/ /g')
                index -n -mb $inputname NUMBER
                column -i $outputname NUMBER AFLAGS >$tempname2
                index -n -mb $tempname2 NUMBER
                jointable -j NUMBER -n $inputname $tempname2 >$tempname3
                column -a -i $tempname3 AFLAGS |compute 'AFLAGS = AFLAGS_2' |column $tnxheader >$tempname4
                mv $tempname4 $inputname
                rm ${inputname}.NUMBER.b
                rm ${tempname2}.NUMBER.b
            fi

            rm $tempname
            rm -f $tempname1
            rm -f $tempname2
            rm -f $tempname3
            rm -f $tempname4
        else
            echo "ERROR: Can not find input file $inputname"
        fi
    fi
done

#  May  2, 2008 Edward J. Los
#  May 20, 2008 Edward J. Los - use X_IMAGE and Y_IMAGE instead of RA an DEC
#  May 29, 2008 Edward J. Los - add b plate support
#  Jun 10, 2008 Edward J. Los - correct jointable numeric bug.
#  Jan 27, 2008 Edward J. Los - split FLAGS into AFLAGS and BFLAGS
#  Sep 29, 2009 Edward J. Los - add multiple exposure support
#  Dec 15, 2009 Edward J. Los - Support flagging of Pickering Wedge plates
#  Feb  4, 2010 Edward J. Los - Look for dsy ghost images
#  Jul 23, 2012 Edward J. Los - Add high background object support
