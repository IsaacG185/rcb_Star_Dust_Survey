# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Combine all of the colorterm files into a single starbase file
# 
# Mar 3, 2009 Edward J. Los - Initial implementation adapted from ingest_matlab2.csh

set imagelist = $1

if ($#argv == 1) then
  set qualifier = ""
else 
  if ($#argv == 2) then
      set qualifier = `echo "_$2"`
  else
     echo "usage: xingest_colorterm.csh list [qualifier]"
     echo "       list      is the list of mosaics to process"
     echo "       qualifier is a string append to output filenames"
      exit
  endif
endif



set listroot = "${imagelist:t:r}"
set catalogdir = $DASCH_CATALOGBIN

# Definition of color term columns
# 3 colunns 1:colorterm 2:errorcolor 3:colorflag



cd $catalogdir
pwd



set listfile = $imagelist:t
set headerflag = ""
date
echo "Processing ${DASCH_SCRIPTS}/$listfile"
echo "Writing: ${listfile:t:r}_xcolorterm.db"

if (-e ${listfile:t:r}_xcolorterm.db) then
  rm ${listfile:t:r}_xcolorterm.db
endif


foreach listentry ( `cat ${DASCH_SCRIPTS}/$listfile` ) 
  set plate = "${catalogdir}/${listentry}"

  set series1 = `echo $listentry | sed 's/_01ww//g' | sed 's/_01r90ww//g' | sed 's/_01r180ww//g' | sed 's/_01r270ww//g' | sed 's/_/ /g' | gawk '{print $1}'`
  set series = `echo $series1 | sed 's/[0-9]//g'`
  set plateNumber =  `echo $series1 | sed 's/[A-Z,a-z]//g'`
#  echo "series1 is $series1 series $series plateNumber $plateNumber"

  set colortermtxt = ${plate:t}${qualifier}.xcolorterm.txt
  if ( -e $colortermtxt ) rm $colortermtxt
  touch $colortermtxt 

  set sbin = 0
  set nbin = 0


  set catalog = `ls ${plate}${qualifier}_a?.xcolorterm.txt` 
  if ($#catalog > 0) then  

    foreach catalog ( ${plate}${qualifier}_a?.xcolorterm.txt ) 
      @ sbin = $sbin + 1 
    end
    echo "mosaic series plateNumber spatial_bin" >! spatial_bin_definition${listroot}${qualifier}.tmp
    echo "------ ------ ----------  -----------" >> spatial_bin_definition${listroot}${qualifier}.tmp
    foreach bin (1 2 3 4 5 6 7 8 9)
      set catalog = "${plate}${qualifier}_a${bin}.xcolorterm.txt"
      if (-e $catalog) then
          @ nbin = $nbin + 1
          # collect the *.para files which record the limiting magnitude in each spatial bin  
         set parameters = `cat ${catalog}`
         echo $parameters >> $colortermtxt
         echo "$listentry $series $plateNumber $bin" >> spatial_bin_definition${listroot}${qualifier}.tmp

      endif
    end
    if ($sbin != $nbin) then
        echo "ERROR number of bins $nbin does not equal the number of files $sbin"
        exit
    endif
    echo $sbin spatial bins found $listentry
    # Just add the bin number to the parameters file
  
    
    echo "colorterm errorcolor colorflag" >! para${listroot}${qualifier}.header
    echo "--------- ---------- ---------" >> para${listroot}${qualifier}.header
    cat para${listroot}${qualifier}.header $colortermtxt >! para${listentry}${qualifier}.2.tmp

  
    paste spatial_bin_definition${listroot}${qualifier}.tmp para${listentry}${qualifier}.2.tmp | gawk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7}' | column $headerflag  >>  ${listfile:t:r}_xcolorterm.db

   set headerflag = "-b"
  


    rm para${listentry}${qualifier}.2.tmp
    rm $colortermtxt
    rm spatial_bin_definition${listroot}${qualifier}.tmp
    rm para${listroot}${qualifier}.header
  else 

    echo "0 spatial bins found: $listentry" 
  endif
end
date 
