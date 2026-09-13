# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#
#  Xcol is the column of the Julian Day
#  Ycol is the column of the magnitude
#  Ecol is the column of the rms
#  if zeropoint is zero, then use the first Xcol entry.
#

set scripts = $DASCH_SCRIPTS

set lc = $1
set period = $2
set zeropoint = $3
set Xcol = $4
set Ycol = $5
set Ecol = $6

if ( $#argv != 6 ) then 
   echo "usage: lc period zeropoint format"
   echo "format = Xcol Ycol Ecol"
   echo "if no errors, Ecol=0"  
   exit
endif

if ( ! -e $lc ) then
  echo "File not found" $lc
  exit
else
  echo $lc
endif


set folded = ${lc}.p=${period}.fold

if ( $Ecol > 0 ) then
  cat $lc | gawk "{print $ $Xcol, $ $Ycol, $ $Ecol}" | sort -n -k 1,1 >! lc_sorted.tmp
else
  cat $lc | gawk "{print $ $Xcol, $ $Ycol, 0.0}" | sort -n -k 1,1 >! lc_sorted.tmp
endif

echo lc_sorted.tmp >! foldpar.tmp
echo $period >> foldpar.tmp
echo $zeropoint >> foldpar.tmp
  
${scripts}/fold < foldpar.tmp | sed 's/\*/9/g' >! $folded
cat ${scripts}/plot_fold_template.wip | sed "s/FILE/${folded}/g" | sed "s/PERIOD/${period}/g" >! ${folded}.wip
#wip ${folded}.wip

rm foldpar.tmp
rm lc_sorted.tmp
