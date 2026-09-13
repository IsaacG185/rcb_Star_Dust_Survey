#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#  compare_REF.csh - compare two files based on Sextractor REF
#
#  Sep 22, 2009 Edward J. Los - initial version
#

if ($#argv != 2) then
    echo "usage: compare_REF.csh <file1> <file2>"
    exit
endif
date
set file1 = $1
set file2 = $2
echo "Comparing $file1" 
echo "  and $file2"
set size1 = `wc -l $file1 | gawk '{print $1}'`
set size2 = `wc -l $file2 | gawk '{print $1}'`
#echo "size1 $size1 size2 $size2"
if ($size1 != $size2) then
    echo "ERROR: size1 $size1 and size2 $size2 do not agree"
endif
set header1 = `head -n 1 $file1`
set header2 = `head -n 1 $file2`
#echo "header1 of size $#header1 is $header1"
#echo "header2 of size $#header2 is $header2"
if ($#header1 != $#header2) then
   echo "ERROR: header1 has $#header1 columns and header2 has $#header2 columns"
endif
set counter = 1;

while ($counter <= $#header1)
if ($header1[$counter] != $header2[$counter]) then
   echo "ERROR: headers are not the same $counter $header1[$counter] $header2[$counter]"
   @ counter = $counter + 1
endif
set joinfile = ${file1}.join.tmp
sorttable -n REF < ${file1} > ${file1}.tmp
sorttable -n REF < ${file2} > ${file2}.tmp
echo "jointable -j REF ${file1}.tmp ${file2}.tmp > $joinfile"
jointable -j REF ${file1}.tmp ${file2}.tmp > $joinfile
echo "join file is $joinfile"
set joinsize = `wc -l $joinfile | gawk '{print $1}'`
if ($size1 != $joinsize) then
   echo "ERROR: size1 $size1 and joinsize $joinsize do not agree"
endif
date
foreach column ($header1)
    if ($column == "REF") continue
    echo "column is $column"
    column -a -i $joinfile result value_1 value_2 | compute "value_1=${column}_1;value_2=${column}_2;result=value_1-value_2"  | column value_1 value_2 result > ${joinfile}.${column}.tmp
    statstable < ${joinfile}.${column}.tmp
end
date
