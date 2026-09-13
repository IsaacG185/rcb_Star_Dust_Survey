# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
# run_fold.csh
#    Parameter 1: fold list consisting of database, star, and period.
#    Parameter 2: 'noflags' to show all points the same way
#                 'flags' to show points identified by detection flags
#    Parameter 3: remove to remove bad points
#                 keep to keep bad points
#
#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#  Mar 24, 2008 Edward J. Los - Plot only local calibration
#

set scripts = $DASCH_SCRIPTS
echo "argv is $#argv"
if ( $#argv != 3 ) then
	echo "usage: list symbols bad"
        echo "       the list contains the subdirectory, star designation and period"
        echo "       symbols is 'noflags' for circles, 'flags' to show bad points"
        echo "       bad is 'remove' to remove bad stars, 'keep' to keep bad points"
        echo "       NOTE: noflags and remove is not a valid combination."
	exit
endif


set starlist = $1
set lightcurvedir = $DASCH_LIGHTCURVES


set flags = `echo "$2" | gawk '/noflags/{print $1}'`
set flags = $#flags
@ flags = 1 - $flags

set keep = `echo "$3" | gawk '/keep/{print $1}'`
set keep = $#keep

echo "Processing $starlist, flags: $flags keep: $keep"
set iterator = 0

foreach listentry (`cat $starlist`)
    if ($iterator == 0) then
        set subdirectory = $listentry
        set iterator = 1
        continue
    endif
    if ($iterator == 1) then
        set star = $listentry
        set iterator = 2
        continue
    endif
    set period = $listentry
    set iterator = 0;
    

    echo "subdirectory $subdirectory star: $star, period: $period"

    date
    cd ${lightcurvedir}/${subdirectory}
    pwd
    echo "star: $star, period: $period"

    if ($flags == 1) then
        ${scripts}/foldflags.csh lc_${star}_all.col $period 0 2 5 6 4

       if ($keep == 0)  then
          cat lc_${star}_all.col.p=${period}.fold | gawk '(($3>0) && ($5<30)) {print $0}' >! good.tmp
          mv good.tmp lc_${star}_all.col.p=${period}.fold

       endif
    else
        ${scripts}/fold.csh lc_${star}_all.col $period 0 2 5 6

       if ($keep == 0)  then
          echo "ERROR: noflags and noremove is not a valid combination"
       endif

     endif



end
date
