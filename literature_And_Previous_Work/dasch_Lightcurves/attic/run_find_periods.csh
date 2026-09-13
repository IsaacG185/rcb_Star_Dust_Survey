# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#  Mar 24, 2008 Edward J. Los - use only the good points from lightcurves.

set scripts = $DASCH_SCRIPTS

if ( $#argv != 5 ) then
	echo "usage: catalog method p1 p2 intF lightcurvedir"
        echo "       method is either scargle or pdm"
        echo "       p1 = starting period in days"
        echo "       p2 = ending period indays"
        echo "       intF = test interval in days"
	exit
endif
set scargle = `echo "$2" | gawk '/scargle/{print $1}'`
set scargle = $#scargle
set pdm = `echo "$2" | gawk '/pdm/{print $1}'`
set pdm = $#pdm
if (($pdm != 1) && ($scargle != 1) ) then
    echo "method must be either scargle and/or pdm"
    exit
endif


set p1 = $3
set p2 = $4
set intF = $5
set lightcurvedir = $DASCH_LIGHTCURVES
set scripts = $DASCH_SCRIPTS

echo "scargle: $scargle, pdm: $pdm, p1: $p1, p2: $p2, intF $intF, lightcurvedir $lightcurvedir";

foreach set ( $1 )
  date
  cd $lightcurvedir
  if ( ! -e $set ) mkdir $set 
  cd ${lightcurvedir}/${set}
  pwd
  foreach lc ( lc*_good.db )
    echo "${scripts}/find_periods.csh ${lc:r}.col $p1 $p2 $intF $scargle $pdm"
    ${scripts}/find_periods.csh ${lc:r}.col $p1 $p2 $intF $scargle $pdm
  end



end
date





