# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#  This script searches a catalog for all stars within a specified tolerance
#  of a given location.
#
#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#  Nov  9, 2007 Edward J. Los - make RA and DEC parameters.
#
set scripts = $DASCH_SCRIPTS
set data = $DASCH_COMPLETED
set catalog = $DASCH_CATALOG
if ($#argv != 3) then
    echo "usage: get_nearest  ra dec tolerance"
    echo "       ra and dec are in degrees"
    echo "       tolerance is in arcsec"
    echo "       catalog is defined by DASCH_CATALOG"
    exit

endif


set ra = $1
set dec = $2 
set radius = $3
# near 3c273
# set ra = 187.28541
# set dec =  2.00572
echo "tolerance $radius arcsec"
echo "ra: $ra"
echo "dec: $dec " 
echo "catalog $catalog"

set conditions = "(ra > ($ra-($radius/3600.))) && (ra < ($ra+($radius/3600.))) && (dec > ($dec-($radius/3600.))) && (dec < ($dec+($radius/3600.)))"

echo "conditions: $conditions"

#echo "row -i $catalog $conditions >! result.tmp"
row -i $catalog $conditions >! result.tmp



column -i result.tmp -a dra ddec drad | compute  "ddec=3600*(dec-$dec); dra=3600*((ra-$ra)*(cos(((dec+$dec)/2)/57.29577951))); drad=sqrt(dra^2+ddec^2)" 
