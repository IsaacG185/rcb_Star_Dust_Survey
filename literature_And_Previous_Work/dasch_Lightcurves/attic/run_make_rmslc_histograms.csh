# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#  Mar 24, 2008 Edward J. Los - Use the summary files in the ingest directory

# run lightcurve script on a bunch of source lists
set platelist = $1
set scripts = $DASCH_SCRIPTS

echo "platelist is $1"

  date
  pwd
#  ${scripts}/extract_lightcurves_by_id5_local.csh ${scripts}/${imagelist} ${scripts}/${set}_REF.db $ingestdir ${scripts}/JulianDates.txt
#  ${scripts}/lc_scatter_flags_local.csh
  ${scripts}/make_rmslc_histograms_local.csh $platelist

date


