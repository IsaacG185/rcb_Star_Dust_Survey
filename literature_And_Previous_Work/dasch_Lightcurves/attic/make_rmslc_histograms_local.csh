# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
#  Mar 24, 2008 Edward J. Los - Use the summary files in the ingest directory
#
set platelist = $1
set dbroot = $platelist:t:r
set db = ${DASCH_INGEST}/id_all_${dbroot}.db
set set = ${DASCH_INGEST}/${dbroot}

echo "make_rmslc_histograms_local.csh with plate database $db"
date

if ( ! -e $set ) mkdir $set 
cd ${set}
pwd


set conditions = 'ngood>10&&rms_local<1'

row -i $db $conditions | histtable -min 0 -max 1 -n 100 rms_local | column -b center count >! rmslc_all_local.hist
row -i $db $conditions | row 'median_local>=7&&median_local<8' | histtable -min 0 -max 1 -n 100 rms_local | column -b center count >! rmslc_b_7-8_local.hist
row -i $db $conditions | row 'median_local>=8&&median_local<9' | histtable -min 0 -max 1 -n 100 rms_local | column -b center count >! rmslc_b_8-9_local.hist
row -i $db $conditions | row 'median_local>=9&&median_local<10' | histtable -min 0 -max 1 -n 100 rms_local | column -b center count >! rmslc_b_9-10_local.hist
row -i $db $conditions | row 'median_local>=10&&median_local<11' | histtable -min 0 -max 1 -n 100 rms_local | column -b center count >! rmslc_b_10-11_local.hist
row -i $db $conditions | row 'median_local>=11&&median_local<12' | histtable -min 0 -max 1 -n 100 rms_local | column -b center count >! rmslc_b_11-12_local.hist
row -i $db $conditions | row 'median_local>=12&&median_local<13' | histtable -min 0 -max 1 -n 100 rms_local | column -b center count >! rmslc_b_12-13_local.hist
row -i $db $conditions | row 'median_local>=13&&median_local<14' | histtable -min 0 -max 1 -n 100 rms_local | column -b center count >! rmslc_b_13-14_local.hist
row -i $db $conditions | row 'median_local>=14&&median_local<15' | histtable -min 0 -max 1 -n 100 rms_local | column -b center count >! rmslc_b_14-15_local.hist

row -i $db $conditions | sorttable -nr rms_local | column -b rms_local | gawk '{print $1,NR}' >! rmslc_all_local.cumulative
row -i $db $conditions | row 'median_local>=7&&median_local<8' | sorttable -nr rms_local | column -b rms_local | gawk '{print $1,NR}' >! rmslc_b_7-8_local.cumulative
row -i $db $conditions | row 'median_local>=8&&median_local<9' |  sorttable -nr rms_local | column -b rms_local | gawk '{print $1,NR}' >! rmslc_b_8-9_local.cumulative
row -i $db $conditions | row 'median_local>=9&&median_local<10' | sorttable -nr rms_local | column -b rms_local | gawk '{print $1,NR}' >! rmslc_b_9-10_local.cumulative
row -i $db $conditions | row 'median_local>=10&&median_local<11' | sorttable -nr rms_local | column -b rms_local | gawk '{print $1,NR}' >! rmslc_b_10-11_local.cumulative
row -i $db $conditions | row 'median_local>=11&&median_local<12' | sorttable -nr rms_local | column -b rms_local | gawk '{print $1,NR}' >! rmslc_b_11-12_local.cumulative
row -i $db $conditions | row 'median_local>=12&&median_local<13' | sorttable -nr rms_local | column -b rms_local | gawk '{print $1,NR}' >! rmslc_b_12-13_local.cumulative
row -i $db $conditions | row 'median_local>=13&&median_local<14' | sorttable -nr rms_local | column -b rms_local | gawk '{print $1,NR}' >! rmslc_b_13-14_local.cumulative
row -i $db $conditions | row 'median_local>=14&&median_local<15' | sorttable -nr rms_local | column -b rms_local | gawk '{print $1,NR}' >! rmslc_b_14-15_local.cumulative



# Compute the Median RMS
set a0 = `row -i $db $conditions | column rms_local | mediantable | column -b`
if ($#a0 == 0) set a0 = 0.0
set a1 = `row -i $db $conditions | row 'median_local>=7&&median_local<8' | column rms_local | mediantable | column -b`
if ($#a1 == 0) set a1 = 0.0
set a2 = `row -i $db $conditions | row 'median_local>=8&&median_local<9' | column rms_local | mediantable | column -b`
if ($#a2 == 0) set a2 = 0.0
set a3 = `row -i $db $conditions | row 'median_local>=9&&median_local<10' | column rms_local | mediantable | column -b`
if ($#a3 == 0) set a3 = 0.0
set a4 = `row -i $db $conditions | row 'median_local>=10&&median_local<11' | column rms_local | mediantable | column -b`
if ($#a4 == 0) set a4 = 0.0
set a5 = `row -i $db $conditions | row 'median_local>=11&&median_local<12' | column rms_local | mediantable | column -b`
if ($#a5 == 0) set a5 = 0.0
set a6 = `row -i $db $conditions | row 'median_local>=12&&median_local<13' | column rms_local | mediantable | column -b`
if ($#a6 == 0) set a6 = 0.0
set a7 = `row -i $db $conditions | row 'median_local>=13&&median_local<14' | column rms_local | mediantable | column -b`
if ($#a7 == 0) set a7 = 0.0
set a8 = `row -i $db $conditions | row 'median_local>=14&&median_local<15' | column rms_local | mediantable | column -b`
if ($#a8 == 0) set a8 = 0.0

# Compute the 90th Percentile RMS (least variable)
set n = `cat rmslc_all_local.cumulative | wc -l | gawk '{print int($1*0.9)}'` 

if ($n > 0) then
  set b0 = `sed -n $n"p" rmslc_all_local.cumulative | gawk '{print $1}'`
else
  set b0 = 0.0
endif

set n = `cat rmslc_b_7-8_local.cumulative | wc -l | gawk '{print int($1*0.9)}'`
if ($n > 0) then
  set b1 = `sed -n $n"p" rmslc_b_7-8_local.cumulative | gawk '{print $1}'`
else
  set b1 = 0.0
endif
set n = `cat rmslc_b_8-9_local.cumulative | wc -l | gawk '{print int($1*0.9)}'`
if ($n > 0) then
  set b2 = `sed -n $n"p" rmslc_b_8-9_local.cumulative | gawk '{print $1}'`
else
  set b2 = 0.0
endif
set n = `cat rmslc_b_9-10_local.cumulative | wc -l | gawk '{print int($1*0.9)}'`
if ($n > 0) then
  set b3 = `sed -n $n"p" rmslc_b_9-10_local.cumulative | gawk '{print $1}'`
else
  set b3 = 0.0
endif
set n = `cat rmslc_b_10-11_local.cumulative | wc -l | gawk '{print int($1*0.9)}'`
if ($n > 0) then
  set b4 = `sed -n $n"p" rmslc_b_10-11_local.cumulative | gawk '{print $1}'`
else
  set b4 = 0.0
endif
set n = `cat rmslc_b_11-12_local.cumulative | wc -l | gawk '{print int($1*0.9)}'`
if ($n > 0) then
  set b5 = `sed -n $n"p" rmslc_b_11-12_local.cumulative | gawk '{print $1}'`
else
  set b5 = 0.0
endif
set n = `cat rmslc_b_12-13_local.cumulative | wc -l | gawk '{print int($1*0.9)}'`
if ($n > 0) then
  set b6 = `sed -n $n"p" rmslc_b_12-13_local.cumulative | gawk '{print $1}'`
else
  set b6 = 0.0
endif
set n = `cat rmslc_b_13-14_local.cumulative | wc -l | gawk '{print int($1*0.9)}'`
if ($n > 0) then
  set b7 = `sed -n $n"p" rmslc_b_13-14_local.cumulative | gawk '{print $1}'`
else
  set b7 = 0.0
endif
set n = `cat rmslc_b_14-15_local.cumulative | wc -l | gawk '{print int($1*0.9)}'`
if ($n > 0) then
  set b8 = `sed -n $n"p" rmslc_b_14-15_local.cumulative | gawk '{print $1}'`
else
  set b8 = 0.0
endif

# Compute the 10th Percentile RMS (most variable)
set n = `cat rmslc_all_local.cumulative | wc -l | gawk '{print int($1*0.1)}'`
if ($n > 0) then
  set c0 = `sed -n $n"p" rmslc_all_local.cumulative | gawk '{print $1}'`
else
  set c0 = 0.0
endif
set n = `cat rmslc_b_7-8_local.cumulative | wc -l | gawk '{print int($1*0.1)}'`
if ($n > 0) then
  set c1 = `sed -n $n"p" rmslc_b_7-8_local.cumulative | gawk '{print $1}'`
else
  set c1 = 0.0
endif
set n = `cat rmslc_b_8-9_local.cumulative | wc -l | gawk '{print int($1*0.1)}'`
if ($n > 0) then
  set c2 = `sed -n $n"p" rmslc_b_8-9_local.cumulative | gawk '{print $1}'`
else
  set c2 = 0.0
endif
set n = `cat rmslc_b_9-10_local.cumulative | wc -l | gawk '{print int($1*0.1)}'`
if ($n > 0) then
  set c3 = `sed -n $n"p" rmslc_b_9-10_local.cumulative | gawk '{print $1}'`
else
  set c3 = 0.0
endif
set n = `cat rmslc_b_10-11_local.cumulative | wc -l | gawk '{print int($1*0.1)}'`
if ($n > 0) then
  set c4 = `sed -n $n"p" rmslc_b_10-11_local.cumulative | gawk '{print $1}'`
else
  set c4 = 0.0
endif
set n = `cat rmslc_b_11-12_local.cumulative | wc -l | gawk '{print int($1*0.1)}'`
if ($n > 0) then
  set c5 = `sed -n $n"p" rmslc_b_11-12_local.cumulative | gawk '{print $1}'`
else
  set c5 = 0.0
endif
set n = `cat rmslc_b_12-13_local.cumulative | wc -l | gawk '{print int($1*0.1)}'`
if ($n > 0) then
  set c6 = `sed -n $n"p" rmslc_b_12-13_local.cumulative | gawk '{print $1}'`
else
  set c6 = 0.0
endif
set n = `cat rmslc_b_13-14_local.cumulative | wc -l | gawk '{print int($1*0.1)}'`
if ($n > 0) then
  set c7 = `sed -n $n"p" rmslc_b_13-14_local.cumulative | gawk '{print $1}'`
else
  set c7 = 0.0
endif
set n = `cat rmslc_b_14-15_local.cumulative | wc -l | gawk '{print int($1*0.1)}'`
if ($n > 0) then
  set c8 = `sed -n $n"p" rmslc_b_14-15_local.cumulative | gawk '{print $1}'`
else
  set c8 = 0.0
endif

#echo "a1: $a1 a2: $a2 a3: $a3 a4: $a4 a5: $a5 a6: $a6 a7: $a7 a8: $a8" 
#echo "b1: $b1 b2: $b2 b3: $b3 b4: $b4 b5: $b5 b6: $b6 b7: $b7 b8: $b8" 
#echo "c1: $c1 c2: $c2 c3: $c3 c4: $c4 c5: $c5 c6: $c6 c7: $c7 c8: $c8" 

# fix to avoid zeros when plotting on log-scale.

foreach hist ( rmslc_b_*-*_local.hist )
  cat $hist | gawk '$2==0{print $1,"0.1"};$2>0{print $0}' >! hist_local.tmp
  mv hist_local.tmp $hist
end

# Write the median and 10th, 90th percentiles for plotting vs magnitude
echo 7.5  $a1 1 $b1 $c1 >! rms_vs_mag_local.col
echo 8.5 $a2 2 $b2 $c2 >> rms_vs_mag_local.col
echo 9.5  $a3 3 $b3 $c3 >> rms_vs_mag_local.col
echo 10.5 $a4 4 $b4 $c4 >> rms_vs_mag_local.col
echo 11.5 $a5 5 $b5 $c5 >> rms_vs_mag_local.col
echo 12.5 $a6 6 $b6 $c6 >> rms_vs_mag_local.col
echo 13.5 $a7 7 $b7 $c7 >> rms_vs_mag_local.col
echo 14.5 $a8 8 $b8 $c8 >> rms_vs_mag_local.col

# Write the global median and 10th, 90th percentiles for plotting over the cumulative plot.
echo $a0 1e-6 $b0 $c0 >! rms_quantiles_local.col
echo $a0 1e6 $b0 $c0 >> rms_quantiles_local.col

# Clean up null plots
set n = `wc -l rmslc_b_7-8_local.cumulative | gawk '{print $1}'`

if ($n == 0) then 
    echo "0.01 1" >! rmslc_b_7-8_local.cumulative
endif
set n = `wc -l rmslc_b_8-9_local.cumulative | gawk '{print $1}'`

if ($n == 0) then 
  echo "0.01 1" >! rmslc_b_8-9_local.cumulative
endif
set n = `wc -l rmslc_b_9-10_local.cumulative | gawk '{print $1}'`

if ($n == 0) then 
  echo "0.01 1" >! rmslc_b_9-10_local.cumulative
endif
set n = `wc -l rmslc_b_10-11_local.cumulative | gawk '{print $1}'`

if ($n == 0) then
  echo "0.01 1" >! rmslc_b_10-11_local.cumulative
endif
set n = `wc -l rmslc_b_11-12_local.cumulative | gawk '{print $1}'`

if ($n == 0) then
  echo "0.01 1" >! rmslc_b_11-12_local.cumulative
endif
set n = `wc -l rmslc_b_12-13_local.cumulative | gawk '{print $1}'`

if ($n == 0) then
  echo "0.01 1" >! rmslc_b_12-13_local.cumulative
endif
set n = `wc -l rmslc_b_13-14_local.cumulative | gawk '{print $1}'`

if ($n == 0) then
  echo "0.01 1" >! rmslc_b_13-14_local.cumulative
endif
set n = `wc -l rmslc_b_14-15_local.cumulative | gawk '{print $1}'`

if ($n == 0) then
  echo "0.01 1" >! rmslc_b_14-15_local.cumulative
endif

echo 

#wip $DASCH_SCRIPTS/plot_rmslc_logscale_local.wip
#wip $DASCH_SCRIPTS/plot_rmslc_logscale_cumulative_local.wip

cd $DASCH_SCRIPTS
pwd

