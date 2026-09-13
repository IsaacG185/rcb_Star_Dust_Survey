# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# 
# script wrapper for contour_plot.f
#  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead

set scripts = $DASCH_SCRIPTS

if ( $#argv < 8 ) then
	echo "usage: file nx ny xtext ytext ncon plotdev nclip scale outfile"
	exit
endif

set infile = $1
set nx = $2
set ny = $3
set xtext = $4
set ytext = $5
set ncontours  = $6
set plotdev = $7
set nclip = $8
set scale = $9
set outfile = $10

if ( -e par${infile} ) rm par${infile}
if ( -e $outfile ) rm $outfile

echo $infile > par${infile}
echo $nx $ny >> par${infile}
echo $xtext >> par${infile}
echo $ytext >> par${infile}
echo $ncontours >> par${infile}
echo $plotdev >> par${infile}
echo $nclip >> par${infile}
echo $scale >> par${infile}
echo $outfile >> par${infile}
echo " " >> par${infile}

echo "${scripts}/contour_plot_clip_smooth_correct2 < par${infile}"
${scripts}/contour_plot_clip_smooth_correct2 < par${infile}

#rm par${infile}

