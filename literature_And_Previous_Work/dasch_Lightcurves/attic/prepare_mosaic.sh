#! /bin/bash
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#
# Check to see if the mosaic exists.  If because of previous processing,
# only the _tnx version exists, recreate the original version from the
# MySQL database using the restoretan function.
#
# if DASCH_COMPLETED is defined as a valid directory, copy the mosaic
# to that directory.
#
# Invokes:
#
# - getdirectory
# - restoretan
# - setScampRMS
#
# Variables used:
#
# - $DASCH_COMPLETED

if [ $# -ne 2 ]; then
    echo "usage: $0 <solutionNumber> <list>"
    exit 1
fi

solutionNumber=$1
imagelist=$2

counter=1
completedFlag=$(echo "$DASCH_COMPLETED" | wc -c)

if [ $solutionNumber -ne 0 ]; then
    # We are done. Multiple solutions work only with headers, not the full fits file.
    exit 0
fi

echo "Processing $imagelist Copy location $DASCH_COMPLETED"

for plate in $(cat $imagelist) ; do
    # Be sure that we are not overriding a better solution with the full mosaics
    # algorithm.
    setScampRMS -c -p $plate -e $solutionNumber
    if [ $? -ne 0 ]; then
        continue
    fi

    datestr=$(date)
    echo "$datestr $counter $plate"
    ((counter+=1))

    directory=$(getdirectory $plate)
    origfilename=$directory/${plate}.fit
    tnxfilename=$directory/${plate}_tnx.fit

    # First look for the file in the target directory and exit if it is there.
    # Next look for the tnx file in the target directory.  Convert it and exit if it is there

    if [ $completedFlag -gt 2 ]; then
        origcomplname=$DASCH_COMPLETED/${plate}.fit
        tnxcomplname=$DASCH_COMPLETED/${plate}_tnx.out

        if [ -e $origcomplname ]; then
            # File is where we expect it.  We are done
            echo "Found (1): $origcomplname"
            continue
        elif [ -e $tnxcomplname ]; then
            # Have the tnx version.  Copy and convert
            echo "Found (2): $tnxcomplname"
            cp $tnxcomplname $origcomplname
            chmod 664 $origcomplname
            restoretan $origcomplname
            continue
        fi
    fi

    # Check to see if the file is in the mosaic directory.  If only the tnx
    #  version is there, then copy and covert it

    if [ -e $origfilename ]; then
        echo "Found (3): $origfilename"
    elif [ -e $tnxfilename ]; then
        echo "Found (4): $tnxfilename"
        chmod 664 $tnxfilename
        mv $tnxfilename $origfilename
        restoretan $origfilename
    else
        echo "Not found: $tnxfilename"
        continue
    fi

    if [ $completedFlag -gt 2 ]; then
        # File is in the original directory.  Copy it to the completed directory
        cp $origfilename $origcomplname
    fi
done

#  Nov  7, 2007 Edward J. Los - First version to replace copylist.perl
#  Nov 22, 2008 Edward J. Los - Change "cp" to "mv" to save disk space
#  Jun  9, 2009 Edward J. Los - change .fits to .out because of changes in refine_wcs.csh
#  Sep 29, 2009 Edward J. Los - add multiple exposure support
#                               NOTE: this routine is not used for multiple exposures; only
#                               for the first solution.
#  Apr  8, 2016 Edward J. Los - Handle read-only mosiac protection
#  Sep 12, 2016 Edward J. Los - Back out read-only change
#  Nov  1, 2016 Edward J. Los - Make sure the mosaic is writable before executin restoretan
