#! /usr/bin/env python3
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Get num of detection in Sextractor table

import sys

from astropy.io import fits

if len(sys.argv) < 2:
    print("Return number of rows in FITS_LDAC catalog")
    print("Usage: sexcat_ndet.py fin")
    print("  fin: input filename")
    sys.exit()

fin = sys.argv[1]
data = fits.getdata(fin, 2)
print(len(data))

# 2010-06-01 Mathieu Servillat
