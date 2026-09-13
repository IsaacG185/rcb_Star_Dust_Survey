#! /usr/bin/env python3
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Filter Sextractor table by limiting the number of detection
# to the N brightest objects
#

import sys, os
import numpy as np
from astropy.io import ascii, fits

import ldac_tools

if len(sys.argv) < 4:
    print("Filter FITS_LDAC catalog. Keep nkeep rows with lower MAG_ISO")
    print("Usage: sexcat_filter.py fin nkeep fout fitsfile")
    print("  fin: input filename")
    print("  nkeep: number of rows to keep")
    print("  fout: output filename")
    print("  fitsfile: fits filename header")
    sys.exit()

fin = sys.argv[1]
nkeep = int(sys.argv[2])
fout = sys.argv[3]
fitsfile = sys.argv[4]

ext = os.path.splitext(fin)[-1]

if ext == ".fits":
    # read fits_ldac table
    # tset = atpy.TableSet(fin, verbose=False)
    # tobj = tset.tables[1]
    # print("  extension is .fits ")
    raise Exception(".fits inputs no longer supported")

# read db and prepare fits_ldac table
tobj = ascii.read(fin, data_start=2)
nobj = len(tobj)

# cols: X_IMAGE Y_IMAGE ERRA_IMAGE ERRB_IMAGE FLUX_ISO FLUXERR_ISO FLAGS X_WORLD Y_WORLD MAG_ISO
tobj.rename_column("BFLAGS", "FLAGS")
tobj.rename_column("ra", "X_WORLD")
tobj.rename_column("dec", "Y_WORLD")

ndet = len(tobj)
if ndet > nkeep:
    order = np.argsort(tobj["MAG_ISO"])
    tobj = tobj[order[:nkeep]]

with fits.open(fitsfile) as hdul:
    header = hdul[0].header

hdu_out = ldac_tools.table_as_ldac(header, tobj)
hdu_out.writeto(fout, overwrite=True)
print(len(tobj), "/", ndet, "sources kept in", fout)

# 2010-06-01 Mathieu Servillat
# 2010-09-27 Edward J. Los - Revise code for .db files
