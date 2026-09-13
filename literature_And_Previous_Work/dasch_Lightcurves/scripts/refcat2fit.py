#!/usr/bin/env python3
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Converts ucac3 table to FITS_LDAC table

import sys
import numpy as np

from astropy.io import ascii, fits
from astropy import units as u

import ldac_tools

if len(sys.argv) < 5:
    print("Convert table from  imcat -dh -c ucac3  to FITS_LDAC table")
    print("Usage: refcat2fit.py fin fout")
    print("  fin: input filename")
    print("  fout: output filename")
    print("  nkeep: n brightest stars to keep")
    print("  mag: which UCAC magnitude for selection (MagB)")
    sys.exit()

fin = sys.argv[1]
fout = sys.argv[2]
nkeep = int(sys.argv[3])
mag = sys.argv[4]

# read table from imcat -dh -c ucac3
data = ascii.read(fin, delimiter=" ", header_start=4, data_start=5)
nref_tot = len(data)
# print nref_tot,'total reference sources'

# filter for objects with 2mass counterpart only
sel = np.where((data["MagJ"] > 6.0) * (data["MagJ"] < 17.0))
data = data[sel]
# filter for objects in magnitude
sel = np.where((data[mag] > 8.0) * (data[mag] < 16.0))
data = data[sel]
nrows = len(data)
# print nrows,'good reference sources'

# keep n brightest stars
if nrows > nkeep:
    order = np.argsort(data[mag])
    data = data[order]
    data = data[list(range(nkeep))]
nref = len(data)

data.rename_column("RA2000", "X_WORLD")
data["X_WORLD"].unit = u.deg

data.rename_column("Dec2000", "Y_WORLD")
data["Y_WORLD"].unit = u.deg

data["ERRA_WORLD"] = (np.zeros(nref, dtype=np.float32) + 0.1) / 3600.0
data["ERRA_WORLD"].unit = u.deg

data["ERRB_WORLD"] = (np.zeros(nref, dtype=np.float32) + 0.1) / 3600.0
data["ERRB_WORLD"].unit = u.deg

data.rename_column(mag, "MAG")
data["MAG"].unit = u.mag

data["MAGERR"] = np.zeros(nref, dtype=np.float32) + 0.1
data["MAGERR"].unit = u.mag

data["FLAGS"] = np.zeros(nref, dtype=np.int32)

# save new file, ready for scamp
hdu_out = ldac_tools.table_as_ldac(fits.PrimaryHDU().header, data)
hdu_out.writeto(fout, overwrite=True)

print(
    len(data),
    "reference sources kept in",
    fout,
    mag,
    np.min(data["MAG"]),
    "to",
    np.max(data["MAG"]),
)

# 2010-06-01 Mathieu Servillat
# 2010-09-27 Edward J. Los Add a zero FLAGS column to avoid SCAMP warnings.
