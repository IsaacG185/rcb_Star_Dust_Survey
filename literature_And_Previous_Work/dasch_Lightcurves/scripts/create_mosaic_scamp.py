#! /usr/bin/env python3
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# create new mosaic and update header info

import sys
import shutil
from unicodedata import normalize

from astropy.io import fits

if len(sys.argv) < 5:
    print("Create _scamp.fit mosaic from initial mosaic and scamp header keywords")
    print("Usage: create_mosaic_scamp.py img_name hd_scamp_name")
    print("  img_name: initial mosaic fits image")
    print("  hd_scamp_name: .head output from scamp")
    print("  img_out: out mosaic file")
    print("  solutionNumber: solution number")
    sys.exit()

img_name = sys.argv[1]
hd_scamp_name = sys.argv[2]
img_out = sys.argv[3]
solutionNumber = sys.argv[4]

hdr_out = hd_scamp_name.split(".sexcat")[0] + "_scamp.hdr"
verbose = 0

# read scamp solution in .head file


def to_ascii(text: str) -> str:
    return normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")


with open(hd_scamp_name, "rt") as fin:
    lines = [to_ascii(l) for l in fin]

end = "END     \n"
if end in lines:
    lines.remove(end)

with open(hd_scamp_name + "_cor", "wt") as fout:
    fout.writelines(lines)

# now read
hd_scamp = fits.Header.fromtextfile(hd_scamp_name + "_cor")

if solutionNumber == "0":
    shutil.copy(img_name, img_out)
    hdul = fits.open(img_out, mode="update")
    hd = hdul[0].header
else:
    print("ERROR: pyfits does not handle headers correctly for " + img_name)
    exit(1)

hd["EQUINOX"] = (2000.0, "Epoch of RA & DEC")

# del previous PV keywords and all astrometry related keywords
if solutionNumber == "0":
    for kw in list(hd.keys()):
        for kwtype in "PV WAT LTM LTV WCSDIM CTYPE CRVAL CRPIX CD".split():
            if kwtype in kw:
                del hd[kw]
else:
    for axisIndex in (1, 2):
        for pvIndex in range(100):
            del hd["PV" + str(axisIndex) + "_" + str(pvIndex)]

for kw in list(hd_scamp.keys()):
    if kw != "HISTORY" and kw != "COMMENT":
        comm = ""
        if " / " in str(hd_scamp[kw]):
            comm = str(hd_scamp[kw]).split(" / ", 1)[-1]
        if verbose == 1:
            if kw in list(hd.keys()):
                print(kw, "=", hd[kw], "-->", hd_scamp[kw])
            else:
                print(kw, "=", hd_scamp[kw])
        hd[kw] = (hd_scamp[kw], comm)

if solutionNumber == "0":
    hdul.close()
else:
    hdu = fits.PrimaryHDU()
    hdu.header = hd
    hdu.writeto(img_out)

# 2009-12-10 Mathieu Servillat
# 2010-09-27 Edward J. Los - do not save the header because it is not compatible with other saved headers.
#                          - for solution numbers > 0, delete PV keywords explicitly and do not delete previous keywords
