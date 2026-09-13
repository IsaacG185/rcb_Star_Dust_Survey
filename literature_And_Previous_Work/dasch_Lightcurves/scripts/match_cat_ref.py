#! /usr/bin/env python3
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Match catalog with reference catalog and save match results

import sys

from astropy.table import Table
import numpy as np

import sourcelists as sl

if len(sys.argv) < 8:
    print(
        "Usage: match_cat.py matchdirectory imagesdir cat_sex_name cat_ref_name cat_out_name radius"
    )
    print("  matchdir:     directory containing catalog files")
    print("  imagesdir:    directory for region files")
    print("  cat_sex_name: Sextractor FITS_LDAC output catalog of detections")
    print("  cat_ref_name: Reference catalog used by scamp in FITS_LDAC format")
    print(
        "  cat_out_name: Output fits file with RA, Dec, RA_REF, Dec_REF, dra, ddec, drad and offset"
    )
    print("  radius: match radius used in cross correlation")
    print("  scale:  plate scale in arcsec per pixel")
    sys.exit()

matchdir = sys.argv[1]
imagesdir = sys.argv[2]
cat_sex_name = sys.argv[3]
cat_ref_name = sys.argv[4]
cat_out_name = sys.argv[5]
radius = float(sys.argv[6])
scale = float(sys.argv[7])

cat_sex = Table.read(matchdir + "/" + cat_sex_name, hdu=2)
print("Sextractor catalog opened:" + matchdir + "/" + cat_sex_name)

ra1 = cat_sex["X_WORLD"]
dec1 = cat_sex["Y_WORLD"]
sl.save_ds9reg(
    imagesdir + "/" + cat_sex_name + ".xy.reg",
    cat_sex["X_IMAGE"],
    cat_sex["Y_IMAGE"],
    errors=4.0,
    header="physical",
    tag="sexcat_xy",
)
sl.save_ds9reg(
    imagesdir + "/" + cat_sex_name + ".reg", ra1, dec1, errors=4.0, tag="sexcat"
)

cat_ref = Table.read(matchdir + "/" + cat_ref_name, hdu=2)
print("Reference catalog opened:" + matchdir + "/" + cat_ref_name)

ra2 = cat_ref["X_WORLD"]
dec2 = cat_ref["Y_WORLD"]
sl.save_ds9reg(
    imagesdir + "/" + cat_ref_name + ".reg",
    ra2,
    dec2,
    errors=4.0,
    color="red",
    tag="ucac3",
)

ind1, ind2, dist12 = sl.match2(
    ra1, dec1, ra2, dec2, errors=radius, radius=radius, onlybest=True
)

if len(ind1) > 0:
    dra = (ra1[ind1] - ra2[ind2]) * np.cos(np.deg2rad(dec2[ind2])) * 3600.0
    ddec = (dec1[ind1] - dec2[ind2]) * 3600.0
    ra_center = np.mean(ra1)
    dec_center = np.mean(dec1)
    offset = sl.separation(ra1[ind1], dec1[ind1], ra_center, dec_center)

    cat_out = Table(
        {
            "RA": ra1[ind1],
            "Dec": dec1[ind1],
            "RA_REF": ra2[ind2],
            "Dec_REF": dec2[ind2],
            "dra": dra,
            "ddec": ddec,
            "drad": dist12,
            "offset": offset,
        }
    )

    cat_out.table_name = matchdir + "/" + cat_out_name
    cat_out.write(matchdir + "/" + cat_out_name, overwrite=True)

    # compute zero-based standard deviation
    rastd = np.sqrt(np.dot(dra, dra) / len(dra))
    clipdra = []
    for j in range(len(dra)):
        if dra[j] < 3 * rastd:
            clipdra = np.append(clipdra, dra[j])
    cliprastd = np.sqrt(np.dot(clipdra, clipdra) / len(clipdra))
    print(
        "len(dra) ",
        len(dra),
        " rastd ",
        rastd / scale,
        "len(clipdra)",
        len(clipdra),
        " cliprastd ",
        cliprastd / scale,
        " for ",
        cat_sex_name,
    )

    decstd = np.sqrt(np.dot(ddec, ddec) / len(ddec))
    clipddec = []
    for j in range(len(ddec)):
        if ddec[j] < 3 * decstd:
            clipddec = np.append(clipddec, ddec[j])
    clipdecstd = np.sqrt(np.dot(clipddec, clipddec) / len(clipddec))
    print(
        "len(ddec) ",
        len(ddec),
        " decstd ",
        decstd / scale,
        "len(clipddec)",
        len(clipddec),
        " clipdecstd ",
        clipdecstd / scale,
        " for ",
        cat_sex_name,
    )

    minlen = min([len(dra), len(ddec)])
    fout = open(imagesdir + "/" + cat_sex_name + ".rms", "w")
    fout.write(
        "length "
        + str(minlen)
        + " RMS RA "
        + str(cliprastd / scale)
        + " RMS DEC "
        + str(clipdecstd / scale)
        + " for "
        + cat_sex_name
        + "\n"
    )
    fout.close()

    print("Match catalog saved", matchdir + "/" + cat_out_name)
else:
    print("ERROR: No match catalog saved", matchdir + "/" + cat_out_name)

# 2009-12-10 Mathieu Servillat
# 2010-09-06 Edward J. Los
# 2015-10-24 Edward J. los Add     cat_out.table_name = matchdir+'/'+cat_out_name
