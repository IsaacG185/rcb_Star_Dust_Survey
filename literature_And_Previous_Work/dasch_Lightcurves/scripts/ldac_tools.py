# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# This was based on:
#
# https://github.com/winter-telescope/mirar/blob/741571cc5612df565ae2cbec49dead833641b4ce/mirar/utils/ldac_tools.py
#
# However the code has been extensively rewritten because it didn't actually
# work.

"""
Functions to convert FITS files or astropy Tables to FITS_LDAC files and
vice versa.
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table


def make_ldac_imhead_hdu(
    header: fits.Header,
) -> fits.BinTableHDU:
    """
    Convert a FITS header into a FITS-LDAC IMHEAD HDU.

    Parameters
    ----------
    header: `astropy.io.fits.Header`
        Header to convert.

    Returns
    -------
    hdu: `astropy.io.fits.BinTableHDU`
        The header as a FITS-LDAC LDAC_IMHEAD HDU.
    """

    # Some Astromatic software searches this header for the text "END      ". It
    # seems that somehow Astropy writes things out such that trailing whitespace
    # ends up disappearing, and the header that Astromatic sees ends only with
    # "END". Hack around this by adding a fake extra COMMENT header, which seems
    # to do the trick.
    header_text = header.tostring("") + "COMMENT" + 73 * " "
    n_headers = len(header) + 2  # make sure to get the END and fake COMMENT
    tblhdr = np.array([header_text], dtype=f"|S{len(header_text)}")
    col1 = fits.Column(
        name="Field Header Card", array=tblhdr, format=f"{len(header_text)}A"
    )
    cols = fits.ColDefs([col1])
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.header["TDIM1"] = f"(80, {n_headers})"
    hdu.header["EXTNAME"] = "LDAC_IMHEAD"
    return hdu


def table_as_ldac(
    imhead: fits.Header,
    table: Table,
) -> fits.HDUList:
    """
    Convert an Astropy table into a FITS-LDAC dataset.

    Parameters
    ----------
    imhead: `astropy.io.fits.Header`
        The associated image header to attach to the data.
    table: `astropy.table.Table`
        The data table to export.

    Returns
    -------
    hdul: `astropy.io.fits.HDUList`
        An HDU list in FITS-LDAC format
    """
    hdul = fits.HDUList()
    hdul.append(fits.PrimaryHDU())
    hdul.append(make_ldac_imhead_hdu(imhead))
    hdu = fits.table_to_hdu(table)
    hdu.header["EXTNAME"] = "LDAC_OBJECTS"
    hdul.append(hdu)
    return hdul
