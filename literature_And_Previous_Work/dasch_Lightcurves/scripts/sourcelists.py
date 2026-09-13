# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Functions to handle source lists and cross-matching of 2 or 3 source lists
#
# 10 Dec 2009: Mathieu Servillat

import calendar
import os
import time

from astropy.io import ascii, fits
from astropy.table import Table
import numpy as np

import ldac_tools

# ---------- GENERAL ----------


def separation(ra1, dec1, ra2, dec2):
    """
    Give the separation in degrees between 2 points (ra1,dec1) or between a source list (ra1,dec1) and a reference point (ra2,dec2)

    Computation:
    given points: 1,2 and let point 3 be the north pole then the spherical astronomy cosine law is:
        cosc=cosa * cosb + sina * sinb * cosC (spherical astronomy, smart pg8)
    where:
        A (pnt 1) angle bewteen pnt 2,3
        B (pnt 2) angle between pnt 1,3
        C (pnt 3) angle between 1,2 just abs(ra1-ra2)
        a - distance  between b,Pole just (90- dec2)
        b - distance  between a,Pole just (90 -dec1)
        c - distance  between a,b  (what we want...)
    """
    ra1 = np.array(ra1)
    dec1 = np.array(dec1)
    anglC = np.deg2rad(ra1 - ra2)
    dista = np.deg2rad(90 - dec1)
    distb = np.deg2rad(90 - dec2)
    cosc = np.cos(dista) * np.cos(distb) + np.sin(dista) * np.sin(distb) * np.cos(anglC)
    distc = np.rad2deg(np.arccos(cosc))
    return distc


def ra_shift(dec1, dist):
    """
    Give RA shift in degrees equivalent to distance dist in degrees for a given ra1 and dec1
    Remark: ra1 is not necessary for this

    Computation:
    given points: 1,2 and let point 3 be the north pole then the spherical astronomy cosine law is:
        cosc=cosa * cosb + sina * sinb * cosC (spherical astronomy, smart pg8)
    where:
        point 2 is at dec2 = dec1 and ra2 = ra1 + ra_shift where ra_shift corresponds to the given distance dist
        A (pnt 1) angle bewteen pnt 2,3
        B (pnt 2) angle between pnt 1,3
        C (pnt 3) angle between 1,2 (what we want...)
        a - distance  between b,Pole just (90- dec1)
        b - distance  between a,Pole just (90 -dec1)
        c - distance  between a,b given as dist
    """
    # dista = np.deg2rad(90-dec1)
    # distb = dista
    # distc = np.deg2rad(dist)
    # cosC = (np.cos(distc)-np.cos(dista)*np.cos(distb)) / (np.sin(dista)*np.sin(distb))
    # anglC = np.rad2deg(np.arccos(cosC))
    # return anglC
    return (dist) / np.cos(np.deg2rad(dec1))


# ---------- READ/WRITE ----------


def save_ds9reg(
    filename,
    ra,
    dec,
    errors=0.5,
    minerror=0.05,
    color="green",
    width=1,
    text="",
    tag="",
    header="fk5",
):
    """
    Save ra and dec as DS9 region file
    """
    n = len(ra)
    if type(errors) == type(1.0):
        errors = np.array(n * [errors])
    sub = np.where(errors < minerror)
    if len(sub) > 0:
        errors[sub] = minerror
    if type(color) == type(""):
        color = n * [color]
    if type(width) == type(1):
        width = n * [width]
    if type(text) == type(""):
        text = n * [text]
    if type(tag) == type(""):
        tag = n * [tag]
    fout = open(filename, "w")
    fout.write(header + "\n")
    for i in range(n):
        t = ""
        if width[i] != 1:
            t += " width=" + str(width[i]) + ""
        if text[i] != "":
            t += " text={" + str(text[i]) + "}"
        if tag[i] != "":
            t += " tag={" + str(tag[i]) + "}"
        fout.write(
            "circle("
            + str(ra[i])
            + ","
            + str(dec[i])
            + ","
            + str(errors[i])
            + '") # color="'
            + color[i]
            + '"'
            + t
            + "\n"
        )
    fout.close()


def save_ds9reg_box(
    filename,
    ra,
    dec,
    h=1.0,
    w=1.0,
    color="green",
    width=1,
    text="",
    tag="",
    header="fk5",
):
    """
    Save ra and dec as DS9 region file with box(), h and w in arcmin
    text is added as a text() in the center of the box
    """
    n = len(ra)
    if type(h) == type(1.0):
        h = np.array(n * [h])
    if type(w) == type(1.0):
        w = np.array(n * [w])
    if type(color) == type(""):
        color = n * [color]
    if type(width) == type(1):
        width = n * [width]
    if type(text) == type(""):
        text = n * [text]
    if type(tag) == type(""):
        tag = n * [tag]
    fout = open(filename, "w")
    fout.write(header + "\n")
    for i in range(n):
        t = ""
        if width[i] != 1:
            t += " width=" + str(width[i]) + ""
        if tag[i] != "":
            t += " tag={" + str(tag[i]) + "}"
        r = ra[i]
        if ":" in r:
            r = hms2deg(r)
        d = dec[i]
        if ":" in d:
            d = dms2deg(d)
        fout.write(
            "box("
            + str(r)
            + ","
            + str(d)
            + ","
            + str(h[i])
            + "',"
            + str(w[i])
            + "') # color=\""
            + color[i]
            + '"'
            + t
            + "\n"
        )
        if text[i] != "":
            fout.write(
                "# text("
                + str(r)
                + ","
                + str(d + (h[i] / 4.0) / 60.0)
                + ") text={"
                + str(text[i])
                + "}\n"
            )
    fout.close()


def save_ds9reg_line(
    filename, ra1, dec1, ra2, dec2, color="green", text="", tag="", width=0
):
    """
    Save ra and dec as DS9 region file with line()
    """
    n = len(ra1)
    if type(color) == type(""):
        color = n * [color]
    if type(text) == type(""):
        text = n * [text]
    if type(tag) == type(""):
        tag = n * [tag]
    if type(width) == type(0):
        width = n * [width]
    fout = open(filename, "w")
    fout.write("fk5\n")
    for i in range(n):
        t = ""
        if text[i] != "":
            t += " text={" + str(text[i]) + "}"
        if tag[i] != "":
            t += " tag={" + str(tag[i]) + "}"
        if width[i] != 0:
            t += " width=" + str(width[i])
        fout.write(
            "line("
            + str(ra1[i])
            + ","
            + str(dec1[i])
            + ","
            + str(ra2[i])
            + ","
            + str(dec2[i])
            + ') # color="'
            + color[i]
            + '"'
            + t
            + "\n"
        )
    fout.close()


def read_ds9reg(filename, shape="circle", radius=True, type="default"):
    """
    Save ra and dec as DS9 region file
    If coordinates are in h:m:s d:m:s format it is converted to degrees
    type='atpy' gives an Astropy table (this was implemented with `atpy` back in the day)
    return ra,dec,r,text
    """
    fin = open(filename)
    ra = []
    dec = []
    r = []
    r2 = []
    text = []
    for line in fin.readlines():
        if shape in line:
            coords = line.split(")")[0].split("(")[-1].split(",")
            # print coords
            if ":" in coords[0]:
                rai = hms2deg(coords[0].strip())
            else:
                rai = float(coords[0])
            if ":" in coords[1]:
                deci = dms2deg(coords[1].strip())
            else:
                deci = float(coords[1])
            ra.append(rai)
            dec.append(deci)
            if '"' in coords[2]:
                ri = float(coords[2][0:-1]) / 3600.0
            else:
                if "'" in coords[2]:
                    ri = float(coords[2][0:-1]) / 60.0
                else:
                    ri = float(coords[2])
            r.append(ri)
            if shape == "ellipse" or shape == "box":
                if '"' in coords[3]:
                    r2i = float(coords[3][0:-1]) / 3600.0
                else:
                    if "'" in coords[3]:
                        r2i = float(coords[2][0:-1]) / 60.0
                    else:
                        r2i = float(coords[3])
                r2.append(r2i)
            if "text={" in line:
                text.append(line.split("text={")[-1].split("}")[0])
            else:
                text.append("")
    if shape == "ellipse":
        r = list(np.mean(list(zip(np.array(r), np.array(r2))), axis=1))
    if type == "atpy":
        t = Table()
        t.add_column("RA", ra)
        t.add_column("Dec", dec)
        t.add_column("err", r)
        t.add_column("text", text)
        return t
    else:
        return ra, dec, r, text


def read_rdb(filename, delimiter="\t", comments="#", data_start=2, verbose=False):
    """
    Read a RDB acsii file (see starbase project)
    """
    tab = Table()
    tab.table_name = filename.split("/")[-1]
    data = ascii.read(
        filename,
        data_start=data_start,
        delimiter=delimiter,
    )
    for c in data.dtype.names:
        tab.add_column(c, data[c])
    if verbose:
        tab.describe()
    return tab


def write_refcat_scamp(filename, ra, dec, erra, errb, mag, magerr):
    data = Table()
    data.add_column("X_WORLD", ra, unit="deg")
    data.add_column("Y_WORLD", dec, unit="deg")
    data.add_column("ERRA_WORLD", erra, dtype="float32", unit="deg")
    data.add_column("ERRB_WORLD", errb, dtype="float32", unit="deg")
    data.add_column("MAG", mag, dtype="float32", unit="mag")
    data.add_column("MAGERR", magerr, dtype="float32", unit="mag")

    hdu_out = ldac_tools.table_as_ldac(fits.PrimaryHDU().header, data)
    hdu_out.writeto(filename, overwrite=True)


def read_obscat_magellan(filename, verbose=False):
    """
    Read a Magellan obs catalog from acsii file into Astropy table
    """
    return ascii.read(
        filename,
        names=[
            "num",
            "name",
            "RA",
            "Dec",
            "epoch",
            "RApm",
            "Decpm",
            "offset",
            "rot",
            "RA_probe1",
            "Dec_probe1",
            "epoch_probe1",
            "RA_probe2",
            "Dec_probe2",
            "epoch_probe2",
        ],
        data_start=0,
        comment="#",
        delimiter=" ",
    )


# ---------- CONVERT COORD ----------


def isarray(obj):
    if hasattr(obj, "shape"):
        s = obj.shape
        if len(s) > 0:
            return True
        else:
            return False
    elif hasattr(obj, "pop"):
        return True
    return False


# Convert HH:MM:SS.SSS into Degrees :
def hms2deg(ra, separator=":"):
    if isarray(ra):
        result = []
        for rai in ra:
            result.append(_hms2deg(rai, separator=separator))
        return np.array(result, dtype="float64")
    else:
        return _hms2deg(ra, separator=separator)


def _hms2deg(ra, separator=":"):
    try:
        sep1 = ra.find(separator)
        hh = int(ra[0:sep1])
        sep2 = ra[sep1 + 1 :].find(separator)
        mm = int(ra[sep1 + 1 : sep1 + sep2 + 1])
        ss = float(ra[sep1 + sep2 + 2 :])
    except:
        raise
    else:
        pass
    return hh * 15.0 + mm / 4.0 + ss / 240.0


# Convert +DD:MM:SS.SSS into Degrees :
def dms2deg(dec, separator=":"):
    if isarray(dec):
        result = []
        for deci in dec:
            result.append(_dms2deg(deci, separator=separator))
        return np.array(result, dtype="float64")
    else:
        return _dms2deg(dec, separator=separator)


def _dms2deg(dec, separator=":"):
    Csign = dec[0]
    if Csign == "-":
        sign = -1.0
        off = 1
    elif Csign == "+":
        sign = 1.0
        off = 1
    else:
        sign = 1.0
        off = 0
    try:
        sep1 = dec.find(separator)
        deg = int(dec[off:sep1])
        sep2 = dec[sep1 + 1 :].find(separator)
        arcmin = int(dec[sep1 + 1 : sep1 + sep2 + 1])
        arcsec = float(dec[sep1 + sep2 + 2 :])
    except:
        raise
    else:
        pass
    return sign * (deg + (arcmin * 5.0 / 3.0 + arcsec * 5.0 / 180.0) / 100.0)


# Convert RA (deg) to H.M.S:
def deg2hms(RAin, separator=":", plus="", digits=3):
    if isarray(RAin):
        result = []
        for RAini in RAin:
            result.append(_deg2hms(RAini, separator=separator, plus=plus))
        return np.array(result)
    else:
        return _deg2hms(RAin, separator=separator, plus=plus)


def _deg2hms(RAin, separator=":", plus="", digits=3):
    if RAin < 0:
        sign = -1
        ra = -RAin
    else:
        sign = 1
        ra = RAin

    h = int(ra / 15.0)
    ra -= h * 15.0
    m = int(ra * 4.0)
    ra -= m / 4.0
    s = ra * 240.0

    if sign == -1:
        out = (
            ("-%02d" % h)
            + separator
            + ("%02d" % m)
            + separator
            + str(round(s, digits)).zfill(3 + digits)
        )
    else:
        out = (
            plus
            + ("%02d" % h)
            + separator
            + ("%02d" % m)
            + separator
            + str(round(s, digits)).zfill(3 + digits)
        )

    return out


# Convert Decl. (deg) to D.M.S:
def deg2dms(Decin, separator=":", plus="", digits=2):
    if isarray(Decin):
        result = []
        for Decini in Decin:
            result.append(_deg2dms(Decini, separator=separator, plus=plus))
        return np.array(result)
    else:
        return _deg2dms(Decin, separator=separator, plus=plus)


def _deg2dms(Decin, separator=":", plus="", digits=2):
    if Decin < 0:
        sign = -1
        dec = -Decin
    else:
        sign = 1
        dec = Decin

    d = int(dec)
    dec -= d
    dec *= 100.0
    m = int(dec * 3.0 / 5.0)
    dec -= m * 5.0 / 3.0
    s = dec * 180.0 / 5.0

    if sign == -1:
        out = (
            ("-%02d" % d)
            + separator
            + ("%02d" % m)
            + separator
            + str(round(s, digits)).zfill(3 + digits)
        )
    else:
        out = (
            plus
            + ("%02d" % d)
            + separator
            + ("%02d" % m)
            + separator
            + str(round(s, digits)).zfill(3 + digits)
        )

    return out


# Convert YYYY-MM-DD into MMJD :
def dateobs2mmjd(dateobs):
    tuple_time = time.strptime(dateobs, "%Y-%m-%d")
    dateobs_as_seconds = time.mktime(tuple_time)
    refepoch = time.strptime("2003-01-01", "%Y-%m-%d")
    refepoch_as_seconds = time.mktime(refepoch)

    return (dateobs_as_seconds - refepoch_as_seconds) / 3600.0 / 24.0 + 1


# Convert YYYY-MM-DDTHH:MM:SS.SSS into MMJD :
def timeobs2mmjd(timeobs):
    day = timeobs[0 : timeobs.find("T")]
    time = timeobs[timeobs.find("T") + 1 :]
    try:
        hour = int(time[0:2])
        min = int(time[3:5])
        sec = float(time[6:])
    except:
        raise Exception(f"Format: {timeobs}")
    else:
        mmjd = dateobs2mmjd(day) + (((sec / 60.0) + min) / 60.0 + hour) / 24.0

    return mmjd


# Convert YYYY-MM-DDTHH:MM:SS.SSS into UTC :
def timeobs2secUT(timeobs):
    try:
        tuple = time.strptime(timeobs[0:19], "%Y-%m-%dT%H:%M:%S")
    except:
        raise Exception(f"Format: {timeobs}")
    else:
        secs = time.mktime(tuple)

    return secs


# MMJD -> Date YYYY-MM-DD & HH:MM
def mmjd2date(mmjd):
    refepoch = time.strptime("2003-01-01", "%Y-%m-%d")
    refepoch_as_seconds = time.mktime(refepoch)
    epoch = refepoch_as_seconds + mmjd * 24.0 * 60.0 * 60.0
    tuple_date = time.gmtime(epoch)
    date = time.strftime("%Y-%m-%dT%H:%M", tuple_date)

    day = date[: date.find("T")]
    hour = date[date.find("T") + 1 :]

    return [day, hour]


# Cut YYYY-MM-DDTHH:MM:SS.SSS into [Y,M,D,H,M,S]:
def splitESOdate(timeobs):
    date = timeobs[0 : timeobs.find("T")]
    time = timeobs[timeobs.find("T") + 1 :]
    try:
        year = int(date[:4])
        mounth = int(date[5:7])
        day = int(date[8:10])
        hour = int(time[0:2])
        min = int(time[3:5])
        sec = float(time[6:])
    except:
        print("! Bad Obs Date Format : %s ." % timeobs)
        raise
    else:
        pass

    return [year, mounth, day, hour, min, sec]


# Return the ObsDate at ESO convention :  date of beginning of night.
# Use GMT to handle change in month :
def getESOnight(timeobs):
    try:
        [y, m, d, h, min, sec] = splitESOdate(timeobs)
    except:
        print("! Not ESO date format : %s " % timeobs)
        raise
    else:
        if h < 12:
            corr = -1.0
        else:
            corr = 0.0

    day = "%s-%s-%s" % (y, m, d)
    tuple = time.strptime(day, "%Y-%m-%d")
    sec = calendar.timegm(tuple) + corr * 24.0 * 60.0 * 60.0
    ESOtime = time.gmtime(sec)
    ESOday = time.strftime("%Y-%m-%d", ESOtime)

    return ESOday


# Return the hour modulus 12. from HH:MM:SS :
def ESOhour(hour):
    try:
        h = int(hour[:2])
        m = int(hour[3:5])
        s = int(hour[6:8])
    except:
        print("! Format error : %s ." % hour)
        raise

    ESOhour = h + m / 60.0 + s / 3600.0

    if ESOhour > 12.0:
        ESOhour -= 24.0

    return ESOhour


# MIDAS MJD @ MMJD=0 :
mmjd2mjd = 52640.0
mmjd2gmt = 12052.0


# ---------- MATCHING ----------


def match2(
    ra1,
    dec1,
    ra2,
    dec2,
    errors=0.2,
    nsigma=1.0,
    radius=0.2,
    onlybest=False,
    bsra=0,
    bsdec=0,
):
    """
    Match 2 (ra,dec) source lists and return the indices and separations in arcsec
    Give radius or individual errors for the second list
    Should be faster if the second sourcelist is shorter than the first

    Input:
    - ra1, dec1: first source list with RA and Dec in degrees
    - ra2, dec2: second source list with RA and Dec in degrees
    - radius: distance in arcsec for the most distant match
    - errors: distance in arcsec for the most distant match for each source in the second list, if given, radius won't be used
    - onlybest: if True, will keep only the best match, if False will list all the matches within the given error or radius
    - bsra, bsdec: boresight correction in arcsec

    Output:
    - ind1, ind2: indices for each source list for successfull matches between the 2 source lists
    - dist: array with same number of elements with the separations in arcsec
    """
    from scipy.spatial import KDTree

    ra1 = np.array(ra1)
    dec1 = np.array(dec1)
    ra2 = np.array(ra2)
    dec2 = np.array(dec2)
    n1 = len(ra1)
    n2 = len(ra2)
    decmd = np.median(np.append(dec1, dec2))
    decmd_rad = np.deg2rad(decmd)
    # ra_factor =  (1/60) / ra_shift(decmd,1/60)
    # apply boresight correction
    bsra_shift = 0.0
    bsdec_shift = 0.0
    if bsra != 0:
        bsra_shift = (bsra / 3600) / np.cos(decmd_rad)
    if bsdec != 0:
        bsdec_shift += bsdec / 3600
    # print '\nMatching 2 lists of',n1,'and',n2,'elements'
    tree1 = KDTree(list(zip(ra1 * np.cos(decmd_rad), dec1)))
    pts = np.array(
        list(zip((ra2 + bsra_shift) * np.cos(decmd_rad), dec2 + bsdec_shift))
    )
    # tree2=KDTree(zip(ra2*ra_factor,dec2))
    if errors == "none":
        errors = radius
    if "float" in str(type(errors)):
        errors = n2 * [float(errors)]
    errors = np.array(errors)
    if onlybest:
        # print 'Keep only best matches'
        ind1 = []
        ind2 = []
        dist12 = []
        for i in range(len(pts)):
            md, mi = tree1.query(pts[i], distance_upper_bound=errors[i] / 3600)
            if mi < n1:
                ind1.append(mi)
                ind2.append(i)
                dist12.append(separation(ra1[mi], dec1[mi], ra2[i], dec2[i]) * 3600)
    else:
        ind1 = []
        ind2 = []
        dist12 = []
        for i in range(len(pts)):
            mi = [0]
            k = 2
            while mi[-1] != n1:
                md, mi = tree1.query(pts[i], k=k, distance_upper_bound=errors[i] / 3600)
                k += 2
            for j in range(len(mi)):
                if mi[j] < n1:
                    ind1.append(mi[j])
                    ind2.append(i)
                    dist12.append(
                        separation(ra1[mi[j]], dec1[mi[j]], ra2[i], dec2[i]) * 3600
                    )

    print("Found", len(ind1), "matches")
    return np.array(ind1), np.array(ind2), np.array(dist12)


# ----------
def match2_boresight(
    ra1,
    dec1,
    ra2,
    dec2,
    errors=0.2,
    nsigma=1.0,
    bsra_init=0.0,
    bsdec_init=0.0,
    outfile=False,
    nsigmasub=None,
    labels=None,
):
    """
    Call match2 and analyse the results to estimate the boresight shift and error from the weighted mean of the difference ra1-ra2 and dec1-dec2

    Input:
    - ra1, dec1: first source list with RA and Dec in degrees
    - ra2, dec2: second source list with RA and Dec in degrees
    - errors: 1 sigma error in arcsec for each source in the second list
    - nsigma: multiply errors by nsigma for matching
    - bsra_init: initial RA correction in arcsec for the second source list (ra2+bsra_init_shift)
    - bsdec_init: initial Dec correction in arcsec for the second source list (dec2+bsdec_init)

    Output:
    - ind1, ind2: indices for each source list for successfull matches between the 2 source lists
    - dist: array with same number of elements with the separations in arcsec
    - bsra, bsdec: boresight corection in arcsec
    - bsrasig, bsdecsig: 1 sigma error on boresight correction
    """
    # correct rashift distance to rashiftd angle !!! Necessary? give angle as input for now
    ra1 = np.array(ra1)
    dec1 = np.array(dec1)
    ra2 = np.array(ra2)
    dec2 = np.array(dec2)
    n1 = len(ra1)
    n2 = len(ra2)
    bsra_init_shift = 0.0
    decmd = np.median(np.append(dec1, dec2))
    decmd_rad = np.deg2rad(decmd)
    if bsra_init != 0:
        bsra_init_shift = (bsra_init / 3600.0) / np.cos(decmd_rad)
    if "float" in str(type(errors)):
        errors = n2 * [float(errors)]
    errors = np.array(errors)
    if nsigmasub == None:
        nsigmasub = nsigma
    # match source lists
    ind1, ind2, dist12 = match2(
        ra1,
        dec1,
        ra2,
        dec2,
        errors=errors * nsigma,
        onlybest=True,
        bsra=bsra_init,
        bsdec=bsdec_init,
    )
    # weighted mean for RA and Dec
    nxi = len(ind2)
    bsra = 0
    bsdec = 0
    bsrasig = 0
    bsdecsig = 0
    if nxi > 0:
        wi = 1.0 / errors[ind2] ** 2
        sumwi = np.sum(wi)
        radist_shift = ra1[ind1] - (ra2[ind2] + bsra_init_shift)
        radist = radist_shift * np.cos(decmd_rad)
        bsra_shift = np.sum(wi * radist_shift) / sumwi
        bsra = bsra_shift * np.cos(decmd_rad)
        decdist = dec1[ind1] - (dec2[ind2] + bsdec_init / 3600.0)
        bsdec = np.sum(wi * decdist) / sumwi
        # sigma for weighted mean
        bsrasig_shift = np.sqrt(
            np.sum(wi * (radist_shift - bsra_shift) ** 2) / ((nxi - 1) * sumwi)
        )
        bsrasig = bsrasig_shift * np.cos(decmd_rad) * 3600.0
        bsdecsig = (
            np.sqrt(np.sum(wi * (decdist - bsdec) ** 2) / ((nxi - 1) * sumwi))
        ) * 3600.0
        # bserr = np.sqrt(bsrasig**2+bsdecsig**2)
        bsra *= 3600.0
        bsdec *= 3600.0
        print(
            "Found boresight correction: "
            + str(round(bsra, 5))
            + " +/- "
            + str(round(bsrasig, 5))
            + ", "
            + str(round(bsdec, 5))
            + " +/- "
            + str(round(bsdecsig, 5))
            + " arcsec"
        )
        if outfile:
            import matplotlib.pyplot as plt

            fig = plt.figure(figsize=(8, 8))
            ax = plt.subplot(111)
            titadd = ""
            if (bsra_init != 0.0) or (bsdec_init != 0):
                titadd = " (bsra=" + str(bsra_init) + ",bsdec=" + str(bsdec_init) + ")"
            ax.set_title(os.path.splitext(outfile)[0] + titadd)
            p1 = ax.plot(radist * 3600.0, decdist * 3600.0, ".b", ms=10)
            if nsigmasub != nsigma:
                ind1s, ind2s, dist12s = match2(
                    ra1[ind1],
                    dec1[ind1],
                    ra2[ind2],
                    dec2[ind2],
                    errors=errors[ind2] * nsigmasub,
                    onlybest=True,
                    bsra=bsra_init,
                    bsdec=bsdec_init,
                )
                p2 = ax.plot(
                    radist[ind1s] * 3600.0, decdist[ind1s] * 3600.0, ".r", ms=10
                )
                ax.legend(
                    (p1, p2),
                    (
                        "< " + str(nsigma) + " sigma matches",
                        "< " + str(nsigmasub) + " sigma matches",
                    ),
                    "upper right",
                    numpoints=1,
                    prop={"size": 11},
                )
            else:
                ax.legend(
                    (p1),
                    ("< " + str(nsigma) + " sigma matches"),
                    "upper right",
                    numpoints=1,
                    prop={"size": 11},
                )
            xymax = (
                round(max(np.append(np.abs(radist), np.abs(decdist))) * 36000.0) / 10.0
                + 0.1
            )
            lim = [-xymax, xymax]
            xl = ax.set_xlabel("RA shift (arcsec)")
            yl = ax.set_ylabel("Dec shift (arcsec)")
            plx = ax.plot([0, 0], lim, "--k")
            ply = ax.plot(lim, [0, 0], "--k")
            xlim = ax.set_xlim(lim)
            ylim = ax.set_ylim(lim)
            g = ax.grid()
            plt.savefig(outfile)

    return ind1, ind2, dist12, bsra, bsdec, bsrasig, bsdecsig, nxi


# ----------
def match2_boresight_loop(
    ra1,
    dec1,
    ra2,
    dec2,
    errors=0.2,
    nsigma=1.0,
    bsra=0.0,
    bsdec=0.0,
    bserr=0.4,
    precision=0.001,
):
    """
    Iterate to get boresight correction with given precision

    Input:
    - ra1, dec1: first source list with RA and Dec in degrees
    - ra2, dec2: second source list with RA and Dec in degrees
    - errors: 1 sigma error in arcsec for each source in the second list
    - nsigma: multiply errors by nsigma for matching
    - bsra: initial RA correction in arcsec for the second source list
    - bsdec: initial Dec correction in arcsec for the second source list
    - bserr: initial boresight 1 signa error in arcsec for the second source list (i.e. the absolute error). To be added to error (which should contain all the other relative errors)
    - precision: refine boresight until the correction is lower than precision

    Output:
    - ind1, ind2: indices for each source list for successfull matches between the 2 source lists
    - dist: array with same number of elements with the separations in arcsec
    - bsra_tot, bsdec_tot: total boresight corection in arcsec
    - bsrasig, bsdecsig: 1 sigma error on boresight correction
    """
    ra1 = np.array(ra1)
    dec1 = np.array(dec1)
    ra2 = np.array(ra2)
    dec2 = np.array(dec2)
    n1 = len(ra1)
    n2 = len(ra2)
    if "float" in str(type(errors)):
        errors = n2 * [float(errors)]
    errors = np.array(errors)
    bsra_tot = bsra
    bsdec_tot = bsdec
    bsra = 1
    bsdec = 1
    i = 0
    while (abs(bsra) > 0.001 or abs(bsdec) > 0.001) and i < 50:
        ind1, ind2, dist12, bsra, bsdec, bsrasig, bsdecsig, nxi = match2_boresight(
            ra1,
            dec1,
            ra2,
            dec2,
            errors=np.sqrt(errors**2 + bserr**2),
            nsigma=nsigma,
            bsra_init=bsra_tot,
            bsdec_init=bsdec_tot,
        )
        bsra_tot += bsra * 0.8
        bsdec_tot += bsdec * 0.8
        print(
            "Total boresight correction: "
            + str(round(bsra_tot, 5))
            + " +/- "
            + str(round(bsrasig, 5))
            + ", "
            + str(round(bsdec_tot, 5))
            + " +/- "
            + str(round(bsdecsig, 5))
            + " arcsec"
        )
        bserr = np.sqrt(bsrasig**2 + bsdecsig**2)
        print("New Boresight error:", bserr)
        i += 1
    return ind1, ind2, dist12, bsra_tot, bsdec_tot, bsrasig, bsdecsig, nxi


# ----------
def match2_boresight_map(
    ra1,
    dec1,
    ra2,
    dec2,
    errors=0.2,
    nsigma=1.0,
    maxsh=5.0,
    incr=0.5,
    bsra=0.0,
    bsdec=0.0,
    outfile="",
):
    """
    create a map of the number of matches with given offsets

    Input:
    - ra1, dec1: first source list with RA and Dec in degrees
    - ra2, dec2: second source list with RA and Dec in degrees
    - errors: 1 sigma error in arcsec for the most distant match for each source in the second list
    - nsigma: multiply errors by nsigma for matching
    - maxsh: maximum shift for matching
    - incr: steps for shifts
    - bsra: boresight correction in RA
    - bsdec: boresight correction in Dec
    - outfile: if given, save plot with this filename

    Output:
    X, Y: coordinates of the map
    Z: number of matches for each coordinate
    """
    shifts = np.array(list(range(2 * maxsh / incr + 1))) * incr - maxsh
    szmap = len(shifts)
    print("Size of map:", szmap, "x", szmap)
    print("Shifts:")
    i = 0
    bmap = np.zeros([szmap, szmap])
    for i, xi in enumerate(shifts):
        for j, yj in enumerate(shifts):
            print("\n", i + 1, "/", szmap, "-", j + 1, "/", szmap)
            i1, i2, d12, bsra_new, bsdec_new, bsrasig, bsdecsig, nxi = match2_boresight(
                ra1,
                dec1,
                ra2,
                dec2,
                errors=errors,
                nsigma=nsigma,
                bsra_init=xi + bsra,
                bsdec_init=yj + bsdec,
            )
            bmap[i, j] = nxi
    X, Y = np.meshgrid(shifts, shifts)
    Z = bmap
    if outfile != "":
        print("Saving " + outfile)
        from matplotlib import cm
        import matplotlib.pyplot as plt

        fig = plt.figure()
        plt.clf()
        imsz = maxsh + 0.5 * incr
        im = plt.imshow(
            Z,
            interpolation="hanning",
            origin="lower",
            extent=(-imsz, imsz, -imsz, imsz),
            cmap=cm.gist_heat,
        )
        plt.grid(c="white")
        cb = plt.colorbar()
        plt.xlabel("RA shift (arcsec)")
        plt.ylabel("Dec shift (arcsec)")
        cb.set_label("N of matches")
        mx = np.max(Z)
        lev = [0.5 * mx, 0.8 * mx, 0.95 * mx]
        CS = plt.contour(X, Y, Z, lev, cmap=cm.Greys)
        # plt.show()
        plt.savefig(outfile)

    return X, Y, Z


# ----------
def match2_separation(ra1, dec1, ra2, dec2, ind):
    """
    Give the separation in arcsec between 2 source lists with the same number of elements and same order

    Input:
    - ra1, dec1: first source list with RA and Dec in degrees
    - ra2, dec2: second source list with RA and Dec in degrees

    Output:
    - dist: array with same number of elements with the separations in arcsec
    """
    ra1 = np.array(ra1)
    dec1 = np.array(dec1)
    ra2 = np.array(ra2)
    dec2 = np.array(dec2)
    dist = []
    for i in range(len(ind)):
        indi = ind[i]
        disti = []
        for ii in indi:
            disti.append(separation(ra1[ii], dec1[ii], ra2[i], dec2[i]) * 3600)
        dist.append(disti)
    return dist


# ----------
def match3(ra1, dec1, ra2, dec2, ra3, dec3, radius=0.6, names=["1", "2", "3"]):
    """
    Match 3 (ra,dec) source lists and return the indices
    Only the best matches withion the circle with given radius are kept

    Input:
    - ra1, dec1: first source list with RA and Dec in degrees
    - ra2, dec2: second source list with RA and Dec in degrees
    - ra3, dec3: third source list with RA and Dec in degrees
    - radius: distance in arcsec for the most distant match

    Output:
    - ind1_123, ind2_123 and ind3_123: indices for each source list for successfull cross match of the 3 source lists
    - dist_123: distance indicator for 3 match sources
    - ind1_12x and ind2_12x: indices for source list 1 and 2 for godd matches only between 1 and 2, and not 3 (1-2 exclusive matches)
    - ind2_23x and ind3_23x: same with source lists 2 and 3
    - ind3_31x and ind1_31x: same with source lists 3 and 1
    - ind1x, ind2x and ind3x: sources without match in each source list
    """
    ra1 = np.array(ra1)
    dec1 = np.array(dec1)
    ra2 = np.array(ra2)
    dec2 = np.array(dec2)
    ra3 = np.array(ra3)
    dec3 = np.array(dec3)
    n1 = len(ra1)
    n2 = len(ra2)
    n3 = len(ra3)
    double = [names[0] + names[1], names[1] + names[2], names[2] + names[0]]
    triple = names[0] + names[1] + names[2]
    print(names)
    print(double)
    print(triple)

    # Match 1-2, 2-3 and 3-1
    print("\nMatching 3 source lists of", str(n1) + ",", n2, "and", n3, "elements")
    print("\n* Cross-match " + names[0] + "-" + names[1])
    ind1_12, ind2_12, dist12 = match2(
        ra1, dec1, ra2, dec2, radius=radius, onlybest="yes"
    )
    print("\n* Cross-match " + names[1] + "-" + names[2])
    ind2_23, ind3_23, dist23 = match2(
        ra2, dec2, ra3, dec3, radius=radius, onlybest="yes"
    )
    print("\n* Cross-match " + names[2] + "-" + names[0])
    ind3_31, ind1_31, dist31 = match2(
        ra3, dec3, ra1, dec1, radius=radius, onlybest="yes"
    )

    ind1_123 = np.array([], dtype="int64")
    ind2_123 = np.array([], dtype="int64")
    ind3_123 = np.array([], dtype="int64")
    dist_123 = np.array([], dtype="int64")
    # Match 12-3
    print("\n* Cross-match " + double[0] + "-" + names[2])
    ra12 = (ra1[ind1_12] + ra2[ind2_12]) / 2
    dec12 = (dec1[ind1_12] + dec2[ind2_12]) / 2
    ind3_12, ind12_3, dist3_12 = match2(
        ra3, dec3, ra12, dec12, radius=radius * 1.5, onlybest="yes"
    )
    ind1_123 = np.append(ind1_123, ind1_12[ind12_3])
    ind2_123 = np.append(ind2_123, ind2_12[ind12_3])
    ind3_123 = np.append(ind3_123, ind3_12)
    dist3_123 = np.column_stack((dist3_12, dist12[ind12_3])).mean(axis=1)
    dist_123 = np.append(dist_123, dist3_123)
    # Match 23-1
    print("\n* Cross-match " + double[1] + "-" + names[0])
    ra23 = (ra2[ind2_23] + ra3[ind3_23]) / 2
    dec23 = (dec2[ind2_23] + dec3[ind3_23]) / 2
    ind1_23, ind23_1, dist1_23 = match2(
        ra1, dec1, ra23, dec23, radius=radius * 1.5, onlybest="yes"
    )
    ind1_123 = np.append(ind1_123, ind1_23)
    ind2_123 = np.append(ind2_123, ind2_23[ind23_1])
    ind3_123 = np.append(ind3_123, ind3_23[ind23_1])
    dist1_123 = np.column_stack((dist1_23, dist23[ind23_1])).mean(axis=1)
    dist_123 = np.append(dist_123, dist1_123)
    # Match 31-2
    print("\n* Cross-match " + double[2] + "-" + names[1])
    ra31 = (ra3[ind3_31] + ra1[ind1_31]) / 2
    dec31 = (dec3[ind3_31] + dec1[ind1_31]) / 2
    ind2_31, ind31_2, dist2_31 = match2(
        ra2, dec2, ra31, dec31, radius=radius * 1.5, onlybest="yes"
    )
    ind1_123 = np.append(ind1_123, ind1_31[ind31_2])
    ind2_123 = np.append(ind2_123, ind2_31)
    ind3_123 = np.append(ind3_123, ind3_31[ind31_2])
    dist2_123 = np.column_stack((dist2_31, dist31[ind31_2])).mean(axis=1)
    dist_123 = np.append(dist_123, dist2_123)

    # print dist_123

    # Remove multiple matches
    print("\n* Clean duplicated triple matches")
    ind3x = Table()
    ind3x.add_column(names[0], ind1_123, dtype="int64")
    ind3x.add_column(names[1], ind2_123, dtype="int64")
    ind3x.add_column(names[2], ind3_123, dtype="int64")
    ind3x.add_column("dist", dist_123, dtype="float64")
    for col in names:
        ind3x.sort([col, "dist"])
        keep = []
        i = 0
        n = len(ind3x)
        while i < n:
            v = ind3x[col][i]
            keep.append(i)
            i += 1
            while i < n and v == ind3x[col][i]:
                i += 1
        ind3x = ind3x.rows(keep)
        print(n, "-->", len(ind3x))
    ind1_123 = ind3x[names[0]]
    ind2_123 = ind3x[names[1]]
    ind3_123 = ind3x[names[2]]
    dist_123 = ind3x["dist"]

    print("\nFound", len(ind1_123), "triple matches")
    print("Mean distance =", np.mean(dist_123))

    # 12 23 31 exclusive groups
    print("\n* Compute exclusive complementary groups of double matches")
    ind1_12x, ind2_12x, dist12x = comp2(
        ind1_123, ind2_123, ind1_12, ind2_12, dist=dist12
    )
    print(names[0] + "-" + names[1] + " matches only:", len(ind1_12x))
    ind2_23x, ind3_23x, dist23x = comp2(
        ind2_123, ind3_123, ind2_23, ind3_23, dist=dist23
    )
    print(names[1] + "-" + names[2] + " matches only:", len(ind2_23x))
    ind3_31x, ind1_31x, dist31x = comp2(
        ind3_123, ind1_123, ind3_31, ind1_31, dist=dist31
    )
    print(names[2] + "-" + names[0] + " matches only:", len(ind3_31x))

    # Remove matches in 12 23 31 exclusive groups if already in 123
    print("\n* Clean double matches that are associated to triple matches")
    ind2x = {}
    ind2x[double[0]] = Table()
    ind2x[double[0]].add_column(names[0], ind1_12x, dtype="int64")
    ind2x[double[0]].add_column(names[1], ind2_12x, dtype="int64")
    ind2x[double[0]].add_column("dist", dist12x, dtype="float64")
    ind2x[double[1]] = Table()
    ind2x[double[1]].add_column(names[1], ind2_23x, dtype="int64")
    ind2x[double[1]].add_column(names[2], ind3_23x, dtype="int64")
    ind2x[double[1]].add_column("dist", dist23x, dtype="float64")
    ind2x[double[2]] = Table()
    ind2x[double[2]].add_column(names[2], ind3_31x, dtype="int64")
    ind2x[double[2]].add_column(names[0], ind1_31x, dtype="int64")
    ind2x[double[2]].add_column("dist", dist31x, dtype="float64")
    for t in double:
        for col in [t[0], t[1]]:
            keep = clean_duplicates(ind2x[t][col], ind3x[col], name=col + " in " + t)
            ind2x[t] = ind2x[t].rows(keep)

    # Remove duplicate matches in 12 23 31 exclusive groups
    # if one index is duplicated keep the double match with min distance only
    print("\n* Clean duplicate double matches")
    for t in double:
        for col in [t[0], t[1]]:
            ind2x[t].sort([col, "dist"])
            keep = []
            i = 0
            n = len(ind2x[t])
            while i < n:
                v = ind2x[t][col][i]
                keep.append(i)
                i += 1
                while i < n and v == ind2x[t][col][i]:
                    i += 1
            ind2x[t] = ind2x[t].rows(keep)
            print(n, "-->", len(ind2x[t]))

    # 1 2 3 exclusive groups
    ind1x = {}
    print("\n* Compute exclusive complementary groups of single sources")
    ind1x[names[0]] = comp(
        np.append(
            np.append(ind1_123, ind2x[double[0]][names[0]]),
            ind2x[double[2]][names[0]],
        ),
        list(range(n1)),
    )
    print("sources only detected in " + names[0] + ":", len(ind1x[names[0]]))
    # ind1x[names[1]]=comp(np.append(np.append(ind2_123,ind2_23x),ind2_12x),range(n2))
    ind1x[names[1]] = comp(
        np.append(
            np.append(ind2_123, ind2x[double[1]][names[1]]),
            ind2x[double[0]][names[1]],
        ),
        list(range(n2)),
    )
    print("sources only detected in " + names[1] + ":", len(ind1x[names[1]]))
    # ind1x[names[2]]=comp(np.append(np.append(ind3_123,ind3_31x),ind3_23x),range(n3))
    ind1x[names[2]] = comp(
        np.append(
            np.append(ind3_123, ind2x[double[2]][names[2]]),
            ind2x[double[1]][names[2]],
        ),
        list(range(n3)),
    )
    print("sources only detected in " + names[2] + ":", len(ind1x[names[2]]), "\n")

    # return ind1_123,ind2_123,ind3_123, dist_123, ind1_12x,ind2_12x,ind2_23x,ind3_23x,ind3_31x,ind1_31x, ind1x,ind2x,ind3x
    return ind3x, ind2x, ind1x


def clean_duplicates(ind, indref, name=""):
    keep = []
    for i, ii in enumerate(ind):
        if ii not in indref:
            keep.append(i)
    print(name + ": " + str(len(ind) - len(keep)) + " duplicates")
    return keep


def comp(ind, allind):
    """
    Return the complement of the array ind in the group allind
    """
    comp = np.array([], dtype="int64")
    for v in allind:
        if not v in ind:
            comp = np.append(comp, v)
    return np.array(comp)


def comp2(ind1, ind2, allind1, allind2, dist=None):
    """
    Return the complement (comp1,comp2) of the array (ind1,ind2) in the group (allind1,allind2)
    """
    from scipy.spatial import KDTree

    comp = np.array(list(zip(allind1, allind2)))
    tree1 = KDTree(list(zip(ind1, ind2)))
    n = len(ind1)
    keep = []
    for i, c in enumerate(comp):
        md, mi = tree1.query(c, distance_upper_bound=0.99)
        if mi == n:
            keep.append(i)
    comp = comp[keep]
    # comp=[]
    # ind=zip(ind1,ind2)
    ##allind=np.column_stack((allind1,allind2))
    # allind=zip(allind1,allind2)
    # comp=allind
    # for i,v in enumerate(ind):
    # if v in comp:
    ##comp.append(list(v))
    # del comp[comp.index(v)]
    comp = np.transpose(np.array(comp))
    if dist != None:
        return (
            np.array(comp[0], dtype="int64"),
            np.array(comp[1], dtype="int64"),
            dist[keep],
        )
    else:
        return np.array(comp[0], dtype="int64"), np.array(comp[1], dtype="int64")
