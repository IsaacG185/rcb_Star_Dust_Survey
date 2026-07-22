# plate_scan_utils.py
# Shared, network-free source-detection and plate-quality logic, PLUS
# Gaia/SIMBAD/paradigm identity-crossmatch and a lightweight APASS
# linearity check. Gaia and APASS both go through raw VizieR TSV cone
# searches, matching how 02_A's Cell 3 does it.

from pathlib import Path
import urllib.request
import urllib.parse

from astropy.io import fits as astrofits
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u
from photutils.detection import DAOStarFinder, find_peaks
from photutils.aperture import CircularAperture, aperture_photometry
from photutils.profiles import RadialProfile

from astroquery.simbad import Simbad

import numpy as np
from scipy import ndimage
from skimage.transform import hough_line, hough_line_peaks
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

FWHM = 3.0
THRESHOLD_SIGMA = 5.0
MIN_SOURCES_DAO = 5
BOX_SIZE = 11
APERTURE_RADIUS = 5.0
NAME_MATCH_SEP_ARCSEC = 3.0
GAIA_QUERY_RADIUS_DEG = 0.15
APASS_QUERY_RADIUS_DEG = 0.15
APASS_MATCH_SEP_ARCSEC = 5.0


def detect_sources(fits_path, aperture_radius=APERTURE_RADIUS):
    """Identical logic to 02_A Cell 2. sources['aperture_flux'] is used
    both for pseudo-label training masks and for the linearity check."""
    data = astrofits.getdata(fits_path).astype(float)
    data = np.nan_to_num(data, nan=np.nanmedian(data))
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    data_sub = data - median

    try:
        daofind = DAOStarFinder(fwhm=FWHM, threshold=THRESHOLD_SIGMA * std, sharpness_range=(0.2, 2.0))
        sources = daofind(data_sub)
        algorithm = "DAOStarFinder"
    except Exception:
        sources = None
        algorithm = "DAOStarFinder_failed"

    if sources is None or len(sources) < MIN_SOURCES_DAO:
        sources = find_peaks(data_sub, threshold=THRESHOLD_SIGMA * std, box_size=BOX_SIZE)
        algorithm = "find_peaks"

    if sources is None or len(sources) == 0:
        return None, algorithm, data, data_sub, std, None, None

    if "x_centroid" in sources.colnames:
        x_col, y_col = "x_centroid", "y_centroid"
    else:
        x_col, y_col = "x_peak", "y_peak"

    positions = np.transpose((sources[x_col], sources[y_col]))
    apertures = CircularAperture(positions, r=aperture_radius)
    phot_table = aperture_photometry(data_sub, apertures)
    sources["aperture_flux"] = phot_table["aperture_sum"]
    return sources, algorithm, data, data_sub, std, x_col, y_col


def detect_plate_errors(data):
    """Verbatim from 02_A: scratches, trailing stars, saturation, dust,
    dead zones, edge artifacts. No network calls."""
    errors = {
        "scratches": [], "trailing": [], "saturation": [], "dust": [],
        "edge": False, "dead_zone_fraction": 0.0, "saturation_area_fraction": 0.0,
    }
    h, w = data.shape
    lo, hi = np.nanpercentile(data, 1), np.nanpercentile(data, 99)
    if hi == lo:
        return errors
    norm = np.clip((data - lo) / (hi - lo), 0, 1)

    try:
        scale = 6
        small = ndimage.zoom(norm, 1 / scale, order=1)
        sobel_h = ndimage.sobel(small, axis=0)
        sobel_v = ndimage.sobel(small, axis=1)
        edges = np.hypot(sobel_h, sobel_v)
        edges = (edges > np.percentile(edges, 97)).astype(np.uint8)
        tested_angles = np.linspace(-np.pi / 2, np.pi / 2, 90, endpoint=False)
        hspace, angles, dists = hough_line(edges, theta=tested_angles)
        peaks = hough_line_peaks(hspace, angles, dists, num_peaks=8, min_distance=20, threshold=0.35 * hspace.max())
        hspace_mean, hspace_std = float(np.mean(hspace)), float(np.std(hspace))
        sig_thresh = hspace_mean + 6 * hspace_std
        sh, sw = small.shape
        for peak_val, angle, dist in zip(*peaks):
            if peak_val < sig_thresh:
                continue
            cos_a, sin_a = np.cos(angle), np.sin(angle)
            if abs(sin_a) > 1e-6:
                x0_s, x1_s = 0, sw - 1
                y0_s = (dist - x0_s * cos_a) / sin_a
                y1_s = (dist - x1_s * cos_a) / sin_a
            else:
                y0_s, y1_s = 0, sh - 1
                x0_s = x1_s = dist / cos_a if abs(cos_a) > 1e-6 else 0
            y0_s = float(np.clip(y0_s, 0, sh - 1)); y1_s = float(np.clip(y1_s, 0, sh - 1))
            x0_s = float(np.clip(x0_s, 0, sw - 1)); x1_s = float(np.clip(x1_s, 0, sw - 1))
            errors["scratches"].append({"x0": x0_s * scale, "y0": y0_s * scale, "x1": x1_s * scale, "y1": y1_s * scale})
    except Exception as e:
        print(f"[SCRATCH DETECT] {e}")

    try:
        bright_mask = (norm > np.percentile(norm, 95)).astype(np.uint8)
        labeled, n_obj = ndimage.label(bright_mask)
        for obj_id in range(1, n_obj + 1):
            region = labeled == obj_id
            area = region.sum()
            if area < 40 or area > 0.005 * h * w:
                continue
            coords = np.argwhere(region)
            if len(coords) < 10:
                continue
            cov = np.cov(coords[:, 1], coords[:, 0])
            eigs = np.sort(np.abs(np.linalg.eigvalsh(cov)))
            if eigs[0] < 1e-6:
                continue
            elongation = np.sqrt(eigs[1] / eigs[0])
            if elongation > 5.0:
                cy_r, cx_r = coords.mean(axis=0)
                angle = 0.5 * np.degrees(np.arctan2(2 * cov[0, 1], cov[0, 0] - cov[1, 1]))
                major, minor = 2 * np.sqrt(eigs[1]), 2 * np.sqrt(eigs[0])
                errors["trailing"].append({"x": float(cx_r), "y": float(cy_r), "w": float(major), "h": float(minor), "angle": float(angle)})
    except Exception as e:
        print(f"[TRAIL DETECT] {e}")

    try:
        sat_thresh = np.percentile(norm, 99.9)
        sat_mask = (norm >= sat_thresh).astype(np.uint8)
        sat_mask = ndimage.binary_dilation(sat_mask, iterations=2).astype(np.uint8)
        errors["saturation_area_fraction"] = float(sat_mask.sum()) / float(h * w)
        labeled, n_obj = ndimage.label(sat_mask)
        for obj_id in range(1, n_obj + 1):
            region = labeled == obj_id
            area = region.sum()
            if area < 200:
                continue
            coords = np.argwhere(region)
            cy_r, cx_r = coords.mean(axis=0)
            ry = (coords[:, 0].max() - coords[:, 0].min()) / 2
            rx = (coords[:, 1].max() - coords[:, 1].min()) / 2
            errors["saturation"].append({"x": float(cx_r), "y": float(cy_r), "r": float(np.hypot(rx, ry))})
    except Exception as e:
        print(f"[SAT DETECT] {e}")

    try:
        bg_size = max(15, min(h, w) // 12)
        local_bg = ndimage.uniform_filter(norm, size=bg_size)
        residual = norm - local_bg
        med_resid = float(np.nanmedian(residual))
        mad = float(np.nanmedian(np.abs(residual - med_resid)))
        sigma_est = max(1.4826 * mad, 1e-3)
        dark_mask = (residual <= (med_resid - 5 * sigma_est)).astype(np.uint8)
        dark_mask = ndimage.binary_opening(dark_mask, iterations=1).astype(np.uint8)
        labeled, n_obj = ndimage.label(dark_mask)
        for obj_id in range(1, n_obj + 1):
            region = labeled == obj_id
            area = region.sum()
            if area < 15 or area > 0.02 * h * w:
                continue
            coords = np.argwhere(region)
            cy_r, cx_r = coords.mean(axis=0)
            ry = (coords[:, 0].max() - coords[:, 0].min()) / 2
            rx = (coords[:, 1].max() - coords[:, 1].min()) / 2
            if rx < 1e-6 or ry < 1e-6:
                continue
            errors["dust"].append({"x": float(cx_r), "y": float(cy_r), "r": float(np.hypot(rx, ry))})
    except Exception as e:
        print(f"[DUST DETECT] {e}")

    try:
        raw_min, raw_max = np.nanmin(data), np.nanmax(data)
        tol = max(1.0, 0.02 * (raw_max - raw_min))
        dead_mask = (data <= raw_min + tol).astype(np.uint8)
        labeled_dead, n_dead_obj = ndimage.label(dead_mask)
        largest_dead_area = 0
        if n_dead_obj > 0:
            sizes = ndimage.sum(dead_mask, labeled_dead, index=range(1, n_dead_obj + 1))
            if len(sizes):
                largest_dead_area = int(np.max(sizes))
        errors["dead_zone_fraction"] = float(largest_dead_area) / float(h * w)
    except Exception as e:
        print(f"[DEAD ZONE DETECT] {e}")

    try:
        border = max(20, int(min(h, w) * 0.05))
        center_med = np.nanmedian(norm[border:h-border, border:w-border])
        center_std = np.nanstd(norm[border:h-border, border:w-border])
        edge_strips = [norm[:border, :], norm[h-border:, :], norm[:, :border], norm[:, w-border:]]
        for strip in edge_strips:
            if abs(np.nanmedian(strip) - center_med) > 3 * center_std:
                errors["edge"] = True
                break
    except Exception as e:
        print(f"[EDGE DETECT] {e}")

    return errors


def classify_plate_quality(errors, n_sources, n_matched=None, lim_mag_apass=None, lim_mag_atlas=None,
                            lim_mag_median_apass=None, lim_mag_median_atlas=None):
    """Defect + limiting-magnitude quality gate for 02_B's purposes (is
    this plate clean enough to trust its detected sources as training
    data). Deliberately simpler than 02_A's APASS-match-fraction version
    -- 02_B doesn't run full field photometry, so there's no match
    fraction to fold in here."""
    n_scratches = len(errors["scratches"])
    n_trailing = len(errors["trailing"])
    n_saturation = len(errors["saturation"])
    n_dust = len(errors["dust"])
    n_edge = 1 if errors["edge"] else 0
    dead_zone_fraction = errors.get("dead_zone_fraction", 0.0)
    saturation_area_fraction = errors.get("saturation_area_fraction", 0.0)
    total_defects = n_scratches + n_trailing + n_saturation + n_dust + n_edge

    if dead_zone_fraction >= 0.08:
        return "defective"
    if saturation_area_fraction >= 0.05:
        return "defective"
    if (n_saturation >= 8) or (errors["edge"] and total_defects >= 6) or (n_sources < 3 and total_defects >= 3):
        return "defective"
    if total_defects >= 8:
        return "too_many_errors"

    defect_based_good = total_defects < 8
    if n_matched is not None:
        match_fraction = (n_matched / n_sources) if n_sources > 0 else 0.0
        defect_based_good = defect_based_good and (match_fraction >= 0.5)

    lim_mag = lim_mag_apass if lim_mag_apass is not None else lim_mag_atlas
    lim_mag_reference = lim_mag_median_apass if lim_mag_apass is not None else lim_mag_median_atlas
    lim_mag_is_deep = (lim_mag is None) or (lim_mag_reference is None) or (lim_mag >= lim_mag_reference)

    if defect_based_good and lim_mag_is_deep:
        return "good_match"
    return "fair"


def annotate_errors(ax, errors):
    """Verbatim from 02_A -- draws error overlays + legend on an axis."""
    legend_handles = []
    for s in errors["scratches"]:
        ax.plot([s["x0"], s["x1"]], [s["y0"], s["y1"]], color="magenta", linewidth=1.2, alpha=0.8, linestyle="--")
    if errors["scratches"]:
        legend_handles.append(Line2D([0], [0], color="magenta", linestyle="--", label=f'Scratch x{len(errors["scratches"])}'))
    for t in errors["trailing"]:
        ax.add_patch(mpatches.Ellipse((t["x"], t["y"]), width=t["w"], height=t["h"], angle=t["angle"],
                                       edgecolor="orange", facecolor="none", linewidth=1.5, alpha=0.85))
    if errors["trailing"]:
        legend_handles.append(mpatches.Patch(edgecolor="orange", facecolor="none", label=f'Trailing x{len(errors["trailing"])}'))
    for s in errors["saturation"]:
        ax.add_patch(plt.Circle((s["x"], s["y"]), s["r"], edgecolor="red", facecolor="none", linewidth=1.5, alpha=0.8, linestyle="-."))
    if errors["saturation"]:
        legend_handles.append(mpatches.Patch(edgecolor="red", facecolor="none", label=f'Saturation x{len(errors["saturation"])}'))
    for d in errors["dust"]:
        ax.add_patch(plt.Circle((d["x"], d["y"]), d["r"], edgecolor="deepskyblue", facecolor="none", linewidth=1.2, alpha=0.8, linestyle=":"))
    if errors["dust"]:
        legend_handles.append(mpatches.Patch(edgecolor="deepskyblue", facecolor="none", label=f'Dust x{len(errors["dust"])}'))
    if legend_handles:
        ax.legend(handles=legend_handles, loc="upper left", fontsize=6, framealpha=0.6, facecolor="black", labelcolor="white", edgecolor="gray")


def categorize_plate(meta):
    quality = meta.get("quality", "fair")
    if meta.get("human_verdict") == "approved" and quality == "fair":
        quality = "good_match"
    is_good = quality in ("good_match", "fair")
    return "ideal" if quality == "good_match" else ("good_no_target" if is_good else "defective_no_target")


def fair_plate_features(errs, n_sources, lim_mag_apass, lim_mag_atlas, lim_mag_median_apass, lim_mag_median_atlas):
    lim = lim_mag_apass if lim_mag_apass is not None else lim_mag_atlas
    ref = lim_mag_median_apass if lim_mag_apass is not None else lim_mag_median_atlas
    lim_z = float(lim - ref) if (lim is not None and ref is not None) else 0.0
    return np.array([
        len(errs.get("scratches", [])), len(errs.get("trailing", [])),
        len(errs.get("saturation", [])), len(errs.get("dust", [])),
        errs.get("dead_zone_fraction", 0.0) * 10.0,
        errs.get("saturation_area_fraction", 0.0) * 10.0,
        lim_z, float(np.log1p(n_sources)),
    ], dtype=float)


# ---------------------------------------------------------------------
# VizieR TSV helpers -- shared by Gaia and APASS queries below.
# ---------------------------------------------------------------------

def _vizier_tsv_query(source, ra, dec, radius_arcsec, out_columns, max_rows=2000, timeout=30):
    params = {
        "-source": source,
        "-c": f"{ra} {dec}",
        "-c.rs": radius_arcsec,
        "-out": out_columns,
        "-out.max": max_rows,
    }
    url = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return response.read().decode("utf-8")

def _parse_tsv_rows(raw_text, n_cols, converters):
    rows = []
    for line in raw_text.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("-"):
            continue
        parts = line.split("\t")
        if len(parts) < n_cols:
            continue
        try:
            row = tuple(conv(p) for conv, p in zip(converters, parts[:n_cols]))
        except ValueError:
            continue
        rows.append(row)
    if rows and not isinstance(rows[0][0], (int, float)):
        rows = rows[1:]
    return rows

def _float_or_nan(s):
    s = s.strip()
    return float(s) if s else np.nan


# ---------------------------------------------------------------------
# Identity crossmatch (Gaia / SIMBAD / paradigm)
# ---------------------------------------------------------------------

def get_plate_catalog_gaia(wcs, shape, gaia_cache, radius_deg=GAIA_QUERY_RADIUS_DEG):
    """Gaia DR3 stars near the plate center via VizieR TSV (I/355/gaiadr3).
    Cached rows are (ra, dec, source_id) tuples."""
    cx, cy = shape[1] / 2, shape[0] / 2
    ra_c, dec_c = wcs.pixel_to_world_values(cx, cy)
    ra_c, dec_c = float(np.array(ra_c)), float(np.array(dec_c))
    key = (round(ra_c, 4), round(dec_c, 4))
    if key in gaia_cache:
        return gaia_cache[key]
    radius_arcsec = radius_deg * 3600.0
    try:
        raw = _vizier_tsv_query("I/355/gaiadr3", ra_c, dec_c, radius_arcsec, "RA_ICRS,DE_ICRS,Source")
        rows = _parse_tsv_rows(raw, 3, [float, float, int])
        gaia_cache[key] = rows
        return rows
    except Exception as e:
        print("[GAIA ERROR]", e)
        gaia_cache[key] = []
        return []


def get_plate_catalog_gaia_point(ra, dec, radius_arcsec=NAME_MATCH_SEP_ARCSEC):
    try:
        raw = _vizier_tsv_query("I/355/gaiadr3", ra, dec, radius_arcsec, "Source", max_rows=1)
        rows = _parse_tsv_rows(raw, 1, [int])
        return rows[0][0] if rows else None
    except Exception as e:
        print("[GAIA POINT ERROR]", e)
        return None


def get_plate_catalog_simbad(wcs, shape, simbad_cache, radius_deg=0.15):
    cx, cy = shape[1] / 2, shape[0] / 2
    ra_c, dec_c = wcs.pixel_to_world_values(cx, cy)
    ra_c, dec_c = float(np.array(ra_c)), float(np.array(dec_c))
    key = (round(ra_c, 4), round(dec_c, 4))
    if key in simbad_cache:
        return simbad_cache[key]
    center = SkyCoord(ra_c * u.deg, dec_c * u.deg)
    try:
        simbad = Simbad()
        simbad.TIMEOUT = 60
        simbad.add_votable_fields("main_id", "ra", "dec")
        result = simbad.query_region(center, radius=radius_deg * u.deg)
        if result is None:
            result = []
        simbad_cache[key] = result
        return result
    except Exception as e:
        print("[SIMBAD ERROR]", e)
        simbad_cache[key] = []
        return []


def match_detected_sources_gaia(ra, dec, catalog, max_sep_arcsec=2.5):
    if catalog is None or len(catalog) == 0:
        return ["Unknown"] * len(ra)
    cat_ra = np.array([row[0] for row in catalog])
    cat_dec = np.array([row[1] for row in catalog])
    cat_id = [row[2] for row in catalog]
    cat_coords = SkyCoord(cat_ra * u.deg, cat_dec * u.deg)
    src_coords = SkyCoord(ra * u.deg, dec * u.deg)
    idx, sep, _ = match_coordinates_sky(src_coords, cat_coords)
    labels = []
    for j in range(len(src_coords)):
        if sep[j].arcsec <= max_sep_arcsec:
            labels.append(f"Gaia {cat_id[idx[j]]}")
        else:
            labels.append("Unknown")
    return labels


def match_detected_sources_simbad(ra, dec, catalog, max_sep_arcsec=5):
    if catalog is None or len(catalog) == 0:
        return ["Unknown"] * len(ra)
    try:
        cat_coords = SkyCoord(catalog["RA"], catalog["DEC"], unit=(u.hourangle, u.deg))
    except Exception as e:
        print("[SIMBAD MATCH]", e)
        return ["Unknown"] * len(ra)
    src_coords = SkyCoord(ra * u.deg, dec * u.deg)
    idx, sep, _ = match_coordinates_sky(src_coords, cat_coords)
    labels = []
    for j in range(len(src_coords)):
        if sep[j].arcsec <= max_sep_arcsec:
            name = str(catalog["MAIN_ID"][idx[j]])
            labels.append("Unknown" if name.startswith("Gaia") else name)
        else:
            labels.append("Unknown")
    return labels


def match_against_paradigm(ra, dec, paradigm_marks, max_sep_arcsec=NAME_MATCH_SEP_ARCSEC):
    all_entries = []
    for entries in paradigm_marks.values():
        for e in entries:
            if e.get("label") and e["label"] != "Unknown" and e.get("ra") is not None and e.get("dec") is not None:
                all_entries.append(e)
    if not all_entries:
        return ["Unknown"] * len(ra)
    cat_coords = SkyCoord([e["ra"] for e in all_entries], [e["dec"] for e in all_entries], unit="deg")
    src_coords = SkyCoord(ra * u.deg, dec * u.deg)
    idx, sep, _ = match_coordinates_sky(src_coords, cat_coords)
    labels = []
    for j in range(len(src_coords)):
        if sep[j].arcsec <= max_sep_arcsec:
            labels.append(all_entries[idx[j]]["label"])
        else:
            labels.append("Unknown")
    return labels


def update_name_cache(gaia_labels, simbad_names, gaia_name_cache):
    changed = False
    for g, s in zip(gaia_labels, simbad_names):
        if s == "Unknown" or not g.startswith("Gaia "):
            continue
        try:
            source_id = int(g.replace("Gaia ", ""))
            gaia_name_cache[source_id] = s
            changed = True
        except Exception:
            pass
    return changed


def resolve_names(gaia_labels, simbad_names, gaia_name_cache):
    labels = []
    for g, s in zip(gaia_labels, simbad_names):
        if s != "Unknown":
            labels.append(s)
            continue
        if g.startswith("Gaia "):
            try:
                source_id = int(g.replace("Gaia ", ""))
                if source_id in gaia_name_cache:
                    labels.append(gaia_name_cache[source_id])
                    continue
            except Exception:
                pass
        labels.append("Unknown")
    return labels


def suggest_name_at(ra, dec, max_sep_arcsec=NAME_MATCH_SEP_ARCSEC):
    coord = SkyCoord(ra * u.deg, dec * u.deg)
    try:
        simbad = Simbad()
        simbad.TIMEOUT = 30
        result = simbad.query_region(coord, radius=max_sep_arcsec * u.arcsec)
        if result is not None and len(result) > 0:
            name = str(result[0]["MAIN_ID"])
            if not name.startswith("Gaia"):
                return name, "SIMBAD"
    except Exception as e:
        print("[SIMBAD lookup]", e)
    try:
        source_id = get_plate_catalog_gaia_point(ra, dec, radius_arcsec=max_sep_arcsec)
        if source_id is not None:
            return f"Gaia {source_id}", "Gaia"
    except Exception as e:
        print("[Gaia lookup]", e)
    return "Unknown", "none"


# ---------------------------------------------------------------------
# Lightweight APASS linearity check -- ported from 02_A's field-wide
# APASS query + linearity fit, used only by 02_B's optional "Linearity
# Check" panel. Not run during prescan() or training; this is purely a
# diagnostic on the algorithm's own detected sources.
# ---------------------------------------------------------------------

def get_plate_catalog_apass(wcs, shape, apass_cache, radius_deg=APASS_QUERY_RADIUS_DEG):
    """APASS DR9 stars with a measured B magnitude near the plate center,
    via VizieR TSV (II/336/apass9). Cached rows are (ra, dec, bmag)
    tuples."""
    cx, cy = shape[1] / 2, shape[0] / 2
    ra_c, dec_c = wcs.pixel_to_world_values(cx, cy)
    ra_c, dec_c = float(np.array(ra_c)), float(np.array(dec_c))
    key = (round(ra_c, 4), round(dec_c, 4))
    if key in apass_cache:
        return apass_cache[key]
    radius_arcsec = radius_deg * 3600.0
    try:
        raw = _vizier_tsv_query("II/336/apass9", ra_c, dec_c, radius_arcsec, "RAJ2000,DEJ2000,Bmag,e_Bmag")
        rows = _parse_tsv_rows(raw, 3, [float, float, _float_or_nan])
        rows = [r for r in rows if not np.isnan(r[2])]
        apass_cache[key] = rows
        return rows
    except Exception as e:
        print("[APASS ERROR]", e)
        apass_cache[key] = []
        return []


def match_detected_sources_apass(ra, dec, catalog, max_sep_arcsec=APASS_MATCH_SEP_ARCSEC):
    """Returns an array of APASS B magnitudes (NaN where no match), same
    length and order as ra/dec."""
    n = len(ra)
    if catalog is None or len(catalog) == 0:
        return np.full(n, np.nan)
    cat_ra = np.array([r[0] for r in catalog])
    cat_dec = np.array([r[1] for r in catalog])
    cat_bmag = np.array([r[2] for r in catalog])
    cat_coords = SkyCoord(cat_ra * u.deg, cat_dec * u.deg)
    src_coords = SkyCoord(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg)
    idx, sep, _ = match_coordinates_sky(src_coords, cat_coords)
    result = np.full(n, np.nan)
    for j in range(n):
        if sep[j].arcsec <= max_sep_arcsec:
            result[j] = cat_bmag[idx[j]]
    return result


def compute_linearity_fit(inst_mag, catalog_mag):
    """Instrumental-vs-catalog linearity fit. Returns
    {"slope","intercept","rms","n_used"} or None if too few points."""
    inst_mag = np.asarray(inst_mag)
    catalog_mag = np.asarray(catalog_mag)
    mask = np.isfinite(inst_mag) & np.isfinite(catalog_mag)
    x = catalog_mag[mask]
    y = inst_mag[mask]
    if len(x) < 5:
        return None
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    rms = float(np.sqrt(np.mean(resid ** 2)))
    return {"slope": float(slope), "intercept": float(intercept), "rms": rms, "n_used": int(len(x))}

def compute_fwhm_for_sources(data, xs, ys, edge_radii=None):
    """Gaussian-fit FWHM (pixels) for each source, via a radial profile
    centered on its detected position. No saturation handling here (02_B
    doesn't run full photometry like 02_A does) -- this is a lightweight
    per-plate diagnostic, not a calibrated measurement."""
    if edge_radii is None:
        edge_radii = np.arange(20)
    n = len(xs)
    fwhm = np.full(n, np.nan)
    for i in range(n):
        try:
            rp = RadialProfile(data, (float(xs[i]), float(ys[i])), edge_radii)
            rp.gaussian_fit
            fwhm[i] = float(rp.gaussian_fwhm)
        except Exception:
            pass
    return fwhm