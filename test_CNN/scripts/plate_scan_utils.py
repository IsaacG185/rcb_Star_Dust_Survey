
from pathlib import Path
from astropy.io import fits as astrofits
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u
from photutils.detection import DAOStarFinder, find_peaks
from photutils.aperture import CircularAperture, aperture_photometry
from astroquery.gaia import Gaia
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

def detect_sources(fits_path, aperture_radius=APERTURE_RADIUS):
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
        for strip in [norm[:border, :], norm[h-border:, :], norm[:, :border], norm[:, w-border:]]:
            if abs(np.nanmedian(strip) - center_med) > 3 * center_std:
                errors["edge"] = True
                break
    except Exception as e:
        print(f"[EDGE DETECT] {e}")
    return errors

def annotate_errors(ax, errors):
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
