#!/usr/bin/env python3
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Plot drad before and after scamp as a function of distance
# if given, plot previous pipeline tnx results

import sys, os
import time

from astropy.table import Table
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

if len(sys.argv) < 5:
    print("Plot drad map and drap vs offset")
    print("Usage: plot_drad.py tab1_name tab2_name tab3_name scale")
    print(
        "  tab1_name: fits table with drad and offset columns (eg before scamp correction)"
    )
    print(
        "  tab2_name: fits table with drad and offset columns (eg after scamp correction)"
    )
    print("  tab3_name: fits table with drad and offset columns (eg from pipeline)")
    print("  scale: size of 1 pixel to overplot 3 pixels limit")
    sys.exit()

outfile_base = sys.argv[1]

# Open match tables
tab1_name = sys.argv[2]
tab2_name = sys.argv[3]
tab3_name = sys.argv[4]
scale = float(sys.argv[5])
if len(sys.argv) > 6:
    testflag = sys.argv[6]
else:
    testflag = ""

plot_initial = False
plot_scamp = False
plot_tnx = False
if os.path.isfile(tab1_name):
    plot_initial = True
if os.path.isfile(tab2_name):
    plot_scamp = True
if os.path.isfile(tab3_name):
    plot_tnx = True

if plot_initial:
    tab1 = Table.read(tab1_name)
if plot_scamp:
    tab2 = Table.read(tab2_name)
if plot_tnx:
    tab3 = Table.read(tab3_name)

if plot_scamp:
    xmax = max(tab2["offset"])
    ymax = max(tab2["drad"])
else:
    if plot_initial:
        xmax = max(tab1["offset"])
        ymax = max(tab1["drad"])
    else:
        print("ERROR: No matches found in " + tab1_name + " or " + tab2_name)
        exit(1)

# plot drad vs offset
plt.clf()
fig = plt.figure()
fig.suptitle(time.ctime() + "\n" + os.path.basename(outfile_base))

if ymax > 260:
    yt = np.arange(ymax / 20.0 + 1) * 20
elif ymax > 40:
    yt = np.arange(ymax / 10.0 + 1) * 10
elif ymax > 20:
    yt = np.arange(ymax / 2.0 + 1) * 2
else:
    yt = np.arange(ymax + 1)

ax = fig.add_subplot(121, position=[0.27, 0.1, 0.65, 0.8])
p_ms = 2.5
p_alpha = 0.5

if plot_initial:
    p1 = ax.plot(
        tab1["offset"],
        tab1["drad"],
        ".r",
        ms=p_ms,
        alpha=p_alpha,
        label="Initial mosaic",
    )
p2 = ax.plot(
    tab2["offset"],
    tab2["drad"],
    ".b",
    ms=p_ms,
    alpha=p_alpha,
    label="Scamp correction",
)

if plot_tnx:
    p3 = ax.plot(
        tab3["offset"] * xmax / max(tab3["offset"]),
        tab3["drad"],
        ".g",
        ms=p_ms,
        alpha=p_alpha,
        label="Pipeline correction",
    )

plim = ax.plot([0, xmax], [scale * 3, scale * 3], "--k", lw=1, label="3 pixels limit")
ax.set_xlabel("Offset [deg]")
ax.set_xlim(0.0, xmax)
ax.set_ylim(0.0, ymax)
ax.set_yticks(yt)
ax.set_yticklabels([])
ax.grid()
font = font_manager.FontProperties(size="x-small")
ax.legend(loc=2, numpoints=1, markerscale=5, prop=font)

ax = fig.add_subplot(122, position=[0.12, 0.1, 0.15, 0.8])

if plot_initial:
    h1 = ax.hist(
        tab1["drad"],
        bins=80,
        density=True,
        histtype="step",
        align="mid",
        orientation="horizontal",
        color="r",
        lw=1.5,
    )

n, bins, patches = ax.hist(
    tab2["drad"],
    bins=80,
    density=True,
    histtype="step",
    align="mid",
    orientation="horizontal",
    color="b",
    lw=1.5,
)

if plot_tnx:
    sel = np.where(tab3["drad"] > 0)
    h3 = ax.hist(
        tab3["drad"][sel],
        bins=80,
        density=True,
        histtype="step",
        align="mid",
        orientation="horizontal",
        color="g",
        lw=1.5,
    )

plim = ax.plot(
    ax.get_xlim(), [scale * 3, scale * 3], "--k", lw=1, label="3 pixels limit"
)

ax.set_ylabel("drad [arcsec]")
ax.set_ylim(0.0, ymax)
ax.set_xlim(ax.get_xlim()[1], 0.0)
ax.set_yticks(yt)
ax.set_xticklabels([])
ax.grid()

plt.savefig(outfile_base + testflag + "_drad_offset.png")

plt.clf()
fig = plt.figure()
fig.suptitle(time.ctime() + "\n" + os.path.basename(outfile_base))

if plot_initial:
    xlm = [max(tab1["RA"]), min(tab1["RA"])]
    ylm = [min(tab1["Dec"]), max(tab1["Dec"])]
else:
    xlm = [max(tab2["RA"]), min(tab2["RA"])]
    ylm = [min(tab2["Dec"]), max(tab2["Dec"])]

asp = 1.0 / np.cos(np.mean(ylm) * np.pi / 180.0)

ax = fig.add_subplot(121, aspect=asp, position=[0.13, 0.08, 0.35, 0.8])
if plot_initial:
    img = ax.scatter(
        tab1["RA_REF"],
        tab1["Dec_REF"],
        c=tab1["drad"],
        alpha=0.5,
        edgecolors=None,
    )
    ax.set_xlabel("RA [deg]")
    ax.set_ylabel("Dec [deg]")
    ax.set_xlim(xlm)
    ax.set_ylim(ylm)
    cbar = fig.colorbar(img, ax=ax, orientation="horizontal")

    for t in cbar.ax.get_xticklabels():
        t.set_fontsize(8)

ax = fig.add_subplot(122, aspect=asp, position=[0.58, 0.08, 0.35, 0.8])
img = ax.scatter(
    tab2["RA_REF"],
    tab2["Dec_REF"],
    c=tab2["drad"],
    alpha=0.5,
    edgecolors=None,
)

ax.set_xlabel("RA [deg]")
# ax.set_ylabel('Dec [deg]')
ax.set_xlim(xlm)
ax.set_ylim(ylm)
cbar = fig.colorbar(img, ax=ax, orientation="horizontal")
for t in cbar.ax.get_xticklabels():
    t.set_fontsize(8)

plt.savefig(outfile_base + testflag + "_drad_map.png")

# 2010-06-01 Mathieu Servillat
# 2011-01-18 Edward J. Los - Reduce markerscale for Fedora 14.
# 2013-03-01 Edward J. Los - Add a test qualifier
