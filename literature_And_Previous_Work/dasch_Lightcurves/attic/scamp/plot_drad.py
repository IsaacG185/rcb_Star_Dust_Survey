#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Plot drad before and after scamp as a function of distance
# if given, plot previous pipeline tnx results
#
# 2010-06-01 Mathieu Servillat


import sys, os
import time
import atpy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

"""
Example input:
outfile_base = 'ac42227_00_01ww'
tab1_name = 'ac42227_00_01ww_match.fits'
tab2_name = 'ac42227_00_01ww_match_scamp.fits'
tab3_name = 'ac42227_00_01ww_match_tnx.fits'
"""

outfile_base = sys.argv[1]

# Open match tables
tab1_name = sys.argv[2]
tab2_name = sys.argv[3]
if len(sys.argv) > 4:
    tab3_name = sys.argv[4]

plot_tnx = False
if os.path.isfile(tab3_name):
    plot_tnx = True

tab1 = atpy.Table(tab1_name, verbose=False)
tab2 = atpy.Table(tab2_name, verbose=False)
if plot_tnx:
    tab3 = atpy.Table(tab3_name, verbose=False)

xmax = max(tab2.data["offset"])
ymax = max(tab2.data["drad"])

# plot and save
# fig = plt.figure()
plt.clf()
plt.title(os.path.basename(outfile_base) + " -- " + time.ctime())
plt.xlabel("Offset [deg]")
plt.ylabel("drad [arcsec]")
plt.plot(tab1.data["offset"], tab1.data["drad"], ".r", ms=1.5)
if plot_tnx:
    scl = xmax / max(tab3.data["offset"])
    plt.plot(tab3.data["offset"] * scl, tab3.data["drad"], ".g", ms=1.5)
plt.plot(tab2.data["offset"], tab2.data["drad"], ".b", ms=1.5)
plt.xlim(0.0, xmax)
plt.ylim(0.0, ymax)
if ymax > 40:
    plt.yticks(np.arange(ymax / 10.0 + 1) * 10)
elif ymax > 10:
    plt.yticks(np.arange(ymax / 2.0 + 1) * 2)
else:
    plt.yticks(np.arange(ymax + 1))
plt.grid()
plt.savefig(outfile_base + "_drad_offset.png")

plt.clf()
fig = plt.figure()
fig.suptitle(os.path.basename(outfile_base) + " -- " + time.ctime())
# ax = Axes3D(fig)
# ax.plot(tab1.data['RA_REF'],tab1.data['Dec_REF'], tab1.data['drad'], '.b', cmap=cm.jet)
ax = fig.add_subplot(121)
img = ax.scatter(
    tab1.data["RA_REF"],
    tab1.data["Dec_REF"],
    c=tab1.data["drad"],
    alpha=0.5,
    edgecolors=None,
)
ax.set_xlabel("RA [deg]")
ax.set_ylabel("Dec [deg]")
ax.set_xlim(min(tab1.data["RA"]), max(tab1.data["RA"]))
ax.set_ylim(min(tab1.data["Dec"]), max(tab1.data["Dec"]))
cbar = fig.colorbar(img, ax=ax, orientation="horizontal")
# cbar.ax.set_xticklabels([str(round(min(tab1.data['drad']))),str(round(max(tab1.data['drad'])))])
# ax.set_zlabel('drad [arcsec]')
ax = fig.add_subplot(122)
img = ax.scatter(
    tab2.data["RA_REF"],
    tab2.data["Dec_REF"],
    c=tab2.data["drad"],
    alpha=0.5,
    edgecolors=None,
)
ax.set_xlabel("RA [deg]")
# ax.set_ylabel('Dec [deg]')
ax.set_xlim(min(tab1.data["RA"]), max(tab1.data["RA"]))
ax.set_ylim(min(tab1.data["Dec"]), max(tab1.data["Dec"]))
cbar = fig.colorbar(img, ax=ax, orientation="horizontal")
# cbar.ax.set_xticklabels([str(round(min(tab2.data['drad']))),str(round(max(tab2.data['drad'])))])
plt.savefig(outfile_base + "_drad_map.png")
