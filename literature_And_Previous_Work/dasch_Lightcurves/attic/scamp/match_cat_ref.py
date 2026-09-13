#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Match catalog with reference catalog and save match results
#
# 2009-12-10 Mathieu Servillat

import sys
import pyfits
import atpy
import shutil, os
import sourcelists as sl
import numpy as np

"""
Examples of input:
cat_ref_name = 'test_mc39048_01_01r270ww/UCAC-3_0841+1937_r256.cat'
cat_sex_name = 'test_mc39048_01_01r270ww/mc39048_01_01r270ww.cat'
cat_sex_name = 'test_mc39048_01_01r270ww/mc39048_01_01r270ww_scamp.cat'
cat_out_name = 'test_mc39048_01_01r270ww/mc39048_01_01r270ww_match.fit'
cat_out_name = 'test_mc39048_01_01r270ww/mc39048_01_01r270ww_scamp_match.fit'
"""

cat_sex_name = sys.argv[1]
cat_ref_name = sys.argv[2]
cat_out_name = sys.argv[3]
radius = float(sys.argv[4])

cat_sex = atpy.Table(cat_sex_name,type='fits',hdu=2,verbose=False)
print 'Sextractor catalog opened',cat_sex_name
ra1=cat_sex.data['ALPHA_J2000']
dec1=cat_sex.data['DELTA_J2000']
sl.save_ds9reg(cat_sex_name+'.xy.reg',cat_sex.data['XWIN_IMAGE'],cat_sex.data['YWIN_IMAGE'],errors=4.,header='physical')
sl.save_ds9reg(cat_sex_name+'.reg',ra1,dec1,errors=4.)

cat_ref = atpy.Table(cat_ref_name,type='fits',hdu=2,verbose=False)
print 'Reference catalog opened',cat_ref_name
ra2=cat_ref.data['X_WORLD']
dec2=cat_ref.data['Y_WORLD']
sl.save_ds9reg(cat_ref_name+'.reg',ra2,dec2,errors=4.,color='red',tag='ucac3')

ind1,ind2,dist12 = sl.match2(ra1,dec1,ra2,dec2,radius=radius,onlybest=True)
#ind1,ind2,dist12 = sl.match2(ra2,dec2,ra2,dec2,radius=3.,onlybest=True)

dra = (ra1[ind1]-ra2[ind2])*np.cos(np.deg2rad(dec2[ind2]))*3600.
ddec = (dec1[ind1]-dec2[ind2])*3600.
ra_center = np.mean(ra1)
dec_center = np.mean(dec1)
offset = sl.separation(ra1[ind1],dec1[ind1],ra_center,dec_center)

cat_out = atpy.Table()
cat_out.add_column('RA',ra1[ind1])
cat_out.add_column('Dec',dec1[ind1])
cat_out.add_column('RA_REF',ra2[ind2])
cat_out.add_column('Dec_REF',dec2[ind2])
cat_out.add_column('dra',dra)
cat_out.add_column('ddec',ddec)
cat_out.add_column('drad',dist12)
cat_out.add_column('offset',offset)

cat_out.write(cat_out_name,overwrite=True,verbose=False)
print 'Match catalog saved',cat_out_name
