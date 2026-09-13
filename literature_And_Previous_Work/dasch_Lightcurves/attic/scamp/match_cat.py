#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Match 2 catalogs and save match results
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

cat_ref = atpy.Table(cat_ref_name,type='fits',hdu=2,verbose=False)
print 'Reference catalog opened',cat_ref_name
ra2=cat_ref.data['ALPHA_J2000']
dec2=cat_ref.data['DELTA_J2000']

ind1,ind2,dist12 = sl.match2(ra1,dec1,ra2,dec2,radius=radius,onlybest=True)
#ind1,ind2,dist12 = sl.match2(ra2,dec2,ra2,dec2,radius=3.,onlybest=True)

sl.save_ds9reg_line(cat_ref_name+'.line.reg',ra1[ind1],dec1[ind1],ra2[ind2],dec2[ind2],color='green',tag='lines')

cat_out = atpy.Table()
cat_out.add_column('RA',ra1[ind1])
cat_out.add_column('Dec',dec1[ind1])
cat_out.add_column('RA_REF',ra2[ind2])
cat_out.add_column('Dec_REF',dec2[ind2])
cat_out.add_column('Dist',dist12)

cat_out.write(cat_out_name,overwrite=True,verbose=False)
print 'Match catalog saved',cat_out_name
