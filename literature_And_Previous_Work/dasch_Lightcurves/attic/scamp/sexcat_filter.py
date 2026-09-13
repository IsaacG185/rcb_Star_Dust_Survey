#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Filter Sextractor table by limiting the number of detection 
# to the N brightest objects
#
# 2010-06-01 Mathieu Servillat

import sys,os
import atpy
import numpy as np

"""
fin = 'ac42227_00_01ww.cat.fits'
nkeep = 10000
"""

fin = sys.argv[1]
nkeep = int(sys.argv[2])

# save
#os.copy(fin,os.path.basename(fin)+'_save')

# read fits_ldac table
tset = atpy.TableSet(fin,verbose=False)
tobj = tset.tables[1]
ndet = len(tobj)
if ndet > nkeep:
    order = np.argsort(tobj.data['MAG_ISO'])
    tobj = tobj.rows(order)
    tobj = tobj.rows(range(nkeep))
    del tset.tables[1]
    tset.append(tobj)
    # save new file, ready for scamp
    tset.write(fin,overwrite=True,verbose=False)
print len(tobj),'/',ndet,'sources kept in',fin
