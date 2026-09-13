#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Get num of detection in Sextractor table 
#
# 2010-06-01 Mathieu Servillat

import sys
import pyfits

"""
fin = 'ac42227_00_01ww.cat.fits'
"""

fin = sys.argv[1]
data = pyfits.getdata(fin,2)
print len(data)