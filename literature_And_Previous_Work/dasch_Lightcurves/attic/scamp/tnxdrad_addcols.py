#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Add columns to tnx_match files generated with
# votable -i #_drad.db -o #_tnx_match.xml
#
# 2010-04-28 Mathieu Servillat
#
# From Ed:
# Set the x axis to  sqrt(1.56*(ix-25)*(ix-25)+(iy-25)*(iy-25)) for mc plates (North along narrow axis) or  sqrt(1.56*(ix-25)*(ix-25)+(iy-25)*(iy-25)) (North along wide axis)
# Set the y axis to sqrt(draMedian*draMedian+ddecMedian*ddecMedian)

import sys
import atpy
import numpy as np

catin = 'test_mc39048_01_01r270ww/mc39048_01_01r270ww_tnx_match.xml'
catin = sys.argv[1]
catout = catin.split('.')[0]+'.fits'

print 'Add offset and drad to file:'
print 'In:',catin
print 'Out:',catout

cat = atpy.Table(catin,verbose=False)

#sqrt(1.56*(ix-25)*(ix-25)+(iy-25)*(iy-25))
xaxis = np.sqrt( 1.56*(cat.data['ix']-25)**2 + (cat.data['iy']-25)**2 )
yaxis = np.sqrt(cat.data['draMedian']**2+cat.data['ddecMedian']**2)

cat.add_column('offset',xaxis)
cat.add_column('drad',yaxis)

cat.write(catout,overwrite=True,verbose=False)
