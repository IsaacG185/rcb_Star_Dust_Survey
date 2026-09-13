#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# Converts ucac3 table to FITS_LDAC table
#
# 2010-06-01 Mathieu Servillat

import sys
import asciitable
import atpy
import numpy as np

fin = sys.argv[1]
fout = sys.argv[2]

# read table from imcat -dh -c ucac3
data = asciitable.read(fin,delimiter=' ',header_start=5,data_start=6)
nrows = len(data)

# create fits_ldac
out  = atpy.TableSet()
# create imhead table
out1 = atpy.Table()
out1.table_name = 'LDAC_IMHEAD'
out1.add_column('Field Header Card',['SIMPLE  =                    T / This is a FITS file'],dtype='|S1680')
out.append(out1)
# create objects table
out2 = atpy.Table()
out2.table_name = 'LDAC_OBJECTS'
out2.add_column('X_WORLD',data['RA2000'],unit='deg')
out2.add_column('Y_WORLD',data['Dec2000'],unit='deg')
out2.add_column('ERRA_WORLD',(np.zeros(nrows)+0.01)/3600.,dtype='float32',unit='deg')
out2.add_column('ERRB_WORLD',(np.zeros(nrows)+0.01)/3600.,dtype='float32',unit='deg')
out2.add_column('MAG',data['MagM'],dtype='float32',unit='mag')
out2.add_column('MAGERR',np.zeros(nrows)+0.1,dtype='float32',unit='mag')
out.append(out2)

# save new file, ready for scamp
out.write(fout,overwrite=True,verbose=False)
print len(out2),'reference sources in',fout
