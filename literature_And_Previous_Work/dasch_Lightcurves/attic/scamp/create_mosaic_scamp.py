#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# create new mosaic and update header info
#
# 2009-12-10 Mathieu Servillat

import sys
import pyfits
import pywcs
import atpy
import shutil, os

img_name = '/media/raid005/ExposureData/Plates/ac/42227_00/ac42227_00_01ww.fit'
hd_scamp_name = '/home/mservillat/Pipeline/scamp/test_ac42227_00_01ww/ac42227_00_01ww.head'

img_name = sys.argv[1]
hd_scamp_name = sys.argv[2]
img_out = img_name.split('.fit')[0]+'_scamp.fit'
hdr_out = hd_scamp_name.split('.cat.head')[0]+'_scamp.hdr'

print '----- Apply Scamp Solution -----'
print 'Image:',img_name
print 'Scamp solution:',hd_scamp_name

# read scamp solution in .head file
print '\n-- Read Scamp solution'
# has to remove the END line of the file
fin = open(hd_scamp_name,"r" )
lines = fin.readlines()
fin.close()
end='END     \n'
if end in lines: lines.remove(end)
fout = open(hd_scamp_name+'_cor',"w" )
fout.writelines(lines)
fout.close()
# now read
hd_scamp = pyfits.Header()
hd_scamp.fromTxtFile(hd_scamp_name+'_cor')
#print hd_scamp.ascard

#copy mosaic
print '\n-- Copy to new mosaic'
if not os.path.isfile(img_out):
    print img_name,'-->',img_out
    shutil.copy(img_name,img_out)
else:
    print 'Use existing file'

# update img header
print '\n-- Update mosaic header'
hdu = pyfits.open(img_out,mode='update')
hd = hdu[0].header

for kw in ['LTV1', 'LTV2', 'LTM1_1' ,'LTM2_2' ,'WAT0_001' ,'WAT1_001' ,'WAT2_001' ,'WCSDIM' ,'EQUINOX']:
    del hd[kw]
hd.update('EQUINOX',2000.,comment='Epoch of RA & DEC')
# del previous PV keywords
for kw in hd.keys():
    if 'PV' in kw: del hd[kw]
    if 'WAT' in kw: del hd[kw]
    if 'LTM' in kw: del hd[kw]
    if 'WCSDIM' in kw: del hd[kw]

for kw in hd_scamp.keys():
    if kw != 'HISTORY' and kw != 'COMMENT':
	comm=''
	if ' / ' in str(hd_scamp.ascard[kw]):
	    comm=str(hd_scamp.ascard[kw]).split(' / ')[-1]
	#print kw,'=',hd_scamp[kw],'/',comm
	if kw in hd.keys():
	    print kw,'=',hd[kw],'-->',hd_scamp[kw]
	else:
	    print kw,'=',hd_scamp[kw]
	hd.update(kw,hd_scamp[kw],comment=comm.split(' / ')[-1])

print ''
hd.toTxtFile(hdr_out,clobber=True)

hdu.close()

#hd_init.toTxtFile(hd_scamp_name+'_all',clobber=True)
#hdu=pyfits.PrimaryHDU(header=hd_init)
#hdu.writeto(hd_scamp_name+'_all',clobber=True)
