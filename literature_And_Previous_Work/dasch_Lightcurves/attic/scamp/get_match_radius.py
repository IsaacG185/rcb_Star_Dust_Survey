#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# return match_radius according to
#
# 2010-06-01 Mathieu Servillat

import sys
import string

mosaic = sys.argv[1]
mosaic_base = mosaic.split('_')[0]
exclude = set(string.digits)
series = ''.join(ch for ch in mosaic_base if ch not in exclude)
exclude = set(string.letters)
num = int(''.join(ch for ch in mosaic_base if ch not in exclude))
match_radius = 10

if series == 'ac':
    if num >= 32289: 
	match_radius = 210
    else:
	match_radius = 500
if series == 'mc':
    match_radius = 10

print match_radius
