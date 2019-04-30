#!/usr/bin/env python
## 
# python script to launch the ncl script calculate_mslp_arguments.ncl that extract 
# mean sea level pressure from wrfout files using command line
# arguments.
# This script goes through all NARCliM simulations. If the mslp field already 
# exists then it goes to the next one.
#
# The way the ncl script is called is better explained here: 
#    https://wiki.c2sm.ethz.ch/Wiki/VisNCLBasicExamples
#
# Author: Alejandro Di Luca
# Created: 30/04/2014

import netCDF4 as nc
import numpy as np
import glob
import os
import copy
import sys
import sh
# Loop parameters
domain = 'd01'
simul_names = ['run1']
ncl_cmd = './runncl.sh calculate_mslp_arguments.ncl'
ftype = 'wrfout'
input_dir = '/short/e14/gs9353/now-cordex24/run'
outdir = '/g/data1/e14/gs9353/POST-PROCESS/wrf_cordex24_with-cpl-sst'
year = sys.argv[1]
months = np.arange(1, 12)

for mm, month in enumerate(months):
    indir = '%s/%04i' % (input_dir, year) 
    print('--> INPUT DIRECTORY: %s' % indir)
    sh.mkdir('mkdir -p '+ outdir)
    print('--> OUTPUT DIRECTORY: %s' % outdir)
    file_out = '%s/WRF_mslp_%s_%04i-%02i.nc' % (outdir, domain, year, mmonth)
    cmd = '%s %s %s %04i %02i' % (ncl_cmd, indir, outdir, year, month)       
    print('--> NCL command: %s' % cmd)
    if os.path.isfile(file_out):
        os.system(cmd)
