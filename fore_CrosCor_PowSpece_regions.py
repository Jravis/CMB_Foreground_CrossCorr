#############################################################################
#Copyright (c) 2017, Sandeep Rana & Tuhin Ghosh
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#   Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer.
#   Redistributions in binary form must reproduce the above copyright notice,
#      this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#   The name of the author may not be used to endorse or promote products
#      derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
#INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
#BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
#OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
#AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
#WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#POSSIBILITY OF SUCH DAMAGE.
#############################################################################

"""
Code for testing  cross-correlation between different foreground maps
for detailed see "Planck intermediate results. XXII. Frequency dependence of
thermal emission from Galactic dust in intensity and
polarization"
"""


import healpy as hp
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import numpy as np
import sys
sys.path.insert(0, '/home/tolstoy/Documents/software/PolSpice_v03-05-01/bin/')
import ispice

region_number = 10
map_Dirname = '/home/tolstoy/Documents/CMB_dust_sync_Cross_Corr/results/regions/region_%d/'%(region_number)


#Polspice routine
#========================================
map_dust   = map_Dirname + 'region_%d_Masked_dust_apod2deg_nside256.fits'%region_number
map_Sync   = map_Dirname + 'region_%d_Masked_Sync_apod2deg_nside256.fits'%region_number
bmask_apod = map_Dirname + 'region_%d_BinaryMask_apod2deg_nside256.fits'%region_number

Nside_map = 256
theta_ap = 2.0

cl_dust = map_Dirname + 'cl_dust_region_%d_%0.1fdeg_apodi.fits' % (region_number, theta_ap)
cl_Sync = map_Dirname + 'cl_Sync_region_%d_%0.1fdeg_apodi.fits' % (region_number, theta_ap)

ispice.ispice(map_dust, cl_dust, nlmax=3*Nside_map-1, weightfile1=bmask_apod,
              thetamax=25.,apodizesigma=1.15*25., apodizetype=1,beam1=60.,label="spice")

ispice.ispice(map_Sync, cl_Sync, nlmax=3*Nside_map-1, weightfile1=bmask_apod,
              thetamax=25.,apodizesigma=1.15*25., apodizetype=1,beam1=60.,label="spice")




