"""
In this I am computing cl on fly without combining regions
"""
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

import numpy as np
import healpy as hp
#from  multiprocessing import Process
#from scipy import special
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import scipy.stats
from tqdm import *
import pymaster as nmt

mpl.rcParams['axes.linewidth'] = 1.5
colombi1_cmap = ListedColormap(np.loadtxt("Planck_Parchment_RGB.txt")/255.)
colombi1_cmap.set_bad("gray")
colombi1_cmap.set_under("white")
cmap1 = colombi1_cmap


# Global values

Nside_map = 256
Nside_ref = 8
#query_radius = hp.nside2resol(Nside_ref, arcmin=False)
query_radius = np.deg2rad(10.0)
plt.style.use("classic")


delta_ell = 12
b = nmt.NmtBin(Nside_map, nlb=delta_ell)
ell_arr_bin = b.get_effective_ells()

beam = hp.gauss_beam(np.radians(60.0/60.0), lmax=3*Nside_map-1)
pix_win = hp.pixwin(Nside_map)
beam_binned = b.bin_cell([beam])[0]
pix_win_binned = b.bin_cell([pix_win[0:3*256]])[0]


#========================================================

def cross_corr(map1, map2):

    RCoeff, temp = scipy.stats.pearsonr(map1, map2)
    return RCoeff

#===============================================================================

def compute_cl(mask_raw, map1):

    aposcale=2.0
    mask_C1=nmt.mask_apodization(mask_raw,aposcale, apotype="C1")
    ell_bin = b.get_effective_ells()
    f_TT=0.0
    f_TT = nmt.NmtField(mask_C1, [map1])

    f_sky = (np.sum(mask_C1)*1.0)/hp.nside2npix(Nside_map)
    nu_b = (2*ell_arr_bin + 1) * (delta_ell*1.0) *f_sky

    cl_tmp = nmt.compute_full_master(f_TT, f_TT, b)

    cl_tt = np.asarray(cl_tmp[0]/(beam_binned**2 * pix_win_binned**2))
    Cosmic_Var = cl_tt * np.sqrt(2./nu_b)

    return cl_tt, Cosmic_Var, mask_C1


#===============================================================================



def main(seq_name, masking):

    """ 
    fits_filename = "../maps/COM_CompMap_%s-GNILC-F857_2048_R2.00.fits" % seq_name[0]
    map_1 = hp.read_map(fits_filename)
    fits_filename = "../maps/haslam408_dsds_Remazeilles2014.fits" # Actual Haslam map Nside 512
    map_2 = hp.read_map(fits_filename)
    
    map_1 = hp.sphtfunc.smoothing(map_1, fwhm=np.radians(56./60.)) # smoothing Dust map to 1 degree
    map_2 = hp.sphtfunc.smoothing(map_2, fwhm=np.radians(56./60.))

    map_1 = hp.ud_grade(map_1, 512)  # degrading GNILC Dust map to Nside 256
    map_2 = hp.ud_grade(map_2, 512)  # degrading Haslam map to Nside 256

    hp.write_map("../maps/GNILC857_ns512_56arcmSmooth_Dust.fits", map_1)
    hp.write_map("../maps/haslam408_ns512_56arcmSmooth_dsds_Remazeilles2014.fits", map_2)

    """


    map_1 = hp.read_map("../maps/GNILC857_ns256_60arcmSmooth_Dust.fits")
    map_2  =hp.read_map("../maps/haslam408_ns256_60arcmSmooth_dsds_Remazeilles2014.fits")

    pixel_indx_map = np.arange(hp.nside2npix(Nside_map), dtype=np.int32) # Nside_ref = 32 pixel index array
    pixel_indx_ref = np.arange(hp.nside2npix(Nside_ref), dtype=np.int32) # Nside_ref = 32 pixel index array

    theta, phi =hp.pixelfunc.pix2ang(Nside_ref, pixel_indx_ref, nest=False, lonlat=False)
    ipixl = hp.pixelfunc.ang2pix(Nside_map, theta, phi, nest=False, lonlat=False)
 

    min_lat = 70.0
    max_lat = 110.0
    lat, lon = hp.pix2ang(Nside_map, pixel_indx_map, lonlat=False)
    gal_mask_map = (np.rad2deg(lat) > min_lat) * (np.rad2deg(lat) < max_lat)
   
    gal_mask_index = pixel_indx_map[gal_mask_map]
    mask  = np.zeros(hp.nside2npix(Nside_map))        # Correlation Coefficient array for actual masking

    selected_pixel = []
    for ipix in tqdm(iterable=ipixl):
        x, y, z =  hp.pix2vec(Nside_map, ipix)
        pix_indx_arr = hp.query_disc(Nside_map, (x, y, z), query_radius, nest=False)
        indx_arr = np.setdiff1d(pix_indx_arr, gal_mask_index)
        
        T1 = map_1[indx_arr]
        T2 = map_2[indx_arr]
        
        if T1.size > 0: 
            coff = cross_corr(T1, T2)
        if coff >= 0.8:
            mask[indx_arr] = 1.0
            selected_pixel.append(ipix)

            cl_dust, err_dust, mask_apod = compute_cl(mask, map_1)
            cl_sync, err_sync, mask_apod  = compute_cl(mask, map_2)

            fname1 = '/jbodstorage/data_sandeep/sandeep/cross-corr_data/nu_GE_0.8_ns256_dns8_qrad_10deg/'\
                     'PseudoCl_nu_GE_0.8_ns256_dns8_qrad_10deg_delta_ell%d_pixindx_%d.txt'% (delta_ell, ipix)
            np.savetxt(fname1, zip(ell_arr_bin, cl_dust, cl_sync, err_dust, err_sync), 
                            delimiter='\t', header="ell_bin\tCl_dust\tCl_sync\tCerr_dust\tCerr_sync")

            hp.write_map('/jbodstorage/data_sandeep/sandeep/cross-corr_data/nu_GE_0.8_ns256_dns8_qrad_10deg/'\
                        'bmask_nu_GE_0.8_ns256_dns8_qrad_10deg_delta_ell%d_pixindx_%d.fits'% (delta_ell, ipix),
                         mask,coord='G',overwrite='True')

        mask  = np.zeros(hp.nside2npix(Nside_map))  

    selected_pixel = np.asarray(selected_pixel)
    np.savetxt('/jbodstorage/data_sandeep/sandeep/cross-corr_data/nu_GE_0.8_ns256_dns8_qrad_10deg/'\
               'PseudoCl_nu_GE_0.8_ns256_dns8_qrad_10deg_pixel_index.txt', selected_pixel, fmt='%d')

#Initialize binning scheme with 4 ells per bandpower
#beam = hp.gauss_beam(np.radians(60.0/60.0), lmax=3*Nside-1)
#pix_win = hp.pixwin(Nside)
#beam_binned = b.bin_cell([beam])[0]
#pix_win_binned = b.bin_cell([pix_win[0:3*256]])[0]
#f_sky = (np.sum(rho_mask)*1.0)/hp.nside2npix(Nside)
#nu_b = (2*ell_arr_bin + 1) * (delta_ell*1.0) *f_sky
#cl_dust = np.asarray(cl_dust[0]/(beam_binned**2 * pix_win_binned**2))

if __name__ == "__main__":

    name1 = 'Dust'
    name2 = 'Synchrotron'
    seq = [name1, name2]
    main(seq, True)














