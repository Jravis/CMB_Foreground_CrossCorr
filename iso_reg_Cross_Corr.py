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


mpl.rcParams['axes.linewidth'] = 3.0
colombi1_cmap = ListedColormap(np.loadtxt("../Planck_Parchment_RGB.txt")/255.)
colombi1_cmap.set_bad("gray")
colombi1_cmap.set_under("white")
cmap1 = colombi1_cmap


# Global values

Nside_map = 256
Nside_ref = 16


#print "Enter which query radius you want (degrees)"
#q_rad = float(raw_input(""))
#print q_rad
query_radius = np.deg2rad(hp.nside2resol(Nside_ref, arcmin=True)/60.0)
q_rad = query_radius
Npix_map = hp.nside2npix(Nside_map)
Npix_ref = hp.nside2npix(Nside_ref)

region_number = 1
plot_Dirname= '/home/tolstoy/Documents/CMB_dust_sync_Cross_Corr/results/GNILC-857-Haslam_region_query_rad_3.6deg/'
map_Dirname = '/home/tolstoy/Documents/CMB_dust_sync_Cross_Corr/results/GNILC-857-Haslam_region_query_rad_3.6deg/'

plt.style.use("classic")


#========================================================

def apodize_mask(mask, delta_c):

    npix = len(mask)
    nside = hp.npix2nside(npix)
    ipix = np.arange(npix)
    vec = np.array(hp.pix2vec(nside, ipix))
    apodized_mask = 1.0*mask
    good_pix = np.where(mask == 1)[0]
    bad_pix = mask == 0

    vec = np.array(hp.pix2vec(nside, ipix))
    vec_bad = vec[:, bad_pix]

    for tmp in good_pix:

        vec_good = vec[:, tmp]
        ctheta = vec_bad[0, :]*vec_good[0] + vec_bad[1, :]*vec_good[1] + vec_bad[2, :]*vec_good[2]
        ctheta[ctheta > 1] = 1.0
        ctheta[ctheta < -1] = -1.0

        theta = np.arccos(ctheta)
        delta_i = np.degrees(np.min(theta))
        if delta_i <= delta_c:
            apodized_mask[tmp] = -1.0 / (2*np.pi) * np.sin(2*np.pi*delta_i / delta_c) + delta_i / delta_c

    return apodized_mask


#========================================================

def cross_corr(vec_arr, map1, map2):
    pix_indx_arr = hp.query_disc(Nside_map, vec_arr, query_radius)
    T1 = map1[pix_indx_arr]
    T2 = map2[pix_indx_arr]
    RCoeff, temp = scipy.stats.pearsonr(T1, T2)
    return RCoeff

#=========================================================

def T_T_Corr(pixel_indx, map1, map2):


    x, y, z =  hp.pix2vec(Nside_ref, pixel_indx[0])
    vec_arr1 = np.array([x, y, z])
    pix_indx_arr1 = hp.query_disc(Nside_map, vec_arr1, query_radius)
    lat, lon = hp.vec2ang(vec_arr1, lonlat=False)

    for ipix in xrange(1, len(pixel_indx)):
        x, y, z =  hp.pix2vec(Nside_ref, pixel_indx[ipix])
        vec_arr1 = np.array([x, y, z])

        pix_indx_arr1 = np.concatenate(( pix_indx_arr1, hp.query_disc(Nside_map,
                                        vec_arr1, query_radius)), axis=0)

#===============================================================================

    bmask = np.ones(hp.nside2npix(Nside_map), dtype=np.float64)
    bmask_binary = np.ones(hp.nside2npix(Nside_map), dtype=np.float64)
    bool_arr = np.ones(hp.nside2npix(256), dtype=bool)

    for ind in tqdm(iterable=pix_indx_arr1):
        bool_arr[ind]=False

#    map1[bool_arr] =  hp.UNSEEN
#    map2[bool_arr] =  hp.UNSEEN
#    map1 = np.ma.masked_values(map1, value=-1.6375e+30)
#    map2 = np.ma.masked_values(map2, value=-1.6375e+30)

    bmask_binary[bool_arr] = 0.0
    bmask[bool_arr] = hp.UNSEEN

    fig = plt.figure(7, figsize=(8, 6))
    hp.mollview(bmask_binary, fig=fig.number, xsize=2000, unit='', nest=False, cmap=cmap1, title='Binary mask')
    hp.graticule()

    #plot_name = plot_Dirname+'region_%d_bmask_GNILC-haslam_nside256_queryRadius_%0.1f.pdf'%(region_number, q_rad)
    plot_name = plot_Dirname+'bmask_GNILC-haslam_nside256_queryRadius_%0.1f.pdf'%(q_rad)
    plt.savefig(plot_name)
    
    #map_name = map_Dirname + 'GNILC_Binary_maskregion_%d_bmask_GNILC-haslam_nside256_queryRadius_%0.1f.fits'%(region_number, q_rad)
    map_name = map_Dirname + 'GNILC_Binary_bmask_GNILC-haslam_nside256_queryRadius_%0.1f.fits'%(q_rad)
    hp.write_map(map_name, bmask_binary, overwrite=True)


def main(seq_name, masking):
    
    """
    Here earlier we were using Component seperation maps like Commander
    Which takes earlier templet map as prior like Haslam 408 MHz map 
    as prior for Commander synchrotron map.
    To see how noise behave instead of using Commander maps we are using 
    actual data maps like GNILC for Dust and Haslam 408 MHz 2014 verion 
    for Cross-Correlation. 
    """
    # reading synchrotron and dust maps
    
    fits_filename = "../CMB_foreground_map/COM_CompMap_%s-GNILC-F857_2048_R2.00.fits" % seq_name[0]
    map_1 = hp.read_map(fits_filename)
    fits_filename = "../CMB_foreground_map/haslam408_dsds_Remazeilles2014.fits" # Actual Haslam map Nside 512
    map_2 = hp.read_map(fits_filename)
    
    map_1 = hp.sphtfunc.smoothing(map_1, fwhm=np.radians(60./60.)) # smoothing Dust map to 1 degree
    map_2 = hp.sphtfunc.smoothing(map_2, fwhm=np.radians(60./60.))

    map_1 = hp.ud_grade(map_1, 256)  # degrading GNILC Dust map to Nside 256
    map_2 = hp.ud_grade(map_2, 256)  # degrading Haslam map to Nside 256

    rho_binary = np.zeros(hp.nside2npix(Nside_ref)) # Correlation Coefficient array for binary masking
    rho = np.zeros(hp.nside2npix(Nside_ref))        # Correlation Coefficient array for actual masking

    pixel_indx = np.arange(hp.nside2npix(Nside_ref)) # Nside_ref = 32 pixel index array

    for ipix in tqdm(xrange(hp.nside2npix(Nside_ref))):
        x, y, z =  hp.pix2vec(Nside_ref, ipix)
        vec = np.array([x, y, z])
        rho_binary[ipix] = cross_corr(vec, map_1, map_2)    # Compute Pearson cross-correlation coefficient
        rho[ipix] = cross_corr(vec, map_1, map_2)           # Compute Pearson cross-correlation coefficient

    # Masking Galactic part
    if masking == True:
        min_lat = 70.0
        max_lat = 110.0
        lat, lon = hp.pix2ang(Nside_ref, pixel_indx, lonlat=False)
        gal_mask = (np.rad2deg(lat) >= min_lat) * (np.rad2deg(lat) <= max_lat)
        rho_binary[gal_mask] = 0.0
        rho[gal_mask] = 0.0


#===================================================================
    # picking up pixel with correlation greater than 0.6

    indx = (rho < 0.6)
    indx1 = (rho >= 0.6)
    rho_binary[indx]  = 0.0
    rho_binary[indx1] = 1.0
    rho[indx] = hp.UNSEEN

    indx1 = (rho_binary == 1)
    pixel_indx = pixel_indx[indx1]

    map_name = map_Dirname + 'bmask_GNILC-haslam_nside16_queryRadius_%0.1f.fits'%(q_rad)
    hp.write_map(map_name, rho_binary, overwrite=True)


    map_name = map_Dirname + 'rho_GNILC-haslam_nside16_queryRadius_%0.1f.fits'%(q_rad)
    hp.write_map(map_name, rho, overwrite=True)

     
    #plot regions with high correlation
    fig = plt.figure(1, figsize=(8, 6))
    hp.mollview(rho_binary, fig=fig.number, xsize=2000, unit=r'$\rho$', nest=False, 
                title='Pearson Correlation binary map', cmap=cmap1)
    hp.graticule()
    plot_name = plot_Dirname+'CrossCor_Haslam_GNILC_rho_binary_galmask_20deg_map_nside16_nu_GE_0.6_queryradius-%0.1f.pdf'%q_rad
    plt.savefig(plot_name, dpi=600)
    
    fig = plt.figure(2, figsize=(8, 6))
    hp.mollview(hp.ma(rho), fig=fig.number, xsize=2000, unit=r'$\rho$', nest=False, 
                title='Pearson Correlation', cmap=cmap1)
    hp.graticule()
    plot_name = plot_Dirname+'CrossCor_Haslam_GNILC_rho_galmask_20deg_map_nside16_queryradius-%0.1f.pdf'%q_rad
    plt.savefig(plot_name, dpi=600)
    

    """
#===================================================================
# out of 11 regions picking up one
#region 10

    lat, lon = hp.pix2ang(Nside_ref, pixel_indx, lonlat=False)
    region1_binary = rho_binary
    
    gal_mask = (np.rad2deg(lat) < 140)
    region1_binary[gal_mask] = 0.0
    rho[gal_mask] = hp.UNSEEN

#    gal_mask = (np.rad2deg(lat) >= 60)
#    region1_binary[gal_mask] = 0.0
#    rho[gal_mask] = hp.UNSEEN

    gal_mask = (np.rad2deg(lon) > 300)
    region1_binary[gal_mask] = 0.0
    rho[gal_mask] = hp.UNSEEN

    gal_mask = (np.rad2deg(lon) < 230)
    region1_binary[gal_mask] = 0.0
    rho[gal_mask] = hp.UNSEEN

    gal_mask = (rho < 0.6)
    rho[gal_mask] = hp.UNSEEN


#===================================================================

    #plotting region 

    fig = plt.figure(3, figsize=(8, 6))
    hp.mollview(region1_binary, fig=fig.number, xsize=2000, unit=r'$\rho$', 
                    nest=False, cmap=cmap1, title='region-%d'%region_number)
    hp.graticule()

    plot_name = plot_Dirname+'region_%d_nside32_nu_GE_0.6_binary_queryRadius_%0.1f.pdf'%(region_number, q_rad)
    plt.savefig(plot_name, dpi=200)

    fig = plt.figure(4, figsize=(8, 6))
    hp.mollview(hp.ma(rho), fig=fig.number, xsize=2000, unit=r'$\rho$', 
                nest=False, cmap=cmap1, title='region-%d'%region_number)
    hp.graticule()
    plot_name = plot_Dirname+'region_%d_nside32_nu_GE_0.6_queryRadius_%0.1f.pdf'%(region_number, q_rad)
    plt.savefig(plot_name, dpi=200)
    indx1 = (region1_binary == 1)
    pixel_indx = pixel_indx[indx1]
    """

    del map_1
    del map_2
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fits_filename = "../CMB_foreground_map/COM_CompMap_%s-GNILC-F857_2048_R2.00.fits" % seq_name[0]
    map_1 = hp.read_map(fits_filename)
    #fits_filename = "../CMB_foreground_map/COM_CompMap_%s-commander_0256_R2.00.fits" % seq_name[1]
    fits_filename = "../CMB_foreground_map/haslam408_dsds_Remazeilles2014.fits" # Actual Haslam map Nside 512
    map_2 = hp.read_map(fits_filename)
    
    map_1 = hp.sphtfunc.smoothing(map_1, fwhm=np.radians(60./60.)) # smoothing Dust map to 1 degree
    map_2 = hp.sphtfunc.smoothing(map_2, fwhm=np.radians(60./60.))
    

    map_1 = hp.ud_grade(map_1, 256)  # degrading GNILC Dust map to Nside 256
    map_2 = hp.ud_grade(map_2, 256)  # degrading Haslam map to Nside 256

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """    
    fig = plt.figure(5, figsize=(8, 6))
    hp.mollview(map_1, fig=fig.number, xsize=2000, norm='hist', unit=r'T($\mu K$)', 
                                        nest=False, cmap=cmap1, title='Dust-GNILC')
    hp.graticule()
    plot_name = plot_Dirname+'GNILC-DustMap_Nside-256_smooth-60arcmin.pdf'
    plt.savefig(plot_name, dpi=200)

    fig = plt.figure(6, figsize=(8, 6))
    hp.mollview(map_2, fig=fig.number, xsize=2000,norm='hist',  unit=r'T($\mu K$)', 
                                        nest=False, cmap=cmap1, title='Synchrotron')
    hp.graticule()
    plot_name = plot_Dirname+'Synchrotron_Nside-256_smooth-60arcmin.pdf'
    plt.savefig(plot_name, dpi=200)
    """
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # Go back to foreground maps and take out correlated regions and also create binary mask
    T_T_Corr(pixel_indx, map_1, map_2)
    

if __name__ == "__main__":

    name1 = 'Dust'
    name2 = 'Synchrotron'
    seq = [name1, name2]
    main(seq, True)





