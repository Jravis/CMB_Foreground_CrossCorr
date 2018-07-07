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

mpl.rcParams['axes.linewidth'] = 3.0
colombi1_cmap = ListedColormap(np.loadtxt("../Planck_Parchment_RGB.txt")/255.)
colombi1_cmap.set_bad("gray")
colombi1_cmap.set_under("white")
cmap1 = colombi1_cmap


# Global values


Nside_map = 256
Nside_ref = 32
query_radius = np.deg2rad(10)

Npix_map = hp.nside2npix(Nside_map)
Npix_ref = hp.nside2npix(Nside_ref)

plt.style.use("classic")

def cross_corr(vec_arr, map1, map2):

    pix_indx_arr = hp.query_disc(Nside_map, vec_arr, query_radius)
    T1 = map1[pix_indx_arr]
    T2 = map2[pix_indx_arr]
    RCoeff, temp = scipy.stats.pearsonr(T1, T2)

    return RCoeff


def T_T_Corr(vec_arr1, map1, map2, count):

    pix_indx_arr1 = hp.query_disc(Nside_map, vec_arr1, query_radius)
    T1 = map1[pix_indx_arr1]
    T2 = map2[pix_indx_arr1]

    T1 = T1-np.mean(T1)
    T2 = T2-np.mean(T2)

    print "Pearson Coeff"
#    print np.corrcoef(T1, T2)
    print scipy.stats.pearsonr(T1, T2)

    lat, lon = hp.vec2ang(vec_arr1, lonlat=False)

    lat = np.rad2deg(lat)
    lon = np.rad2deg(lon)

    bool_arr = np.ones(hp.nside2npix(256), dtype=bool)
    for ind in pix_indx_arr1:
        bool_arr[ind]=False

    map1[bool_arr] =  hp.UNSEEN
    map2[bool_arr] =  hp.UNSEEN

    map1 = np.ma.masked_values(map1, value=-1.6375e+30)
    map2 = np.ma.masked_values(map2, value=-1.6375e+30)

    name = "../plots/dust_Synch_gnomeview_nside32_r5_30degGalCut/Cross_Corr_gnom_map_r5_dust-Synchrotron_%d.png" % count


    fig_size_inch = 10, 10
    fig = plt.figure(2, figsize=fig_size_inch)
    hp.gnomview(map1, fig=fig.number, rot=(lon, 90.-lat), xsize=1600, cmap=cmap1, flip="astro", unit=r'$T_{dust}$', nest=False,
               title='Dust', sub=(2, 2, 1))
    hp.gnomview(map2, fig=fig.number, rot=(lon, 90-lat), xsize=1600, cmap=cmap1, flip="astro", unit=r'$T_{Sync}$', nest=False,
               title='Synchrotron', sub=(2, 2, 2))
    plt.subplot(2, 2, 3)
    plt.plot(T1, T2, mec='m',marker='.', ms=3, linestyle='None')
    plt.minorticks_on()
    plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
    plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
    plt.xlabel(r"$\delta T_{d}$", fontsize=18, fontweight='bold', fontstyle='oblique')
    plt.ylabel(r"$\delta T_{Sync}$", fontsize=18, fontweight='bold', fontstyle='oblique')
    plt.savefig(name, dpi=800, bbox_inches="tight")
    plt.close()
    del bool_arr


def main(seq_name, masking):

    fits_filename = "../CMB_foreground_map/COM_CompMap_%s-commander_0256_R2.00.fits" % seq_name[0]
    map_1 = hp.read_map(fits_filename)
    fits_filename = "../CMB_foreground_map/COM_CompMap_%s-commander_0256_R2.00.fits" % seq_name[1]
    map_2 = hp.read_map(fits_filename)

    rho = np.zeros(hp.nside2npix(Nside_ref))
    pixel_indx = np.arange(hp.nside2npix(Nside_ref))

    for ipix in xrange(hp.nside2npix(Nside_ref)):
        x, y, z =  hp.pix2vec(Nside_ref, ipix)
        vec = np.array([x, y, z])
        rho[ipix] = cross_corr(vec, map_1, map_2)

    if masking == True:

        min_lat = 60.0
        max_lat = 120.0
        lat, lon = hp.pix2ang(Nside_ref, pixel_indx, lonlat=False)
        gal_mask = (np.rad2deg(lat) >= min_lat) * (np.rad2deg(lat) <= max_lat)
        rho[gal_mask] = 0.0

    indx = (rho < 0.6)
    indx1 = (rho >= 0.6)
    rho[indx]  = 0.0
    rho[indx1] = 1.0

    """
    pixel_indx = pixel_indx[indx1]
    icount = 1
    for ipix in xrange(len(pixel_indx)):

        fits_filename = "../CMB_foreground_map/COM_CompMap_%s-commander_0256_R2.00.fits" % seq_name[0]
        map_1 = hp.read_map(fits_filename, verbose=False)
        fits_filename = "../CMB_foreground_map/COM_CompMap_%s-commander_0256_R2.00.fits" % seq_name[1]
        map_2 = hp.read_map(fits_filename, verbose=False)

        x, y, z =  hp.pix2vec(Nside_ref, pixel_indx[ipix])
        vec = np.array([x, y, z])
        T_T_Corr(vec, map_1, map_2, icount)
        icount+=1
        del vec
    """

    titl = '%s-%s'%(seq_name[0], seq_name[1])
    #dpi1 = 800
    fig = plt.figure(1, figsize=(8, 6))
    hp.mollview(rho, fig=fig.number, xsize=2000, unit=r'$\rho$', nest=False, title =titl, cmap=cmap1)
    hp.graticule()
    #name = "../plots/dust_synch/rho_map_nside32_30degGalCut_cut0.6_r5"+titl+".png"
    #plt.savefig(name, dpi=dpi1, bbox_inches="tight")
    plt.show()

if __name__ == "__main__":

    name1 = 'dust'
    name2 = 'Synchrotron'
    seq = [name1, name2]
    main(seq, True)






