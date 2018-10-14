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
Nside=256
#query_radius = hp.nside2resol(Nside_ref, arcmin=False)
plt.style.use("classic")


delta_ell = 20
b = nmt.NmtBin(Nside, nlb=delta_ell)
ell_arr_bin = b.get_effective_ells()

beam = hp.gauss_beam(np.radians(40.0/60.0), lmax=3*Nside-1) # Planck synchroton has 40arcmin resolution
pix_win = hp.pixwin(Nside)
beam_binned = b.bin_cell([beam])[0]
pix_win_binned = b.bin_cell([pix_win[0:3*256]])[0]


#===============================================================================

def compute_master(f_a,f_b,wsp) :
    cl_coupled=nmt.compute_coupled_cell(f_a,f_b)
    cl_decoupled=wsp.decouple_cell(cl_coupled)
    return cl_decoupled

#===============================================================================


def compute_Polcl(mask_C1, map_q, map_u):


    ell_bin = b.get_effective_ells()

    # Create Spin-2 field

    f2_pure = 0.0
    f2_npure = 0.0
    
    f2_npure = nmt.NmtField(mask_C1, [map_q, map_u])
    f2_pure = nmt.NmtField(mask_C1,[map_q,map_u],purify_e=True,purify_b=True)

    w_np = nmt.NmtWorkspace() 
    w_np.compute_coupling_matrix(f2_npure, f2_npure,b)

    w_yp=nmt.NmtWorkspace() 
    w_yp.compute_coupling_matrix(f2_pure, f2_pure,b)

    clnp = compute_master(f2_npure, f2_npure, w_np)
    clyp = compute_master(f2_pure, f2_pure, w_yp)

    cl_EEp = np.asarray(clyp[0]/(beam_binned**2 * pix_win_binned**2))
    cl_BBp = np.asarray(clyp[3]/(beam_binned**2 * pix_win_binned**2))

    cl_EEnp = np.asarray(clnp[0]/(beam_binned**2 * pix_win_binned**2))
    cl_BBnp = np.asarray(clnp[3]/(beam_binned**2 * pix_win_binned**2))

    f_sky = (np.sum(mask_C1)*1.0)/hp.nside2npix(Nside)
    nu_b = (2*ell_arr_bin + 1) * (delta_ell*1.0) *f_sky

    Cos_Var_EEp  = cl_EEp * np.sqrt(2./nu_b)
    Cos_Var_BBp  = cl_BBp * np.sqrt(2./nu_b)
    Cos_Var_BBnp = cl_BBnp * np.sqrt(2./nu_b)
    Cos_Var_EEnp = cl_EEnp * np.sqrt(2./nu_b)

    Data = np.zeros((8, len(cl_EEp)), dtype=np.float64)

    Data[0,:] = cl_EEp
    Data[1,:] = Cos_Var_EEp 
    Data[2,:] = cl_BBp 
    Data[3,:] = Cos_Var_BBp 
    Data[4,:] = cl_EEnp 
    Data[5,:] = Cos_Var_EEnp 
    Data[6,:] = cl_BBnp 
    Data[7,:] = Cos_Var_BBnp
    
    return Data

#===============================================================================

def main():
    Dir_name = '/jbodstorage/data_sandeep/sandeep/cross-corr_data/nu_GE_0.8_ns256_dns8_qrad_10deg/'
    pixels = np.genfromtxt(Dir_name+'PseudoCl_nu_GE_0.8_ns256_dns8_qrad_10deg_pixel_index.txt',usecols=0)
    delta_ell = 20

    fits_filename = '../maps/COM_CompMap_QU-synchrotron-commander_2048_R3.00_hm1.fits'

    Sync_Q = hp.read_map(fits_filename, field=0)
    Sync_U = hp.read_map(fits_filename, field=1)
    Sync_Q = hp.sphtfunc.smoothing(Sync_Q, fwhm=np.radians(40./60.)) # smoothing Dust map to 1 degree
    Sync_U = hp.sphtfunc.smoothing(Sync_U, fwhm=np.radians(40./60.))

    Sync_Q = hp.ud_grade(Sync_Q, 256)  # degrading GNILC Dust map to Nside 256
    Sync_U = hp.ud_grade(Sync_U, 256)  # degrading Haslam map to Nside 256

#    fits_filename = '../maps/COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits'
#    Dust_Q = hp.read_map(fits_filename, field=2)
#    Dust_U = hp.read_map(fits_filename, field=3)

#    Dust_Q = hp.sphtfunc.smoothing(Dust_Q, fwhm=np.radians(40./60.)) # smoothing Dust map to 1 degree
#    Dust_U = hp.sphtfunc.smoothing(Dust_U, fwhm=np.radians(40./60.))

#    Dust_Q = hp.ud_grade(Dust_Q, 256)  # degrading GNILC Dust map to Nside 256
#    Dust_U = hp.ud_grade(Dust_U, 256)  # degrading Haslam map to Nside 256


    for ipix in tqdm(xrange(len(pixels))):    
        mask = hp.read_map(Dir_name+'mask_dir/bmask_nu_GE_0.8_ns256_dns8_qrad_10deg'\
                                            '_delta_ell12_pixindx_%d.fits'% 
                                    ( pixels[ipix]), verbose=False)
        data =  compute_Polcl(mask, Sync_Q, Sync_U)
#        data =  compute_Polcl(mask, Dust_Q, Dust_U)
        fname1 = '/jbodstorage/data_sandeep/sandeep/cross-corr_data/nu_GE_0.8_ns256_dns8_qrad_10deg/'\
                 'Sync_Pol/PseudoCl_Sync_nu_GE_0.8_ns256_dns8_qrad_10deg_delta_ell%d_pixindx_%d.txt'% (delta_ell, pixels[ipix])
        np.savetxt(fname1, data.T, delimiter='\t', fmt='%0.6e', header="cl_EEp\tCos_Var_EEp\tcl_BBp\t"\
                                                    "Cos_Var_BBp\tcl_EEnp\tCos_Var_EEnp\t"\
                                                    "cl_BBnp\tCos_Var_BBnp")

if __name__ == "__main__":
    main()














