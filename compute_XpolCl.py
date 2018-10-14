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
from numpy import *
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import scipy.stats
from tqdm import *
import pymaster as nmt
import xpol
import matplotlib.gridspec as gridspec
import lmfit

mpl.rcParams['axes.linewidth'] = 1.5
colombi1_cmap = ListedColormap(np.loadtxt("Planck_Parchment_RGB.txt")/255.)
colombi1_cmap.set_bad("gray")
colombi1_cmap.set_under("white")
cmap1 = colombi1_cmap


# Global values
Nside=256
#query_radius = hp.nside2resol(Nside_ref, arcmin=False)
plt.style.use("classic")

lmax=3*Nside-1
tagnames = ['TT','EE','BB','TE','TB','EB']
delta_ell=12

#===========================================================

def Plaw(x, amp, k):
    "two Parameters: amplitude (A), and exponent (k), in: A x**k"
    return (amp * x**k)

#===========================================================

#    pars = lmfit.Parameters()
#    pars.add('amp', value=10**2, min=10, max=10**4)
#    pars.add('k', value=2.5, min=1.0, max=3.0)



def best_fit(data_fit):


    def residual(params):

        amp = params['amp'].value
        k = params['k'].value

        model = amp*(data_fit[0, :])**(-k)

        return (data_fit[1, :]-model)/data_fit[2,:]


    params = lmfit.Parameters()
    params.add_many(('amp', (10**5)), ('k', (2.5)))
   
    mini = lmfit.Minimizer(residual, params)

    # first solve with Nelder-Mead
    out1 = mini.minimize(method='Nelder')

    lmfit.printfuncs.report_fit(out1.params, min_correl=0.5)
    # then solve with Levenberg-Marquardt using the
    # Nelder-Mead solution as a starting point
    
    out2 = mini.minimize(method='leastsq', params=out1.params)
    lmfit.printfuncs.report_fit(out2.params, min_correl=0.5)
    
    return out2.params['amp'].value, out2.params['amp'].stderr, out2.params['k'].value, out2.params['k'].stderr 


#===============================================================================

beam = hp.gauss_beam(np.radians(60.0/60.0), lmax=3*Nside-1)

def compute_Xpolcl(mask_raw, dT):

    aposcale = 2.0
    mask_C1 = nmt.mask_apodization(mask_raw, aposcale, apotype="C1")

    #generate binning from l=2 to lmax with deltal=1
    b = xpol.Bins.fromdeltal( 2, lmax, delta_ell)
    ell_arr_bin = b.lbin

    #generate xpol class
    xp = xpol.Xpol(mask_C1, b, polar=False)
    #compute spectra from map dT
    pcl_TT, cl_TT = xp.get_spectra( dT, bell=beam, pixwin=True)

    f_sky = (np.sum(mask_C1)*1.0)/hp.nside2npix(Nside)
    nu_b = (2*ell_arr_bin + 1) * (delta_ell*1.0) *f_sky
    Cosmic_TT  = (cl_TT) * np.sqrt(2./nu_b)
    
    Data = np.zeros((4, len(cl_TT))) 
    Data[0,:] = ell_arr_bin
    Data[1,:] = b.bin_spectra(pcl_TT)
    Data[2,:] = (cl_TT)
    Data[3,:] = Cosmic_TT

    return Data

#===============================================================================

def main():

    Dir_name = '/jbodstorage/data_sandeep/sandeep/cross-corr_data/nu_GE_0.8_ns256_dns8_qrad_8deg/'
    Xpol_dir = '/jbodstorage/data_sandeep/sandeep/cross-corr_data/Xpol_qrad_8deg/dust_Pol/'
    pixels = np.genfromtxt(Dir_name+'PseudoCl_nu_GE_0.8_ns256_dns8_qrad_8deg_pixel_index.txt',usecols=0)
    delta_ell = 12

    #map_1  =hp.read_map("../maps/haslam408_ns256_60arcmSmooth_dsds_Remazeilles2014.fits")
    map_1  =hp.read_map("../maps/GNILC857_ns256_60arcmSmooth_Dust.fits")

    for ipix in tqdm(xrange(len(pixels))):    
        mask = hp.read_map(Dir_name+'mask_dir/bmask_nu_GE_0.8_ns256_dns8_qrad_8deg'\
                                            '_delta_ell12_pixindx_%d.fits'% 
                                    ( pixels[ipix]), verbose=False)

        data =  compute_Xpolcl(mask, map_1)

        fname1 =Xpol_dir + 'XpolCl_dust_nu_GE_0.8_ns256_dns8_'\
                            'qrad_8deg_delta_ell%d_pixindx_%d.txt'% (delta_ell, pixels[ipix])

        np.savetxt(fname1, data.T, delimiter='\t', fmt='%0.6e', header="ell_bin\tpcl_TT\tcl_TT\tCos_Var_EEp")

#=============================================================================================
    
        index = (data[0,:]>=40) * (data[0,:]<=160)*(data[2,:]>0)

        ell_arr_bin = data[0,:][index]
        cl_TT = data[2,:][index]
        err_TT = data[3,:][index]

        Data_fit = np.zeros((3, len(ell_arr_bin)), dtype=np.float64)

        factor = ell_arr_bin*(ell_arr_bin+1)/2./np.pi

        Data_fit[0,:] = (ell_arr_bin)
        Data_fit[1,:] = factor*cl_TT
        Data_fit[2,:] = factor*err_TT

        Amp, Amp_err, kk, kk_err = best_fit(Data_fit)
        
        fig = plt.figure(ipix+1, figsize=(8, 6.5))
        gs = gridspec.GridSpec(1, 1, hspace=0, wspace=0)
        ax1 = plt.subplot(gs[0, 0])

        ax1.errorbar(Data_fit[0,:], Data_fit[1,:], yerr=Data_fit[2,:], marker='o', 
                    mfc='none', markersize=9, mew=2,mec='b', 
                    ecolor='b',elinewidth=1.5, capsize=1, 
                    color='b', label=r'$unbiased-C_{\ell}^{Dust}$', linestyle='', alpha=0.7)

        ax1.plot(Data_fit[0,:], Amp*Data_fit[0,:]**(-kk), color='r', linestyle='--', 
                    linewidth=2, label=r'$index(k) = %0.3f\pm%0.3f$'%(kk, kk_err) )

        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.legend()
        ax1.legend(loc=1, fontsize=16, frameon=False, fancybox=False, shadow=False)
        ax1.minorticks_on()
        ax1.set_xlim(40, 160)
        #ax1.set_ylim(0.1, 50)
        plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=18)
        plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=18)
        #ax1.spines['top'].set_visible(2.0)
        #ax1.spines['right'].set_visible(2.0)
        ax1.spines['bottom'].set_linewidth(2.0)
        ax1.spines['left'].set_linewidth(2.0)
        ax1.spines['top'].set_linewidth(2.0)
        ax1.spines['right'].set_linewidth(2.0)
        ax1.set_xlabel(r"$\ell$", fontsize=20, weight='bold')
        ax1.set_ylabel(r"${\cal D}_{\ell}/(B_{\ell}^{2}W_{\ell}^{2})$", fontsize=20, weight='bold')
        plt.title("pixel index = %d"%pixels[ipix])
        plt.tight_layout()
        plt.savefig(Xpol_dir+'plots/Dust_GE_0.8_ns256_dns8_qrad_8deg'\
                        '_delta_ell%d_pixindx_%d.pdf'% (delta_ell, pixels[ipix]))
        #plt.show()
        plt.close()
     


if __name__ == "__main__":
    main()














