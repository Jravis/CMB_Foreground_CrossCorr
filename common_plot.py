"""
Comman Plot for NaMaster and Xpol
with Cosmic error
"""
import numpy as np
import lmfit
import healpy as hp
#from  multiprocessing import Process
#from scipy import special
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import scipy.stats
from tqdm import *
import pymaster as nmt
import matplotlib.gridspec as gridspec
import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')

mpl.rcParams['axes.linewidth'] = 1.5
colombi1_cmap = ListedColormap(np.loadtxt("Planck_Parchment_RGB.txt")/255.)
colombi1_cmap.set_bad("gray")
colombi1_cmap.set_under("white")
cmap1 = colombi1_cmap


def main():
    

    Xpol_dir = '/jbodstorage/data_sandeep/sandeep/cross-corr_data/Xpol_qrad_8deg/dust_Pol/'
    Dir_name = '/jbodstorage/data_sandeep/sandeep/cross-corr_data/nu_GE_0.8_ns256_dns8_qrad_8deg/'
    pixels = np.genfromtxt(Dir_name+'PseudoCl_nu_GE_0.8_ns256_dns8_qrad_8deg_pixel_index.txt',usecols=0)
    delta_ell = 12

    
    for ipix in tqdm(xrange(2, 3)):#len(pixels))):

        fname1 = Dir_name+'PseudoCl_nu_GE_0.8_ns256_dns8'\
                          '_qrad_8deg_delta_ell%d_pixindx_%d.txt'%(delta_ell, pixels[ipix])

        fname2 = Xpol_dir + 'XpolCl_dust_nu_GE_0.8_ns256_dns8_'\
                            'qrad_8deg_delta_ell%d_pixindx_%d.txt'% (delta_ell, pixels[ipix])


        ell_arr_master = np.genfromtxt(fname1, usecols=0)
        master_cl      = np.genfromtxt(fname1, usecols=1)
        err_master     = np.genfromtxt(fname1, usecols=3)

        ell_arr_Xpol   = np.genfromtxt(fname2, usecols=0)
        cl_Xpol        = np.genfromtxt(fname2, usecols=2)
        err_Xpol       = np.genfromtxt(fname2, usecols=3)


        index = (ell_arr_master>=40) * (ell_arr_master<=160)*(master_cl>0)
        ell_arr_master =   ell_arr_master[index] 
        master_cl      =   master_cl[index]      
        err_master     =   err_master[index]     

        factor_master = ell_arr_master*(ell_arr_master+1)/2./np.pi

        index = (ell_arr_Xpol>=40) * (ell_arr_Xpol<=160)*(cl_Xpol>0)
        ell_arr_Xpol   = ell_arr_Xpol[index]   
        cl_Xpol        = cl_Xpol[index]        
        err_Xpol       = err_Xpol[index]       


        factor_Xpol = ell_arr_Xpol*(ell_arr_Xpol+1)/2./np.pi

        fig = plt.figure(1, figsize=(8, 6.5))
        gs = gridspec.GridSpec(1, 1, hspace=0, wspace=0)
        ax1 = plt.subplot(gs[0, 0])

        ax1.errorbar(ell_arr_master, factor_master*master_cl, yerr=factor_master*err_master, marker='o', 
                    mfc='none', markersize=9, mew=2,mec='b', 
                    ecolor='b',elinewidth=1.5, capsize=1, 
                    color='b', label=r'$Dust-master$', linestyle='', alpha=0.7)

        ax1.errorbar(ell_arr_Xpol, factor_Xpol*cl_Xpol, yerr=factor_Xpol*err_Xpol, marker='o', 
                    mfc='none', markersize=9, mew=2,mec='m', 
                    ecolor='m',elinewidth=1.5, capsize=1, 
                    color='m', label=r'$Dust-Xpol$', linestyle='', alpha=0.7)

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
        plt.tight_layout()

#        plt.savefig(Dir_name+'plots/Dust_GE_0.8_ns256_dns8_qrad_8deg'\
#                        '_delta_ell%d_pixindx_%d.pdf'% (delta_ell, pixels[ipix]))
        
        plt.show()
        #plt.close()
        

if __name__ == "__main__":
    main()

