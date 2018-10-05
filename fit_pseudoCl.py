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

Dir_name = '/jbodstorage/data_sandeep/sandeep/cross-corr_data/nu_GE_0.8_ns256_dns8_qrad_10deg/'
pixels = np.genfromtxt(Dir_name+'PseudoCl_nu_GE_0.8_ns256_dns8_qrad_10deg_pixel_index.txt',usecols=0)
delta_ell = 12
#b = nmt.NmtBin(256, nlb=delta_ell)
#ell_arr_bin = b.get_effective_ells()
#print ell_arr_bin

#===========================================================

def Plaw(x, amp, k):
    "two Parameters: amplitude (A), and exponent (k), in: A x**k"
    return (amp * x**k)

#===========================================================

#    pars = lmfit.Parameters()
#    pars.add('amp', value=10**2, min=10, max=10**4)
#    pars.add('k', value=2.5, min=1.0, max=3.0)



def best_fit(data):


    def residual(params):

        amp = params['amp'].value
        k = params['k'].value

        model = amp*(data[0, :])**(-k)

        return (data[1, :]-model)/data[2,:]


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

#===========================================================

def mask_plot():
    
    map_2 = hp.read_map(Dir_name+'bmask_nu_GE_0.8_ns256_dns8_qrad_10deg'\
                                            '_delta_ell%d_pixindx_%d.fits'% 
                                    (delta_ell, pixels[ipix]), verbose=False)
    fig = plt.figure(1, figsize=(8, 6))
    hp.mollview(map_2, fig=fig.number, xsize=2000,norm='hist',  unit='', 
                                        nest=False, cmap=cmap1, title='CPix-%d'%pixels[ipix])
    plt.savefig(Dir_name+'bmask_nu_GE_0.8_ns256_dns8_qrad_10deg'\
                        '_delta_ell%d_pixindx_%d.pdf'% (delta_ell, pixels[ipix]))
    plt.close()

#===========================================================

def main():
    
    #mask_plot()
    #gmodel = Model(Plaw)

    for ipix in tqdm(xrange(len(pixels))):

        fname1 = Dir_name+'PseudoCl_nu_GE_0.8_ns256_dns8'\
                          '_qrad_10deg_delta_ell%d_pixindx_%d.txt'%(delta_ell, pixels[ipix])

        ell_arr_bin = np.genfromtxt(fname1, usecols=0)
        cl_dust     = np.genfromtxt(fname1, usecols=1)
        cl_sync     = np.genfromtxt(fname1, usecols=2)
        err_dust    = np.genfromtxt(fname1, usecols=3)
        err_sync    = np.genfromtxt(fname1, usecols=4)

        index = (ell_arr_bin>=40) * (ell_arr_bin<=160)*(cl_sync>0)

        ell_arr_bin =  ell_arr_bin[index] 
 #       cl_dust     =  cl_dust [index] 
        cl_sync     =  cl_sync[index] 
 #       err_dust    =  err_dust[index] 
        err_sync    =  err_sync[index] 
       
        Data = np.zeros((3, len(ell_arr_bin)), dtype=np.float64)

        factor = ell_arr_bin*(ell_arr_bin+1)/2./np.pi

        Data[0,:] = (ell_arr_bin)
        Data[1,:] = ( factor*cl_sync)
        Data[2,:] = (factor*err_sync)
    
        Amp, Amp_err, kk, kk_err = best_fit(Data)
        
        fig = plt.figure(ipix+1, figsize=(8, 6.5))
        gs = gridspec.GridSpec(1, 1, hspace=0, wspace=0)
        ax1 = plt.subplot(gs[0, 0])

        ax1.errorbar(Data[0,:], Data[1,:], yerr=Data[2,:], marker='o', 
                    mfc='none', markersize=9, mew=2,mec='teal', 
                    ecolor='teal',elinewidth=1.5, capsize=1, 
                    color='teal', label=r'$Sync$', linestyle='')#, alpha=0.7)

        ax1.plot(Data[0,:], Amp*Data[0,:]**(-kk), color='peru', linestyle='--', 
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
        plt.tight_layout()
#        plt.savefig(Dir_name+'plots/Dust_GE_0.8_ns256_dns8_qrad_10deg'\
#                        '_delta_ell%d_pixindx_%d.pdf'% (delta_ell, pixels[ipix]))

        plt.savefig(Dir_name+'plots/Sync_GE_0.8_ns256_dns8_qrad_10deg'\
                        '_delta_ell%d_pixindx_%d.pdf'% (delta_ell, pixels[ipix]))
        
        plt.close()


if __name__ == "__main__":
    main()
#plt.show()

