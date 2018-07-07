import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import pymaster as nmt
from scipy.interpolate import interp1d

plt.style.use("classic")

Nside = 256
map_Dirname = '/home/sandeep/PycharmProjects/Psuedo_cl_analysis/maps/'

Sync_map = hp.read_map(map_Dirname+'COM_CompMap_Synchrotron-commander_0256_R2.00.fits', verbose=False) 
dust_map = hp.read_map(map_Dirname+'COM_CompMap_dust-commander_0256_R2.00.fits', verbose=False) 

#Initialize binning scheme with 4 ells per bandpower
b = nmt.NmtBin(Nside, nlb=4)
ell_arr = b.get_effective_ells()
ell = np.linspace(min(ell_arr), max(ell_arr), 200)

for i in xrange(1, 11):
    mask_Dirname = '/home/sandeep/PycharmProjects/Psuedo_cl_analysis/region_%d/'% i
    plot_Dirname = '/home/sandeep/PycharmProjects/Psuedo_cl_analysis/region_%d/plots/'% i
    print mask_Dirname
    print plot_Dirname
    mask = hp.read_map(mask_Dirname+'region_%d_BinaryMask_apod2deg_nside256.fits'%i, verbose=False) 

    f_Sync = nmt.NmtField(mask, [Sync_map])
    f_dust = nmt.NmtField(mask, [dust_map])
    cl_dust = nmt.compute_full_master(f_dust, f_dust, b)
    cl_sync = nmt.compute_full_master(f_Sync, f_Sync, b)

    #f1 = interp1d(ell_arr, cl_dust, kind='cubic')
    #f2 = interp1d(ell_arr, cl_sync, kind='cubic')

    #y_dust = f1(ell)
    #y_sync = f2(ell)

    plt.figure(i, figsize=(8,6))
    plt.plot(ell_arr, cl_dust[0],'o-', linewidth=2, color='red', alpha=0.7,label='dust region %d'%i)
    #plt.plot(ell, y_dust[0],'-', linewidth=2, color='r',label='fit')
    plt.loglog() 
    plt.legend(loc=3, fontsize=16, frameon=False, fancybox=False, ncol=2, shadow=False)
    plt.minorticks_on()
    plt.tick_params(axis='both', which='minor', length=8, width=1.5, labelsize=18)
    plt.tick_params(axis='both', which='major', length=12, width=1.5, labelsize=18)
    plt.gca().spines['bottom'].set_linewidth(2.0)
    plt.gca().spines['left'].set_linewidth(2.0)
    plt.gca().spines['top'].set_linewidth(2.0)
    plt.gca().spines['right'].set_linewidth(2.0)
    plt.xlabel("$\\ell$", fontsize=25, fontstyle='oblique', weight='bold')
    plt.ylabel('$C_\\ell$', fontsize=25, fontstyle='oblique', weight='bold')
    plt.tight_layout()
    plt.savefig(plot_Dirname+'region_%d_Pseudo_Cl_dust.pdf'%i, dpi=600)
    plt.close()

    plt.figure(i+1, figsize=(8,6))
    plt.plot(ell_arr, cl_sync[0], 'o-', linewidth=2, color='peru', label='Sync region %d'%i)
    #plt.plot(ell, y_sync[0],'-', linewidth=2, color='r',label='fit')
    plt.loglog() 
    plt.legend(loc=3, fontsize=16, frameon=False, fancybox=False, ncol=2, shadow=False)
    plt.minorticks_on()
    plt.tick_params(axis='both', which='minor', length=8, width=1.5, labelsize=18)
    plt.tick_params(axis='both', which='major', length=12, width=1.5, labelsize=18)
    plt.gca().spines['bottom'].set_linewidth(2.0)
    plt.gca().spines['left'].set_linewidth(2.0)
    plt.gca().spines['top'].set_linewidth(2.0)
    plt.gca().spines['right'].set_linewidth(2.0)
    plt.xlabel("$\\ell$", fontsize=25, fontstyle='oblique', weight='bold')
    plt.ylabel('$C_\\ell$', fontsize=25, fontstyle='oblique', weight='bold')
    plt.tight_layout()
    plt.savefig(plot_Dirname+'region_%d_Pseudo_Cl_Sync.pdf'%i, dpi=600)
    plt.close()
plt.show()


