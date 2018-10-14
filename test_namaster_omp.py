import healpy as hp
import lmfit
import numpy as np
import matplotlib.pyplot as plt
import pymaster as nmt

from matplotlib.colors import ListedColormap
import matplotlib as mpl
import scipy.stats

import matplotlib.gridspec as gridspec
mpl.rcParams['axes.linewidth'] = 1.5
colombi1_cmap = ListedColormap(np.loadtxt("Planck_Parchment_RGB.txt")/255.)
colombi1_cmap.set_bad("gray")
colombi1_cmap.set_under("white")
cmap1 = colombi1_cmap

def compute_master(f_a,f_b,wsp) :
    cl_coupled=nmt.compute_coupled_cell(f_a,f_b)
    cl_decoupled=wsp.decouple_cell(cl_coupled)

    return cl_decoupled


plt.style.use("classic")

Nside=256

fits_filename = "../maps/haslam408_ns256_60arcmSmooth_dsds_Remazeilles2014.fits" 
map_1 = hp.read_map(fits_filename, verbose=False)
clarr = hp.anafast(map_1, lmax=3*Nside-1)

mask_raw = np.ones((hp.nside2npix(256)))
min_lat = 70.0
max_lat = 110.0
pixel_indx = np.arange(hp.nside2npix(256)) # Nside_ref = 256 pixel index array
lat, lon = hp.pix2ang(256, pixel_indx, lonlat=False)
gal_mask = (np.rad2deg(lat) >= min_lat) * (np.rad2deg(lat) <= max_lat)
mask_raw[gal_mask] = 0.0

delta_ell = 12
b = nmt.NmtBin(256, nlb=delta_ell)
ell_arr_bin = b.get_effective_ells()
w=nmt.NmtWorkspace()

aposcale=5.0
mask_C1 = nmt.mask_apodization(mask_raw, aposcale, apotype="C1")


f_dust = nmt.NmtField(mask_C1, [map_1])
#cl_dust = nmt.compute_full_master(f_dust, f_dust, b)
w.compute_coupling_matrix(f_dust,f_dust,b)
cl_dust=compute_master(f_dust,f_dust,w)[0]




cw=nmt.NmtCovarianceWorkspace()
cw.compute_coupling_coefficients(w,w) #<- This is the time-consuming operation
covar=nmt.gaussian_covariance(cw,clarr)#,clarr,clarr,clarr)

print covar[0,:]
print ""
print covar[10,:]

f_sky = (np.sum(mask_C1)*1.0)/hp.nside2npix(Nside)
nu_b = (2*ell_arr_bin + 1) * (delta_ell*1.0) *f_sky


print cl_dust * np.sqrt(2./nu_b)

plt.figure(); plt.imshow(covar,origin='lower',interpolation='nearest')
plt.show()

"""
print cl_dust.shape
print cl_noise.shape

#Initialize binning scheme with 4 ells per bandpower
beam = hp.gauss_beam(np.radians(60.0/60.0), lmax=3*Nside-1)
pix_win = hp.pixwin(Nside)
beam_binned = b.bin_cell([beam])[0]
pix_win_binned = b.bin_cell([pix_win[0:3*256]])[0]
cl_dust_val = np.asarray(cl_dust[0])#/(beam_binned**2 * pix_win_binned**2))
cl_dust_err = np.asarray(cl_noise[0])#/(beam_binned**2 * pix_win_binned**2))

fig = plt.figure(1, figsize=(8, 6.5))
gs = gridspec.GridSpec(1, 1, hspace=0, wspace=0)
ax1 = plt.subplot(gs[0, 0])

ax1.errorbar(ell_arr_bin, cl_dust_val, yerr=cl_dust_err, marker='o', 
            mfc='none', markersize=9, mew=2,mec='teal', 
            ecolor='teal',elinewidth=1.5, capsize=1, 
            color='teal', label=r'$Sync$', linestyle='')#, alpha=0.7)


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
"""





