import matplotlib.pyplot as plt
import numpy as np
import healpy as hp


plt.style.use("classic")

#==========================================
region_number = 10
plot_Dirname= '/home/tolstoy/Documents/CMB_dust_sync_Cross_Corr/plots/regions/region_%d/'%(region_number)
map_Dirname = '/home/tolstoy/Documents/CMB_dust_sync_Cross_Corr/results/regions/region_%d/'%(region_number)

theta_ap = 2.0

f_name3 = map_Dirname + 'cl_dust_region_%d_%0.1fdeg_apodi.fits' % (region_number, theta_ap)
f_name4 = map_Dirname + 'cl_Sync_region_%d_%0.1fdeg_apodi.fits' % (region_number, theta_ap)


NSIDE=256
beam = hp.gauss_beam(np.radians(60.0/60.0), lmax=3*NSIDE-1)
pixel_window = hp.pixwin(256, pol=False)
print len(pixel_window), len(beam)

spice_cl_dust=  hp.fitsfunc.read_cl(f_name3)
spice_cl_Sync=  hp.fitsfunc.read_cl(f_name4)
plt.style.use("classic")

LMAX = 3*NSIDE-1
ell = np.arange(0, LMAX+1)


fig = plt.figure(1, figsize=(8, 6))

plt.plot(ell, ell * (ell + 1) * (spice_cl_Sync*beam**2*pixel_window[0 : 3*256])/np.pi/2., 
        '-o', color='teal', linewidth=2, label='Spice Sync' )

plt.yscale("log")
plt.xscale("log")
#plt.ylim(1e12,)
#plt.xlim(2,400)
plt.xlim(2,256)
plt.xlabel(r'$\ell$', fontsize=25, fontstyle='italic', weight='extra bold')
plt.ylabel(r'$\ell(\ell+1)C_{\ell}b_{\ell}^{2}W_{\ell}/(2\pi)$', fontsize=25, fontstyle='italic', weight='extra bold')

plt.legend(loc=1, fontsize=18, frameon=False, fancybox=False, ncol=1, shadow=False)
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=8, width=1.5, labelsize=16)
plt.tick_params(axis='both', which='major', length=15, width=1.5, labelsize=16)
plt.gca().spines['bottom'].set_linewidth(1.5)
plt.gca().spines['left'].set_linewidth(1.5)
plt.gca().spines['top'].set_linewidth(1.5)
plt.gca().spines['right'].set_linewidth(1.5)
plt.tight_layout()
name = plot_Dirname+'region_%d_Cl_apd2deg_Sync.pdf'%region_number
plt.savefig(name, dpi=800)


fig = plt.figure(2, figsize=(8, 6))

plt.plot(ell, ell * (ell + 1) * (spice_cl_dust*beam**2*pixel_window[0 : 3*256])/np.pi/2., 
        '-o', color='peru', linewidth=2, label='Spice dust' )

plt.yscale("log")
plt.xscale("log")
#plt.ylim(1,)
plt.xlim(2,256)
plt.xlabel(r'$\ell$', fontsize=25, fontstyle='italic', weight='extra bold')
plt.ylabel(r'$\ell(\ell+1)C_{\ell}b_{\ell}^{2}W_{\ell}/(2\pi)$', fontsize=25, fontstyle='italic', weight='extra bold')
plt.legend(loc=1, fontsize=18, frameon=False, fancybox=False, ncol=1, shadow=False)
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=8, width=1.5, labelsize=16)
plt.tick_params(axis='both', which='major', length=15, width=1.5, labelsize=16)
plt.gca().spines['bottom'].set_linewidth(1.5)
plt.gca().spines['left'].set_linewidth(1.5)
plt.gca().spines['top'].set_linewidth(1.5)
plt.gca().spines['right'].set_linewidth(1.5)
plt.tight_layout()

name = plot_Dirname+'region_%d_Cl_apd2deg_dust.pdf'%region_number
plt.savefig(name, dpi=800)





plt.show()






