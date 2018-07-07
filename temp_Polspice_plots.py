"""
Generalized plotting routine for 
Cross-Corr analysis
@author Sandeep Rana
"""



import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from matplotlib.colors import ListedColormap
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = 3.0
colombi1_cmap = ListedColormap(np.loadtxt("../Planck_Parchment_RGB.txt")/255.)
colombi1_cmap.set_bad("gray")
colombi1_cmap.set_under("white")
cmap1 = colombi1_cmap

plt.style.use("classic")


#This Part is for plotting common binary mask plot.
#==========================================

def bmask_plots():

    fig = plt.figure(1, figsize=(15, 15))
    plot_Dirname= '/home/tolstoy/Documents/CMB_dust_sync_Cross_Corr/plots/regions/'

    for region_number in xrange(1, 11):

        map_Dirname = '/home/tolstoy/Documents/CMB_dust_sync_Cross_Corr/results/regions/region_%d/'%(region_number)
        bmask_apod_name = map_Dirname + 'region_%d_BinaryMask_apod2deg_nside256.fits'%region_number
        bmask_apod = hp.read_map(bmask_apod_name, verbose=False)
        hp.mollview(bmask_apod, fig=fig.number, xsize=2000, unit=r'$\rho$', nest=False, cmap=cmap1,
                    sub=(4,3, region_number) ,title='region-%d'%region_number)

    hp.graticule()
    plot_name = plot_Dirname+'common_bmask_2deg.pdf'
    plt.savefig(plot_name, dpi=800, bbox_inches='tight', transparent=True)
    


#This Part is for plotting Polspice results
#==========================================
def Spice_plots():



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



def main():


    key_spice=True
    key_binary=True
    print "This routine is for various plotting"
    print ""
    print "=========================================="

    print "Do you want Spice plots[Y/n]"
    
    k = raw_input("")
    if k=='n':
        key_spice = False

    if key_spice==True:
        print "plotting spice"
        Spice_plots()        

    print "Do you want common binary mask plots[Y/n]"
    
    k = raw_input("")
    if k=='n':
        key_binary = False

    if key_binary==True:
        bmask_plots()        



if __name__ == "__main__":
    main()


plt.show()

