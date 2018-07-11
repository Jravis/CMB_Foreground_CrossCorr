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
import matplotlib.gridspec as gridspec
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

    region_number = 1
    plot_Dirname= '/home/tolstoy/Documents/CMB_dust_sync_Cross_Corr/plots/regions/'
    theta_ap = 2.0
    NSIDE=256
    beam = hp.gauss_beam(np.radians(60.0/60.0), lmax=3*NSIDE-1)
    pixel_window = hp.pixwin(256, pol=False)
    LMAX = 3*NSIDE-1
    ell = np.arange(0, LMAX+1)

    fig = plt.figure(2, figsize=(15, 15))
    gs = gridspec.GridSpec(4, 3, hspace=0.3, wspace=0.3)

    for indx in xrange(0, 4):
        for indy in xrange(0, 3):

            if region_number>10:
                break
            ax1 = plt.subplot(gs[indx, indy])
            map_Dirname = '/home/tolstoy/Documents/CMB_dust_sync_Cross_Corr/results/regions/region_%d/'%(region_number)
            #f_name = map_Dirname + 'cl_dust_region_%d_%0.1fdeg_apodi.fits' % (region_number, theta_ap)
            f_name = map_Dirname + 'cl_Sync_region_%d_%0.1fdeg_apodi.fits' % (region_number, theta_ap)

            spice_cl=  hp.fitsfunc.read_cl(f_name)

            ax1.plot(ell, ell * (ell + 1) * (spice_cl*beam**2*pixel_window[0 : 3*256])/np.pi/2., 
                    '-o', color='peru', linewidth=2, label='Spice Sync' )

            ax1.set_title("region-%d"%region_number)
            ax1.set_xscale("log")
            ax1.set_yscale("log")
            #ax1.set_ylim(10.5, 14.1)
            ax1.set_xlim(2, 256)
            ax1.set_xlabel(r'$\ell$', fontsize=16, fontstyle='italic', weight='extra bold')
            ax1.set_ylabel(r'$\ell(\ell+1)C_{\ell}b_{\ell}^{2}W_{pix}/(2\pi)$', fontsize=16, 
                                     fontstyle='italic', weight='bold')
            ax1.minorticks_on()

            plt.legend(loc=1, fontsize=12, frameon=False, fancybox=False, ncol=1, shadow=False)
            plt.tick_params(axis='both', which='minor', length=8, width=1.5, labelsize=14)
            plt.tick_params(axis='both', which='major', length=15, width=1.5, labelsize=14)
            plt.gca().spines['bottom'].set_linewidth(1.5)
            plt.gca().spines['left'].set_linewidth(1.5)
            plt.gca().spines['top'].set_linewidth(1.5)
            plt.gca().spines['right'].set_linewidth(1.5)
            region_number+=1

    name = plot_Dirname+'common_Cl_apd2deg_Sync.pdf'
    plt.savefig(name, dpi=800, bbox_inches='tight')


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

