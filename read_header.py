import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib as mpl

from matplotlib.colors import ListedColormap
mpl.rcParams['axes.linewidth'] = 1.5
colombi1_cmap = ListedColormap(np.loadtxt("Planck_Parchment_RGB.txt")/255.)
colombi1_cmap.set_bad("gray")
colombi1_cmap.set_under("white")
cmap1 = colombi1_cmap

#'======================================'

def apodize_mask(mask, delta_c):

    npix = len(mask)
    nside = hp.npix2nside(npix)

    ipix = np.arange(npix)
    vec = np.array(hp.pix2vec(nside, ipix))


    apodized_mask = 1.0*mask

    good_pix = np.where(mask == 1)[0]
    bad_pix = mask == 0

    vec = np.array(hp.pix2vec(nside, ipix))
    vec_bad = vec[:, bad_pix]

    for tmp in good_pix:

        vec_good = vec[:, tmp]

        ctheta = vec_bad[0, :]*vec_good[0] + vec_bad[1, :]*vec_good[1] + vec_bad[2, :]*vec_good[2]
        ctheta[ctheta > 1] = 1.0
        ctheta[ctheta < -1] = -1.0

        theta = np.arccos(ctheta)
        delta_i = np.degrees(np.min(theta))
        if delta_i <= delta_c:
            apodized_mask[tmp] = -1.0 / (2*np.pi) * np.sin(2*np.pi*delta_i / delta_c) + delta_i / delta_c

    return apodized_mask



print hp.nside2resol(256,arcmin=True)

#fits_filename = 'lambda_chipass_healpix_r10.fits'
fits_filename = '../maps/COM_CompMap_QU-synchrotron-commander_2048_R3.00_hm1.fits'
hdul = fits.open(fits_filename)
hdul.info()
print '========================'
print repr(hdul[1].header)


map_Q = hp.read_map(fits_filename, field=0)
map_U = hp.read_map(fits_filename, field=1)
hp.mollview(map_Q, coord=['G'], cmap=cmap1, norm='hist', title='Pol-Q')
hp.mollview(map_U, coord=['G'], cmap=cmap1, norm='hist', title='Pol-U')


fits_filename = '../maps/lambda_chipass_healpix_r10.fits'
chipass = hp.read_map(fits_filename, field=0)
index = (chipass <0)
chipass[index] = hp.UNSEEN

hp.mollview(hp.ma(chipass), coord=['G'], cmap=cmap1, norm='hist')

plt.show()





