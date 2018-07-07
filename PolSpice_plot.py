import lmfit
import corner
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp


plt.style.use("classic")

def fit_line(data, fig_num):
    
    pars = lmfit.Parameters()
    pars.add_many(('a', 4.), ('b', 1.3))

    def residual(p):
        a = p['a'].value
        b = p['b'].value

        return (data[:, 1] - (a+data[:, 0]*(-b)))

    mini = lmfit.Minimizer(residual, pars)

    out1 = mini.minimize(method='Nelder')

    out2 = mini.minimize(method='leastsq', params=out1.params)

    lmfit.report_fit(out2.params, min_correl=0.5)

    ci, trace = lmfit.conf_interval(mini, out2, sigmas=[0.68, 0.95],
                                trace=True, verbose=False)

    lmfit.printfuncs.report_ci(ci)
    cx, cy, grid = lmfit.conf_interval2d(mini, out2, 'a', 'b', 30, 30)
    plt.figure(fig_num, figsize=(6.5, 6))
    plt.contour(cx, cy, grid, levels=[0.68, 0.95], colors=['teal', 'peru'], linewidths=1.5)
    plt.xlabel('a', fontsize=18)
    plt.ylabel('b', fontsize=18)
    plt.legend(loc=2, fontsize=16, frameon=False)
    plt.minorticks_on()
    plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=16)
    plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=16)
    plt.tight_layout()

    return out2.params['a'].value, out2.params['b'].value


#==========================================
region_number = 1
plot_Dirname= '/home/tolstoy/Documents/CMB_dust_sync_Cross_Corr/plots/regions/region_%d/'%(region_number)
map_Dirname = '/home/tolstoy/Documents/CMB_dust_sync_Cross_Corr/results/regions/region_%d/'%(region_number)

theta_ap = 2.0

#f_name3 = map_Dirname + 'cl_dust_region_%d_%0.1fdeg_apodi.fits' % (region_number, theta_ap)
#f_name4 = map_Dirname + 'cl_Sync_region_%d_%0.1fdeg_apodi.fits' % (region_number, theta_ap)

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



print("Statistics for ell > 15.0 and ell<=300 range")
#==============================================================

y = ell * (ell + 1)*(spice_cl_dust*beam**2*pixel_window[0:3*256])/np.pi/2.
x = ell

index = (x >= 15)*(x <= 200)
data1 = np.zeros((len(x[index]), 2), dtype=np.float64)
data1[:, 0] = np.log10(x[index])
data1[:,1] = np.log10(y[index])

dust_a, dust_b = fit_line(data1, 1)

#================================================================
print"================================================================"

y = ell * (ell + 1)*(spice_cl_Sync*beam**2*pixel_window[0:3*256])/np.pi/2.
x = ell

index = (x >= 15)*(x <= 200)
data2 = np.zeros((len(x[index]), 2), dtype=np.float64)
data2[:, 0] = np.log10(x[index])
data2[:,1] = np.log10(y[index])


Sync_a, Sync_b = fit_line(data2, 2)


#================================================================
fig = plt.figure(3, figsize=(8, 6))

#plt.plot(ell, ell * (ell + 1) * (spice_cl_Sync*beam**2*pixel_window[0 : 3*256])/np.pi/2., '-o', color='teal', linewidth=2, label='Spice Sync' )

plt.plot(data2[:,0], data2[:,1], '-o', color='peru', linewidth=2, label='Spice Sync' )
plt.plot(data2[:,0], Sync_a+data2[:,0]*(-Sync_b), '-', color='teal', linewidth=2, label='bestfit' )


#plt.yscale("log")
#plt.xscale("log")
#plt.ylim(1e9,)
#plt.xlim(15,400)
plt.xlabel(r'$\ell$', fontsize=25, fontstyle='italic', weight='extra bold')
plt.ylabel(r'$\ell(\ell+1)C_{\ell}$', fontsize=25, fontstyle='italic', weight='extra bold')

plt.legend(loc=1, fontsize=18, frameon=False, fancybox=False, ncol=1, shadow=False)
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=8, width=1.5, labelsize=16)
plt.tick_params(axis='both', which='major', length=15, width=1.5, labelsize=16)
plt.gca().spines['bottom'].set_linewidth(1.5)
plt.gca().spines['left'].set_linewidth(1.5)
plt.gca().spines['top'].set_linewidth(1.5)
plt.gca().spines['right'].set_linewidth(1.5)
plt.tight_layout()

fig = plt.figure(4, figsize=(8, 6))

#plt.plot(ell, ell * (ell + 1) * (spice_cl_dust*beam**2*pixel_window[0:3*256])/np.pi/2., '-o', color='peru', linewidth=2, label='Spice dust' )
plt.plot(data1[:,0], data1[:,1], '-o', color='peru', linewidth=2, label='Spice dust' )
plt.plot(data1[:,0], dust_a+data1[:,0]*(-dust_b), '-', color='teal', linewidth=2, label='bestfit' )

#plt.yscale("log")
#plt.xscale("log")
#plt.ylim(0.1,1e4)
#plt.xlim(15,400)
plt.xlabel(r'$\ell$', fontsize=25, fontstyle='italic', weight='extra bold')
plt.ylabel(r'$\ell(\ell+1)C_{\ell}$', fontsize=25, fontstyle='italic', weight='extra bold')

plt.legend(loc=1, fontsize=18, frameon=False, fancybox=False, ncol=1, shadow=False)
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=8, width=1.5, labelsize=16)
plt.tick_params(axis='both', which='major', length=15, width=1.5, labelsize=16)
plt.gca().spines['bottom'].set_linewidth(1.5)
plt.gca().spines['left'].set_linewidth(1.5)
plt.gca().spines['top'].set_linewidth(1.5)
plt.gca().spines['right'].set_linewidth(1.5)
plt.tight_layout()


#name = plot_Dirname+'region_%d_Cl_apd2deg_dust.pdf'%region_number
#plt.savefig(name, dpi=800)
#=============================================================================================


fig = plt.figure(5, figsize=(8, 6))

plt.plot(ell, ell * (ell + 1) * (spice_cl_Sync*beam**2*pixel_window[0 : 3*256])/np.pi/2., '-o', color='teal', linewidth=2, label='Spice Sync' )

plt.yscale("log")
plt.xscale("log")
#plt.ylim(1e9,)
plt.xlim(2,400)
plt.xlabel(r'$\ell$', fontsize=25, fontstyle='italic', weight='extra bold')
plt.ylabel(r'$\ell(\ell+1)C_{\ell}$', fontsize=25, fontstyle='italic', weight='extra bold')

plt.legend(loc=1, fontsize=18, frameon=False, fancybox=False, ncol=1, shadow=False)
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=8, width=1.5, labelsize=16)
plt.tick_params(axis='both', which='major', length=15, width=1.5, labelsize=16)
plt.gca().spines['bottom'].set_linewidth(1.5)
plt.gca().spines['left'].set_linewidth(1.5)
plt.gca().spines['top'].set_linewidth(1.5)
plt.gca().spines['right'].set_linewidth(1.5)
plt.tight_layout()






plt.show()

























