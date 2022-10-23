from read_data import read_optical_data
import numpy as np
from matplotlib import pyplot as plt



def star_galaxy_classification(self):
    print('Seperating the probable stellar sources from the optical dataset which is filtered\
     for nan values')
    optical_data=self.optical_data
    ps_phot = self.read_optical_data()
    print('Using psf-kron criteria to seperate stars and galaxies')
    ps_ra, e_ps_ra, ps_dec, e_ps_dec, gpsf, rpsf, ipsf, zpsf, ypsf, e_gpsf, e_rpsf, e_ipsf, e_zpsf, e_ypsf, gkron, rkron, ikron, zkron, ykron = ps_phot
    sgc_index = np.where((gpsf - gkron < 0.05) & (rpsf - rkron < 0.05) & (ipsf - ikron < 0.05) & (zpsf - zkron < 0.05) & (ypsf - ykron < 0.05))[0]
    galaxy_index = np.where((gpsf - gkron > 0.05) & (rpsf - rkron > 0.05) & (ipsf - ikron > 0.05) & (zpsf - zkron > 0.05) & (ypsf - ykron > 0.05))[0]

    #gpsf, rpsf, ipsf, zpsf, ypsf = gpsf[sgc_index], rpsf[sgc_index], ipsf[sgc_index], zpsf[sgc_index], ypsf[sgc_index]
    #e_gpsf, e_rpsf, e_ipsf, e_zpsf, e_ypsf = e_gpsf[sgc_index], e_rpsf[sgc_index], e_ipsf[sgc_index], e_zpsf[sgc_index], e_ypsf[sgc_index]

    #gkron, rkron, ikron, zkron, ykron = gkron[sgc_index], rkron[sgc_index], ikron[sgc_index], zkron[sgc_index], ykron[sgc_index]

    print('Plotting the CCD which shows stars and galaxies as magenta and black points respectively')
    print(len(ipsf), len(rpsf), len(gpsf), len(sgc_index))
    plt.clf()
    plt.scatter((gpsf[sgc_index] - rpsf[sgc_index]), (rpsf[sgc_index] - ipsf[sgc_index]), s=5, color= 'm', alpha = 0.3)
    plt.scatter((gpsf[galaxy_index] - rpsf[galaxy_index]), (rpsf[galaxy_index] - ipsf[galaxy_index]), s=5, color = 'k', alpha = 0.3)
    plt.xlabel('$(g-r)$')
    plt.ylabel('$(r-i)$')
    plt.savefig('ccd_stars_and_galaxies_seperated.png')
    plt.clf()

    print('Plotting the (ipsf-ikron) vs (ikron) scatter plot which shows stars and galaxies as magenta and black points respectively')


    plt.clf()
    plt.scatter(ipsf[sgc_index], (ipsf[sgc_index] - ikron[sgc_index]), s=5, color='m', alpha = 0.3)
    plt.scatter(ipsf[galaxy_index], (ipsf[galaxy_index] - ikron[galaxy_index]), s=5, color='k', alpha = 0.3)
    plt.xlabel('$i_{psf}$')
    plt.ylabel('$i_{psf}-i_{kron}$')
    plt.savefig('psf_vs_kron_stars_and_galaxies_seperated.png')
    plt.clf()

    psf_phot = ps_ra[sgc_index], e_ps_ra[sgc_index], ps_dec[sgc_index], e_ps_dec[sgc_index], gpsf[sgc_index], rpsf[sgc_index], ipsf[sgc_index], zpsf[sgc_index], ypsf[sgc_index], e_gpsf[sgc_index], e_rpsf[sgc_index], e_ipsf[sgc_index], e_zpsf[sgc_index], e_ypsf[sgc_index]

    print('Created an input optical catalogue of stellar sources')

    return psf_phot
