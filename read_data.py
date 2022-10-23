import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits

print('Reading input datafiles')

def read_optical_data(self):
    print('Reading optical data')
    optical_data = self.optical_data

    #p1 = np.genfromtxt('/mnt/c/Users/sshah/Documents/itcc/irgsc/script/ps_tf6.csv', delimiter = ',')
    p1 = optical_data#np.genfromtxt(optical_data, delimiter = ',', skip_header=1)

    ps_ra = p1[:,0]
    ps_dec = p1[:,1]

    ps_ra_err = (p1[:,2])
    ps_dec_err = (p1[:,3])

    gpsf = p1[:,4]
    gkron = p1[:,6]

    e_gpsf = p1[:,5]
    e_gkron = p1[:,7]

    rpsf = p1[:,8]
    rkron = p1[:,10]

    e_rpsf = p1[:,9]
    e_rkron = p1[:,11]

    ipsf = p1[:,12]
    ikron = p1[:,14]

    e_ipsf = p1[:,13]
    e_ikron = p1[:,15]

    zpsf = p1[:,16]
    zkron = p1[:,18]
    e_zpsf = p1[:,17]
    e_zkron = p1[:,19]

    ypsf = p1[:,20]
    ykron = p1[:,22]
    e_ypsf = p1[:,21]
    e_ykron = p1[:,23]

    print('Now filtering the optical data for nan values')

    ind = np.where((gpsf!= -999) & (ipsf!= -999) & (rpsf != -999) & (zpsf != -999) & (ypsf != -999) & (e_gpsf != -999) & (e_rpsf != -999) & (e_ipsf != -999) & (e_zpsf != -999) & (e_ypsf != -999) & (gkron!= -999) & (ikron!= -999) & (zkron != -999) & (ykron != -999) & (rkron != -999) & (e_gpsf < 0.2) & (e_rpsf < 0.2) & (e_ipsf < 0.2) & (e_zpsf < 0.2) & (e_ypsf < 0.2))[0]

    print('Number of Sources in the Optical Dataset=', len(ind))

    print('Plotting the psf-kron vs psf scatter plot of the sources in the optical data')

    plt.clf()
    plt.scatter(ipsf[ind], ipsf[ind] - ikron[ind], color='green', s=5, alpha=0.5)
    plt.xlabel('$i_{psf}$')
    plt.ylabel('$i_{kron}$')
    plt.grid()
    plt.xlim(10,24)
    plt.ylim(-0.8,2)
    plt.savefig('psf_vs_kron_whole_dataset.png')
    plt.clf()

    print('Plotting the CCD of the sources in the optical data')

    plt.clf()
    plt.scatter((gpsf[ind] - rpsf[ind]), (rpsf[ind] - ipsf[ind]), s=5,color='green', alpha = 0.3, label='Number of sources='+str(len(ind)))
    plt.xlabel('$(g-r)$')
    plt.ylabel('$(r-i)$')
    plt.grid()
    plt.xlim(-1.5,1.5)
    plt.ylim(-1.2,2.2)
    plt.savefig('ccd_whole_dataset.png')
    plt.clf()

    filtered_ra = ps_ra[ind]
    filtered_err_ra = ps_ra_err[ind]
    filtered_dec = ps_dec[ind]
    filtered_err_dec = ps_dec_err[ind]

    filtered_gpsf = gpsf[ind]
    filtered_rpsf = rpsf[ind]
    filtered_ipsf = ipsf[ind]
    filtered_zpsf = zpsf[ind]
    filtered_ypsf = ypsf[ind]

    filtered_e_gpsf = e_gpsf[ind]
    filtered_e_rpsf = e_rpsf[ind]
    filtered_e_ipsf = e_ipsf[ind]
    filtered_e_zpsf = e_zpsf[ind]
    filtered_e_ypsf = e_ypsf[ind]

    filtered_gkron = gkron[ind]
    filtered_rkron = rkron[ind]
    filtered_ikron = ikron[ind]
    filtered_zkron = zkron[ind]
    filtered_ykron = ykron[ind]

    raw_optical_data = filtered_ra, filtered_err_ra, filtered_dec, filtered_err_dec, filtered_gpsf, filtered_rpsf, filtered_ipsf, filtered_zpsf, filtered_ypsf, filtered_e_gpsf, filtered_e_rpsf, filtered_e_ipsf, filtered_e_zpsf, filtered_e_ypsf, filtered_gkron, filtered_rkron, filtered_ikron, filtered_zkron, filtered_ykron
    return raw_optical_data

def read_nir_data(self):
    print('vd=', self.validate)
    if self.validate == True:
        if self.validating_data == None:
            print('Error: Validating data not provided')
        print('Now reading the NIR survey data for validation')
        #hdulist = fits.open('/mnt/c/Users/sshah/Documents/itcc/irgsc/script/tf6.fits',  memmap=True)
        validating_data = self.validating_data
        hdulist = validating_data#fits.open(str(validating_data), memap=True)

        p8 = hdulist[1].data

        petro_j = p8['JPETROMAG']
        e_petro_j = p8['jPetroMagErr']
        petro_h = p8['HPETROMAG']
        e_petro_h = p8['hPetroMagErr']
        petro_k = p8['KPETROMAG']
        e_petro_k = p8['kPetroMagErr']
        t1_ra = (p8['RA'])*(180.0/np.pi)
        t1_dec = p8['DEC']*(180.0/np.pi)

        nir_index = np.where((petro_j != -9.99999500e+08) & (e_petro_j != -9.99999500e+08) & (petro_h != -9.99999500e+08) & (e_petro_h != -9.99999500e+08) & (petro_k != -9.99999500e+08) & (e_petro_k != -9.99999500e+08) & (e_petro_j < 0.1) & (e_petro_h < 0.1) & (e_petro_k < 0.1))[0]
        print('Number of Stars in the NIR data = ', len(nir_index))

        filtered_petro_j = petro_j[nir_index]
        e_petro_h = e_petro_h[nir_index]
        filtered_petro_h = petro_h[nir_index]
        e_petro_k = e_petro_k[nir_index]
        filtered_petro_k = petro_k[nir_index]
        e_petro_j = e_petro_j[nir_index]
        filtered_t1_ra = t1_ra[nir_index]
        filtered_t1_dec = t1_dec[nir_index]

        plt.clf()
        plt.hist(filtered_petro_j)
        plt.xlabel('$J_{petro}$')
        plt.ylabel('Bins')
        plt.savefig('hist_ukidss_J.png')
        plt.clf()

        plt.clf()
        plt.hist(filtered_petro_h)
        plt.xlabel('$H_{petro}$')
        plt.ylabel('Bins')
        plt.savefig('hist_ukidss_H.png')
        plt.clf()

        plt.clf()
        plt.hist(filtered_petro_k)
        plt.xlabel('$K_{petro}$')
        plt.ylabel('Bins')
        plt.savefig('hist_ukidss_K.png')
        plt.clf()

        nir_data =  filtered_petro_j, filtered_petro_h, filtered_petro_k, e_petro_j, e_petro_h, e_petro_k, filtered_t1_ra, filtered_t1_dec
        return nir_data
