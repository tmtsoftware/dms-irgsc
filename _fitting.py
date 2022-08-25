import numpy as np
from _sam import read_sam_file, select_kurucz_models, select_phoenix_models
from _extinction_correction import extinction_corrected_photometry
from matplotlib import pyplot as plt
import sys



def stdv(sfavg, v1, v2, v3, v4, v5): # function to calculate standard deviation
    mu = sfavg
    n=5
    sig = (((mu - v1)**2 + (mu - v2)**2 + (mu-v3)**2 + (mu - v4)**2 + (mu - v5)**2)/5)**0.5
    return sig

def find_nearest(self, array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def calc_sf(self, j, om, e_om, sm, index_minv):
    print('Calculating the scale factor while it is not set as a free parameter')

    ec_gmag, ec_rmag, ec_imag, ec_zmag, ec_ymag = om
    e_ec_gmag, e_ec_rmag, e_ec_imag, e_ec_zmag, e_ec_ymag = e_om
    sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k = sm

    sf_mean = ((ec_gmag[j] - sam_g[index_minv]) + (ec_rmag[j] - sam_r[index_minv]) + (ec_imag[j] - sam_i[index_minv]) + (ec_zmag[j] - sam_z[index_minv]) + (ec_ymag[j] - sam_y[index_minv]))/5.0

    e_sf_mean = (1/5)*np.sqrt(e_ec_gmag[j]**2 + e_ec_rmag[j]**2 + e_ec_imag[j]**2  + e_ec_zmag[j]**2  + e_ec_ymag[j]**2 )

    cj = sf_mean + self.aj + sam_j[index_minv] - 0.91 #0.91 is the conversion constant from J_AB to J_Vega
    ch = sf_mean + self.ah + sam_h[index_minv] - 1.39 #1.39 is the conversion constant from H_AB to H_Vega
    ck = sf_mean + self.ak + sam_k[index_minv] - 1.85 #1.85 is the conversion constant from Ks_AB to Ks_Vega

    e_cj = np.sqrt(e_sf_mean**2)# + (self.e_aj)**2)
    e_ch = np.sqrt(e_sf_mean**2)# + (self.e_ah)**2)
    e_ck = np.sqrt(e_sf_mean**2)# + (self.e_ak)**2)

    return sf_mean, e_sf_mean, cj, e_cj, ch, e_ch, ck, e_ck


def computed_reduced_chi2(self, j, oc, mc, e_oc):
    print('Compute the reduced chi2 by matching the observed and model colours')

    obs_gr, obs_gi, obs_gz, obs_gy, obs_ri, obs_ry, obs_rz, obs_iz, obs_iy, obs_zy = oc
    sam_gr, sam_gi, sam_gz, sam_gy, sam_ri, sam_rz, sam_ry, sam_iz, sam_iy, sam_zy = mc
    e_obs_gr, e_obs_gi, e_obs_gz, e_obs_gy, e_obs_ri, e_obs_ry, e_obs_rz, e_obs_iz, e_obs_iy, e_obs_zy = e_oc
    chi2_gr = ((obs_gr[j] - sam_gr)/(e_obs_gr[j]))**2
    chi2_gi = ((obs_gi[j] - sam_gi)/(e_obs_gi[j]))**2
    chi2_gz = ((obs_gz[j] - sam_gz)/(e_obs_gz[j]))**2
    chi2_gy = ((obs_gy[j] - sam_gy)/(e_obs_gy[j]))**2
    chi2_ri = ((obs_ri[j] - sam_ri)/(e_obs_ri[j]))**2
    chi2_rz = ((obs_rz[j] - sam_rz)/(e_obs_rz[j]))**2
    chi2_ry = ((obs_ry[j] - sam_ry)/(e_obs_ry[j]))**2
    chi2_iz = ((obs_iz[j] - sam_iz)/(e_obs_iz[j]))**2
    chi2_iy = ((obs_iy[j] - sam_iy)/(e_obs_iy[j]))**2
    chi2_zy = ((obs_zy[j] - sam_zy)/(e_obs_zy[j]))**2
    len_photometric_points = 10.0
    len_sam_params = 3.0
    reduced_chi2 = (1.0/(len_photometric_points - len_sam_params))*(chi2_gr + chi2_gi + chi2_gz + chi2_gy + chi2_ri + chi2_rz \
                            + chi2_ry + chi2_iz + chi2_iy + chi2_zy)
    return reduced_chi2, np.min(reduced_chi2),  chi2_gr, chi2_gi, chi2_gz, chi2_gy, chi2_ri, chi2_rz, chi2_ry, chi2_iz, chi2_iy, chi2_zy


def reduced_chi2_all_stars(self):

    dvf = []; model_params=[]; cat_ps_ra = []; cat_e_ps_ra = [];\
cat_ps_dec = []; cat_e_ps_dec = []; cat_ec_gmag = []; cat_e_ec_gmag = [];\
cat_ec_rmag = []; cat_e_ec_rmag = []; cat_ec_imag = []; cat_e_ec_imag = [];\
cat_ec_zmag = []; cat_e_ec_zmag = []; cat_ec_ymag = []; cat_e_ec_ymag = [];\
cat_sf = []; cat_e_sf = []; cat_minv = []; cat_sam_g = []; cat_sam_r = [];\
cat_sam_i = []; cat_sam_z = []; cat_sam_y = []; cat_teff = []; cat_logg = [];\
cat_feh = []; cat_computed_j = []; cat_e_computed_j = []; cat_computed_h = [];\
cat_e_computed_h = []; cat_computed_k = []; cat_e_computed_k = []


    if self.use_reduced_chi2_leq_2 is not False:
        print('Computing the NIR magnitudes for all stars by computing chi2r after matching the colours')

        ps_ra, err_ps_ra, ps_dec, err_ps_dec, ec_gmag, ec_rmag, ec_imag, ec_zmag, ec_ymag, e_ec_gmag, e_ec_rmag, e_ec_imag, e_ec_zmag, e_ec_ymag, de_reddened_gr, de_reddened_ri, de_reddened_gi, de_reddened_gz, de_reddened_gy, de_reddened_rz, de_reddened_ry, de_reddened_iz, de_reddened_iy, de_reddened_zy, e_gr, e_ri, e_gi, e_gz, e_gy, e_rz, e_ry, e_iz, e_iy, e_zy = self.extinction_corrected_photometry()

        if self.use_sam == None:
            print("Please enter the name of the Stellar Atmospheric Model to be used")

        if self.use_sam == 'Kurucz':
                print('Using Kurucz Model templates')
                model_params = self.select_kurucz_models()
                teff, logg, feh, sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k = model_params

        elif self.use_sam == 'Phoenix':
                print('Using Phoenix Model templates')
                model_params = self.select_phoenix_models()
                teff, logg, feh, sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k = model_params

        sam_gr = sam_g - sam_r; sam_ri = sam_r - sam_i; sam_gi = sam_g - sam_i;\
         sam_gz = sam_g - sam_z; sam_gy = sam_g - sam_y; sam_ry = sam_r - sam_y;\
         sam_rz = sam_r - sam_z; sam_iz = sam_i - sam_z; sam_iy = sam_i - sam_y;\
         sam_zy = sam_z - sam_y

        observed_colours = de_reddened_gr, de_reddened_gi, de_reddened_gz, de_reddened_gy, de_reddened_ri, de_reddened_ry, de_reddened_rz, de_reddened_iz, de_reddened_iy, de_reddened_zy
        model_colours = sam_gr, sam_gi, sam_gz, sam_gy, sam_ri, sam_rz, sam_ry, sam_iz, sam_iy, sam_zy
        e_observed_colours = e_gr, e_ri, e_gi, e_gz, e_gy, e_rz, e_ry, e_iz, e_iy, e_zy

        observed_optical_magnitudes = ec_gmag, ec_rmag, ec_imag, ec_zmag, ec_ymag
        e_observed_optical_magnitudes = e_ec_gmag, e_ec_rmag, e_ec_imag, e_ec_zmag, e_ec_ymag
        sam_magnitudes = sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k

        with open('catalogue.txt', 'w') as file0:
            for j in range(len(de_reddened_gr)):
                dvf, min_dvf, _, _, _, _, _, _, _, _, _, _ = self.computed_reduced_chi2(j, oc = observed_colours, mc = model_colours, e_oc = e_observed_colours)
                minv = self.find_nearest(dvf, min_dvf)
                index_best_fit_sam = np.where(minv == (dvf))[0]
                sf, e_sf, computed_j, e_computed_j, computed_h, e_computed_h, computed_k, e_computed_k = self.calc_sf(j=j, om = observed_optical_magnitudes, e_om = e_observed_optical_magnitudes, sm = sam_magnitudes, index_minv = index_best_fit_sam)
                #print('sf=', sf, e_sf, computed_j, e_computed_j, computed_h, e_computed_h, computed_k, e_computed_k)
                file0.write('%0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %s\n' %(ps_ra[j], err_ps_ra[j], ps_dec[j], err_ps_dec[j], ec_gmag[j], e_ec_gmag[j], ec_rmag[j], e_ec_rmag[j], ec_imag[j], e_ec_imag[j], ec_zmag[j], e_ec_zmag[j], ec_ymag[j], e_ec_ymag[j], sf, e_sf, minv, sam_g[index_best_fit_sam], sam_r[index_best_fit_sam], sam_i[index_best_fit_sam], sam_z[index_best_fit_sam], sam_y[index_best_fit_sam], teff[index_best_fit_sam], logg[index_best_fit_sam], feh[index_best_fit_sam], computed_j, e_computed_j, computed_h, e_computed_h, computed_k, e_computed_k, self.use_sam))
                cat_ps_ra = np.append(cat_ps_ra, ps_ra[j])
                cat_e_ps_ra = np.append(cat_e_ps_ra, err_ps_ra[j])
                cat_ps_dec = np.append(cat_ps_dec, ps_dec[j])
                cat_e_ps_dec = np.append(cat_e_ps_dec, err_ps_dec[j])
                cat_ec_gmag = np.append(cat_ec_gmag, ec_gmag[j])
                cat_e_ec_gmag = np.append(cat_e_ec_gmag, e_ec_gmag[j])
                cat_ec_rmag = np.append(cat_ec_rmag, ec_rmag[j])
                cat_e_ec_rmag = np.append(cat_e_ec_rmag, e_ec_rmag[j])
                cat_ec_imag = np.append(cat_ec_imag, ec_imag[j])
                cat_e_ec_imag = np.append(cat_e_ec_imag, e_ec_imag[j])
                cat_ec_zmag = np.append(cat_ec_zmag, ec_zmag[j])
                cat_e_ec_zmag = np.append(cat_e_ec_zmag, e_ec_zmag[j])
                cat_ec_ymag = np.append(cat_ec_ymag, ec_ymag[j])
                cat_e_ec_ymag = np.append(cat_e_ec_ymag, e_ec_ymag[j])
                cat_sf = np.append(cat_sf, sf)
                cat_e_sf = np.append(cat_e_sf, e_sf)
                cat_minv = np.append(cat_minv, minv)
                cat_sam_g = np.append(cat_sam_g, sam_g[index_best_fit_sam])
                cat_sam_r = np.append(cat_sam_r, sam_r[index_best_fit_sam])
                cat_sam_i = np.append(cat_sam_i, sam_i[index_best_fit_sam])
                cat_sam_z = np.append(cat_sam_z, sam_z[index_best_fit_sam])
                cat_sam_y = np.append(cat_sam_y, sam_y[index_best_fit_sam])
                cat_teff = np.append(cat_teff, teff[index_best_fit_sam])
                cat_logg = np.append(cat_logg, logg[index_best_fit_sam])
                cat_feh = np.append(cat_feh, feh[index_best_fit_sam])
                cat_computed_j = np.append(cat_computed_j, computed_j)
                cat_e_computed_j = np.append(cat_e_computed_j, e_computed_j)
                cat_computed_h = np.append(cat_computed_h, computed_h)
                cat_e_computed_h = np.append(cat_e_computed_h, e_computed_h)
                cat_computed_k = np.append(cat_computed_k, computed_k)
                cat_e_computed_k = np.append(cat_e_computed_k, e_computed_k)
        catalogue_params = cat_ps_ra, cat_e_ps_ra, cat_ps_dec, cat_e_ps_dec,\
         cat_ec_gmag, cat_e_ec_gmag, cat_ec_rmag, cat_e_ec_rmag, cat_ec_imag, cat_e_ec_imag,\
         cat_ec_zmag, cat_e_ec_zmag, cat_ec_ymag, cat_e_ec_ymag, cat_sf, cat_e_sf,\
         cat_minv, cat_sam_g, cat_sam_r, cat_sam_i, cat_sam_z, cat_sam_y, cat_teff, cat_logg, cat_feh,\
         cat_computed_j, cat_e_computed_j, cat_computed_h, cat_e_computed_h, cat_computed_k, cat_e_computed_k
        return catalogue_params


def reduced_chi2_leq_2(self):
    dvf = []; model_params=[]; cat_ps_ra = []; cat_e_ps_ra = [];\
cat_ps_dec = []; cat_e_ps_dec = []; cat_ec_gmag = []; cat_e_ec_gmag = [];\
cat_ec_rmag = []; cat_e_ec_rmag = []; cat_ec_imag = []; cat_e_ec_imag = [];\
cat_ec_zmag = []; cat_e_ec_zmag = []; cat_ec_ymag = []; cat_e_ec_ymag = [];\
cat_sf = []; cat_e_sf = []; cat_minv = []; cat_sam_g = []; cat_sam_r = [];\
cat_sam_i = []; cat_sam_z = []; cat_sam_y = []; cat_teff = []; cat_logg = [];\
cat_feh = []; cat_computed_j = []; cat_e_computed_j = []; cat_computed_h = [];\
cat_e_computed_h = []; cat_computed_k = []; cat_e_computed_k = []


    if self.use_reduced_chi2_leq_2 is not False:
        print('Computing the NIR magnitudes for those stars whose chi2_r is less than or \
        equal to 2 after matching the observed and model colours')

        ps_ra, err_ps_ra, ps_dec, err_ps_dec, ec_gmag, ec_rmag, ec_imag, ec_zmag, ec_ymag, e_ec_gmag, e_ec_rmag, e_ec_imag, e_ec_zmag, e_ec_ymag, de_reddened_gr, de_reddened_ri, de_reddened_gi, de_reddened_gz, de_reddened_gy, de_reddened_rz, de_reddened_ry, de_reddened_iz, de_reddened_iy, de_reddened_zy, e_gr, e_ri, e_gi, e_gz, e_gy, e_rz, e_ry, e_iz, e_iy, e_zy = self.extinction_corrected_photometry()

        if use_sam == None:
            print("Please enter the name of the Stellar Atmospheric Model to be used")
        elif use_sam == 'Kurucz':
                model_params = self.select_kurucz_models()
        elif use_sam == 'Phoenix':
                model_params = self.select_phoenix_models()
        teff, logg, feh, sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k = model_params

        sam_gr = sam_g - sam_r; sam_ri = sam_r - sam_i; sam_gi = sam_g - sam_i;\
         sam_gz = sam_g - sam_z; sam_gy = sam_g - sam_y; sam_ry = sam_r - sam_y;\
         sam_rz = sam_r - sam_z; sam_iz = sam_i - sam_z; sam_iy = sam_i - sam_y;\
         sam_zy = sam_z - sam_y

        observed_colours = de_reddened_gr, de_reddened_gi, de_reddened_gz, de_reddened_gy, de_reddened_ri, de_reddened_ry, de_reddened_rz, de_reddened_iz, de_reddened_iy, de_reddened_zy
        model_colours = sam_gr, sam_gi, sam_gz, sam_gy, sam_ri, sam_rz, sam_ry, sam_iz, sam_iy, sam_zy
        e_observed_colours = e_gr, e_ri, e_gi, e_gz, e_gy, e_rz, e_ry, e_iz, e_iy, e_zy

        with open('catalogue.txt', 'w') as file0:
            for j in range(len(de_reddened_gr)):
                print('j=', j)
                dvf, min_dvf, _, _, _, _, _, _, _, _, _, _ = self.computed_reduced_chi2(j, oc = observed_colours, mc = model_colours, e_oc = e_observed_colours)
                minv = self.find_nearest(dvf, min_dvf)
                if minv <= 2.0:
                    index_best_fit_sam = np.where(minv == (dvf))[0]
                    sf, e_sf, computed_j, e_computed_j, computed_h, e_computed_h, computed_k, e_computed_k = self.calc_sf(j, index_minv = index_best_fit_sam)
                    file0.write('%0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %s\n' %(ps_ra[j], err_ps_ra[j], ps_dec[j], err_ps_dec[j], ec_gmag[j], e_ec_gmag[j], ec_rmag[j], e_ec_rmag[j], ec_imag[j], e_ec_imag[j], ec_zmag[j], e_ec_zmag[j], ec_ymag[j], e_ec_ymag[j], sf, e_sf, minv, sam_g[index_best_fit_sam], sam_r[index_best_fit_sam], sam_i[index_best_fit_sam], sam_z[index_best_fit_sam], sam_y[index_best_fit_sam], teff[index_best_fit_sam], logg[index_best_fit_sam], feh[index_best_fit_sam], computed_j, e_computed_j, computed_h, e_computed_h, computed_k, e_computed_k, self.use_sam))
                    cat_ps_ra = np.append(cat_ps_ra, ps_ra[j])
                    cat_e_ps_ra = np.append(cat_e_ps_ra, err_ps_ra[j])
                    cat_ps_dec = np.append(cat_ps_dec, ps_dec[j])
                    cat_e_ps_dec = np.append(cat_e_ps_dec, err_ps_dec[j])
                    cat_ec_gmag = np.append(cat_ec_gmag, ec_gmag[j])
                    cat_e_ec_gmag = np.append(cat_e_ec_gmag, e_ec_gmag[j])
                    cat_ec_rmag = np.append(cat_ec_rmag, ec_rmag[j])
                    cat_e_ec_rmag = np.append(cat_e_ec_rmag, e_ec_rmag[j])
                    cat_ec_imag = np.append(cat_ec_imag, ec_imag[j])
                    cat_e_ec_imag = np.append(cat_e_ec_imag, e_ec_imag[j])
                    cat_ec_zmag = np.append(cat_ec_zmag, ec_zmag[j])
                    cat_e_ec_zmag = np.append(cat_e_ec_zmag, e_ec_zmag[j])
                    cat_ec_ymag = np.append(cat_ec_ymag, ec_ymag[j])
                    cat_e_ec_ymag = np.append(cat_e_ec_ymag, e_ec_ymag[j])
                    cat_sf = np.append(cat_sf, sf)
                    cat_e_sf = np.append(cat_e_sf, e_sf)
                    cat_minv = np.append(cat_minv, minv)
                    cat_sam_g = np.append(cat_sam_g, sam_g[index_best_fit_sam])
                    cat_sam_r = np.append(cat_sam_r, sam_r[index_best_fit_sam])
                    cat_sam_i = np.append(cat_sam_i, sam_i[index_best_fit_sam])
                    cat_sam_z = np.append(cat_sam_z, sam_z[index_best_fit_sam])
                    cat_sam_y = np.append(cat_sam_y, sam_y[index_best_fit_sam])
                    cat_teff = np.append(cat_teff, teff[index_best_fit_sam])
                    cat_logg = np.append(cat_logg, logg[index_best_fit_sam])
                    cat_feh = np.append(cat_feh, feh[index_best_fit_sam])
                    cat_computed_j = np.append(cat_computed_j, computed_j)
                    cat_e_computed_j = np.append(cat_e_computed_j, e_computed_j)
                    cat_computed_h = np.append(cat_computed_h, computed_h)
                    cat_e_computed_h = np.append(cat_e_computed_h, e_computed_h)
                    cat_computed_k = np.append(cat_computed_k, computed_k)
                    cat_e_computed_k = np.append(cat_e_computed_k, e_computed_k)
        catalogue_params = cat_ps_ra, cat_e_ps_ra, cat_ps_dec, cat_e_ps_dec,\
         cat_ec_gmag, cat_e_ec_gmag, cat_ec_rmag, cat_e_ec_rmag, cat_ec_imag, cat_e_ec_imag,\
         cat_ec_zmag, cat_e_ec_zmag, cat_ec_ymag, cat_e_ec_ymag, cat_sf, cat_e_sf,\
         cat_minv, cat_sam_g, cat_sam_r, cat_sam_i, cat_sam_z, cat_sam_y, cat_teff, cat_logg, cat_feh,\
         cat_computed_j, cat_e_computed_j, cat_computed_h, cat_e_computed_h, cat_computed_k, cat_e_computed_k
        return catalogue_params

def compute_flux(self, mag):
    return 10**(-0.4*(48.6 + mag))

def calc_sf2(self, initial_guess, smag, obmag, e_obmag, e_ebv, reddening_constants):

    if self.use_compute_nir_by_keeping_sf_and_reddening_free is True:
        print('Computing the NIR magnitudes by keeping the scale factor and reddening free')

        rg, rr, ri, rz, ry = reddening_constants
        sf, ebv_fp = initial_guess

        sgmag, srmag, simag, szmag, symag = smag
        gobmag, robmag, iobmag, zobmag, yobmag = obmag
        e_gobmag, e_robmag, e_iobmag, e_zobmag, e_yobmag = e_obmag

        gobflux = self.compute_flux(gobmag); robflux = self.compute_flux(robmag); \
        iobflux = self.compute_flux(iobmag); zobflux = self.compute_flux(zobmag); \
        yobflux = self.compute_flux(yobmag)

        e_gobflux = gobflux*10**(0.4*(e_gobmag)); e_robflux = robflux*10**(0.4*(e_robmag)); e_iobflux = iobflux*10**(0.4*(e_iobmag));\
        e_zobflux = zobflux*10**(0.4*(e_zobmag)); e_yobflux = yobflux*10**(0.4*(e_yobmag))


        gsflux =  self.compute_flux(sgmag); rsflux =  self.compute_flux(srmag);\
        isflux =  self.compute_flux(simag); zsflux =  self.compute_flux(szmag);\
        ysflux =  self.compute_flux(symag)

        chi2r = (1/2)*(((gobflux - pow(10,sf)*gsflux*pow(10,0.4*(ebv_fp*rg)))/(e_gobflux))**2 + ((robflux - pow(10,sf)*rsflux*pow(10,0.4*(ebv_fp*rr)))/(e_robflux))**2 + ((iobflux - pow(10,sf)*isflux*pow(10,0.4*(ebv_fp*ri)))/(e_iobflux))**2 + ((zobflux - pow(10,sf)*zsflux*pow(10,0.4*(ebv_fp*rz)))/(e_zobflux))**2 + ((yobflux - pow(10,sf)*ysflux*pow(10,0.4*(ebv_fp*ry)))/(e_yobflux))**2)
        return chi2r

import scipy.optimize
#function to calculate the NIR magnitudes by keeping the scale factor and reddening free parameters
def compute_nir2(self):

    bnds = [(-np.inf,np.inf), (0,2)]
    #ps_ra, e_ps_ra, ps_dec, e_ps_dec, ec_gmag, ec_rmag, ec_imag, ec_zmag, ec_ymag, e_ec_gmag, e_ec_rmag, e_ec_imag, e_ec_zmag, e_ec_ymag, de_reddened_gr, de_reddened_gi, de_reddened_ri, de_reddened_gy, de_reddened_gz, de_reddened_ry, de_reddened_rz, de_reddened_iy, de_reddened_iz, de_reddened_zy, e_ec_gr, e_ec_gi, e_ec_gz, e_ec_gy, e_ec_ri, e_ec_rz, e_ec_ry, e_ec_iz, e_ec_iy, e_ec_zy = self.extinction_corrected_photometry()
    ps_ra, e_ps_ra, ps_dec, e_ps_dec, ps1g, ps1r, ps1i, ps1z, ps1y, ps1_eg, ps1_er, ps1_ei, ps1_ez, ps1_ey = self.star_galaxy_classification()

    gobmag = ps1g; robmag = ps1r; iobmag = ps1i; zobmag = ps1z; yobmag = ps1y
    e_gobmag = ps1_eg; e_robmag = ps1_er; e_iobmag = ps1_ei; e_zobmag = ps1_ez; e_yobmag = ps1_ey

    #reddening constants adopted from Tonry et.al. 2012
    Rg = (3.613 - 0.0972*(ps1g - ps1i) + 0.01*(ps1g - ps1i)**2)
    Rr = (2.585 - 0.0315*(ps1g - ps1i))
    Ri = (1.908 - 0.0152*(ps1g - ps1i))
    Rz = (1.499 - 0.0023*(ps1g - ps1i))
    Ry = (1.251 - 0.0027*(ps1g - ps1i))

    if use_sam == None:
            print("Please enter the name of the Stellar Atmospheric Model to be used")
    elif use_sam == 'Kurucz':
        model_params = self.select_kurucz_models()
    elif use_sam == 'Phoenix':
        model_params = self.select_phoenix_models()
    teff, logg, feh, sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k = model_params

    with open('catalogue.txt', 'w') as file0:

        for j in range(len(ps_ra)):
            chi2_arr = []
            sf_arr = []
            ebv_arr = []
            for j0 in range(len(teff)):
                if use_sam == None:
                        print("Please enter the name of the Stellar Atmospheric Model to be used")

                smag = sam_g[j0], sam_r[j0], sam_i[j0], sam_z[j0], sam_y[j0]
                obmag = gobmag[j], robmag[j], iobmag[j], zobmag[j], yobmag[j]
                e_obmag = e_gobmag[j], e_robmag[j], e_iobmag[j], e_zobmag[j],e_yobmag[j]
                reddening_constants = [Rg[j], Rr[j], Ri[j], Rz[j], Ry[j]]
                op = scipy.optimize.minimize(self.calc_sf2, self.starting_guess, args=(smag, obmag, e_obmag, self.err_ebv, reddening_constants), bounds=bnds)
                res = op.x
                sf_arr = np.append(sf_arr, res[0])
                ebv_arr = np.append(ebv_arr, res[1])
                chi2_arr = np.append(chi2_arr, self.calc_sf2(res, smag, obmag, e_obmag, self.err_ebv, reddening_constants))
            minv = np.min(chi2_arr)
            index_minv = np.where(chi2_arr==minv)[0]
            best_sf = sf_arr[index_minv]
            kfj = self.compute_flux(sam_j[index_minv])
            kfh = self.compute_flux(sam_h[index_minv])
            kfk = self.compute_flux(sam_k[index_minv])
            computed_j = -48.6 - 2.5*np.log10(pow(10,best_sf) * kfj) + self.aj - 0.91
            computed_h = -48.6 - 2.5*np.log10(pow(10,best_sf) * kfh) + self.ah - 1.39
            computed_k = -48.6 - 2.5*np.log10(pow(10,best_sf) * kfk) + self.ak -1.85
            file0.write('%0.16f %0.8f %0.8f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.5f %0.5f %s\n' %(best_sf, actual_ps_ra[j], actual_ps_dec[j], actual_err_ps_ra[j], actual_err_ps_dec[j], kurucz_ukidss_j[index_minv], actual_gmag[j], actual_rmag[j], actual_imag[j], actual_zmag[j], actual_ymag[j], teff[index_minv], logg[index_minv], feh[index_minv], actual_e_gmag[j], actual_e_rmag[j], actual_e_imag[j], actual_e_zmag[j], actual_e_ymag[j], kurucz_ps_g[index_minv], kurucz_ps_r[index_minv], kurucz_ps_i[index_minv], kurucz_ps_z[index_minv], kurucz_ps_y[index_minv], minv, err_obs[j], computed_j, computed_h, computed_k, kurucz_ukidss_h[index_minv], kurucz_ukidss_k[index_minv], ebv_arr[index_minv], err_ebv, %s))
