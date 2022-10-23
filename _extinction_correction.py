from _sgc import star_galaxy_classification
import numpy as np


def extinction_corrected_photometry(self):
    print('Correcting the optical photometry of the probable stellar sources for reddening and extinction.')
    optical_data = self.optical_data
    ps_ra, e_ps_ra, ps_dec, e_ps_dec, ps1g, ps1r, ps1i, ps1z, ps1y, ps1_eg, ps1_er, ps1_ei, ps1_ez, ps1_ey = self.star_galaxy_classification()

    #error in observed colours
    e_gr = np.sqrt(ps1_eg**2 + ps1_er**2); e_gi = np.sqrt(ps1_eg**2 + ps1_ei**2); e_gz = np.sqrt(ps1_eg**2 + ps1_ez**2)
    e_gy = np.sqrt(ps1_eg**2 + ps1_ey**2); e_ri = np.sqrt(ps1_er**2 + ps1_ei**2); e_rz = np.sqrt(ps1_er**2 + ps1_ez**2)
    e_ry = np.sqrt(ps1_er**2 + ps1_ey**2); e_iy = np.sqrt(ps1_ei**2 + ps1_ey**2); e_iz = np.sqrt(ps1_ei**2 + ps1_ez**2)
    e_zy = np.sqrt(ps1_ez**2 + ps1_ey**2)

    #extinction in ps1 filters taken from Tonry et.al. 2012
    ag = (self.ebv)*0.88*(3.613 - 0.0972*(ps1g - ps1i) + 0.01*(ps1g - ps1i)**2)
    ar = (self.ebv)*0.88*(2.585 - 0.0315*(ps1g - ps1i)); ai = (self.ebv)*0.88*(1.908 - 0.0152*(ps1g - ps1i))
    az = (self.ebv)*0.88*(1.499 - 0.0023*(ps1g - ps1i)); ay = (self.ebv)*0.88*(1.251 - 0.0027*(ps1g - ps1i))

    #error in extinction
    e_ag = ((self.err_ebv)*(3.613 - 0.0972*(ps1g - ps1i) + 0.02*(ps1g - ps1i)**2))+\
            ((self.ebv)*((-0.0972*e_gi)+(0.02*(ps1g - ps1i)*e_gi)))
    e_ar = (self.err_ebv)*(2.585 - 0.0315*(ps1g - ps1i)) + (self.ebv)*(-0.0315*e_gi)
    e_ai = (self.err_ebv)*(1.908 - 0.0152*(ps1g - ps1i)) + (self.ebv)*(-0.0152*e_gi)
    e_az = (self.err_ebv)*(1.499 - 0.0023*(ps1g - ps1i)) + (self.ebv)*(-0.0023*e_gi)
    e_ay = (self.err_ebv)*(1.251 - 0.0027*(ps1g - ps1i)) + (self.ebv)*(-0.0027*e_gi)

    #prefix ec_ stands for extinction corrected magnitudes and e_ec_ stands for error in them
    ec_gmag = ps1g - ag; ec_rmag = ps1r - ar; ec_imag = ps1i - ai; ec_zmag = ps1z - az
    ec_ymag = ps1y - ay

    e_ec_gmag = np.sqrt((ps1_eg)**2 + (e_ag)**2)
    e_ec_rmag = np.sqrt((ps1_er)**2 + (e_ar)**2)
    e_ec_imag = np.sqrt((ps1_ei)**2 + (e_ai)**2)
    e_ec_zmag = np.sqrt((ps1_ez)**2 + (e_az)**2)
    e_ec_ymag = np.sqrt((ps1_ey)**2 + (e_ay)**2)

    e_ec_gr = np.sqrt(e_ec_gmag**2 + e_ec_rmag**2)
    e_ec_gi = np.sqrt(e_ec_gmag**2 + e_ec_imag**2)
    e_ec_gz = np.sqrt(e_ec_gmag**2 + e_ec_zmag**2)
    e_ec_gy = np.sqrt(e_ec_gmag**2 + e_ec_ymag**2)
    e_ec_ri = np.sqrt(e_ec_rmag**2 + e_ec_imag**2)
    e_ec_rz = np.sqrt(e_ec_rmag**2 + e_ec_zmag**2)
    e_ec_ry = np.sqrt(e_ec_rmag**2 + e_ec_ymag**2)
    e_ec_iz = np.sqrt(e_ec_imag**2 + e_ec_zmag**2)
    e_ec_iy = np.sqrt(e_ec_imag**2 + e_ec_ymag**2)
    e_ec_zy = np.sqrt(e_ec_zmag**2 + e_ec_ymag**2)


    de_reddened_gr = ec_gmag - ec_rmag
    de_reddened_ri = ec_rmag - ec_imag
    de_reddened_gi = ec_gmag - ec_imag
    de_reddened_gy = ec_gmag - ec_ymag
    de_reddened_gz = ec_gmag - ec_zmag
    de_reddened_ry = ec_rmag - ec_ymag
    de_reddened_rz = ec_rmag - ec_zmag
    de_reddened_iy = ec_imag - ec_ymag
    de_reddened_iz = ec_imag - ec_zmag
    de_reddened_zy = ec_zmag - ec_ymag


    psf_phot = ps_ra, e_ps_ra, ps_dec, e_ps_dec, ec_gmag, ec_rmag, ec_imag, ec_zmag, ec_ymag, e_ec_gmag, e_ec_rmag, e_ec_imag, e_ec_zmag, e_ec_ymag, de_reddened_gr, de_reddened_gi, de_reddened_ri, de_reddened_gy, de_reddened_gz, de_reddened_ry, de_reddened_rz, de_reddened_iy, de_reddened_iz, de_reddened_zy, e_ec_gr, e_ec_gi, e_ec_gz, e_ec_gy, e_ec_ri, e_ec_rz, e_ec_ry, e_ec_iz, e_ec_iy, e_ec_zy
    return psf_phot
