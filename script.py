import os
import sys
import numpy as np
import matplotlib

sys.path.append('/mnt/c/Users/sshah/Documents/itcc/irgsc/irgsc_class/') # change this to the location where you have installed IRGSC

from __init__ import IRGSC


ps1 = np.genfromtxt('ps_tf6.csv', delimiter = ',') # optical data from PANSTARRS
ukidss = np.genfromtxt('tf6.fits') # NIR observed data from UKIDSS to validate the computed NIR magnitudes

ic = IRGSC(optical_data = ps1)

ic.validate=True #False if there is no data to validate
ic.use_kurucz = True #Use Kurucz stellar atmospheric models
ic.use_phoenix = True #use Phoenix stellar atmospheric models. User can only use Kurucz models if NIR magnitudes in H and Ks bands need to be more accurate

self.validate = False

self.aj = 0.033 #extinction coefficient in J band
self.ah = 0.025 #extinction coefficient in H band
self.ak = 0.017 #extinction coefficient in Ks band
self.ebv = 0.04 #mean value of reddening
self.e_ebv = 0.00016 #sigma reddening
#extinctions in the optical bands are calculated using the relations from Tonry et.al. 2012
#in the _extinction_correction.py file

ic.use_interpolated_kurucz_models_greater_than_4000k = True

ic.reduced_chi2_all_stars() #outputs catalogue.txt file which contains computed NIR magnitudes for all the stars. See _fitting.py file for more details.
#ic.reduced_chi2_leq_2() #outputs catalogue.txt file which contains computed NIR magnitudes for stars having chi2_r <=2.0
#use self.use_compute_nir_by_keeping_sf_and_reddening_free is True to use ic.compute_nir2()
#ic.compute_nir2() #outputs catalogue.txt file which contains the computed NIR magnitudes for stars by keeping scale factor flux and reddening as free parameters.
#ic.validation_ukidss()
#ic.plot_validation_plots()
