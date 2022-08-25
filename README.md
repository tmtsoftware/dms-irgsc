# dms-irgsc
IR Guide Star Catalog
####################IRGSC##########################

1. Introduction: 
This piece of code generates a catalogue of the NIR magnitudes
of the stars across the TMT's observable sky. It uses the optical data
from the PANSTARRS and computes the NIR magnitudes of the stars by
using Kurucz and Phoenix/NextGen stellar atmospheric models.
The computed magnitudes are validated by using the UKIDSS data which has
observed NIR magnitudes of the stars in the same field. Details of the methods can be found in Shah et.al.


2. Synthetic Photometry:
A catalogue of synthetic magnitudes using the stellar atmospheric models
that are convolved with the PANSTARRS and UKIDSS response functions
are generated and supplied with this code. Kurucz and Phoenix stellar atmospheric models are used to generate this catalogue. The catalogue is then interpolated for magnitudes along the T_eff, logg and [Fe/H] by reducing the spacing between them.
The synthetic photometry files are 'interpolated_kurucz.txt' and 'interpolated_phoenix.txt'. User can select which ever file to use for modeling in the script file.

3. Preparing the data:
The files read_data.py, _sgc.py, extinction_correction.py and _sam.py read the optical and NIR data, seperate the stars and galaxies using the (psf-kron) criteria, correct the stars for reddening and extinction and read the stellar atmospheric models (sam) respectively. The following things to be noted here are:
    a. The PS1 photometry downloaded contains psf and kron magnitudes in all the five optical filters for the sources.
    b. The nan value is denoted by -999 in the PS1 dataset and -9.99999500e+08 in the UKIDSS dataset.
    c. The user can select whether he wants to apply the whole set of Kurucz/Phoenix sams or only models that have Teff>4000K/Teff<4000K. This has to be specified in the script.

4. Fitting of models to the data:
In the _fitting.py code, there are various functions to fit the models as describe in Shah et.al.:
a. reduced_chi2_for_all_stars: In this function, the reduced chi2 for all the stars is calculated by matching the observed colours to the model colours.
b. reduced_chi2_leq_2: In the function only stars that have reduced chi2 value leq 2 are selected to generate the catalogue. This function also uses the colour matching method.
c. compute_nir2: In this function, the reduced chi2 is computed for all the stars by keeping the extinction and scale factor free and minimizing the chi2. Raw i.e. optical photometry that is not corrected for reddening is supplied to this function.

The user can specify which method to use in the scipt file. The computed NIR magnitudes are converted to Vega system from AB system in these functions itself.

5. Validation:
In this file, the computed NIR magnitudes are valiated in the validation_ukidss function. This can be done only when the UKIDSS data is supplied to the code. The function matches the position of the stars with the sources in the UKIDSS data to 1" and the finds the accuracy of the computed NIR magnitudes. In the plot_validation_plots function, the scatter plots of the difference in the observed and computed NIR magnitudes vs the observed NIR magnitudes are then plotted. A user can specify in the script file if he/she wants to validate the computed magnitudes.
