import numpy as np

def read_sam_file(self):
    if self.use_sam is not None:
        if self.use_sam == 'Kurucz':
            print('Reading Interpolated Kurucz SAMs')
            p2 = np.genfromtxt('interpolated_kurucz.txt')
        elif self.use_sam == 'Phoenix':
            print('Reading Interpolated Phoenix SAMs')
            p2 = np.genfromtxt('interpolated_phoenix.txt')

        teff = p2[:,0]
        logg = p2[:,2]
        feh = p2[:,1]
        sam_g = p2[:,3]
        sam_r = p2[:,4]
        sam_i = p2[:,5]
        sam_z = p2[:,6]
        sam_y = p2[:,7]
        sam_j = p2[:,8]
        sam_h = p2[:,9]
        sam_k = p2[:,10]

        sam_params =  teff, logg, feh, sam_g, sam_r, sam_i, sam_z, sam_y, sam_j, sam_h, sam_k
        return sam_params
    else:
        return ('tellar Atmospheric Model files not found')

def select_kurucz_models(self):#, use_interpolated_kurucz_models_greater_than_4000K = False, use_interpolated_kurucz_models_greater_than_4000K = False):
    teff, logg, feh, kg, kr, ki, kz, ky, kj, kh, kk = self.read_sam_file()
    if self.use_interpolated_kurucz_models_greater_than_4000K is True:
            indp = np.where(teff>4000)[0]
            teff_g4k, logg_g4k, feh_g4k, kg_g4k, kr_g4k, ki_g4k, kz_g4k, ky_g4k, kj_g4k, kh_g4k, kk_g4k =  teff[indp], logg[indp], feh[indp], kg[indp], kr[indp], ki[indp], kz[indp], ky[indp], kj[indp], kh[indp], kk[indp]
            sam_params = teff_g4k, logg_g4k, feh_g4k, kg_g4k, kr_g4k, ki_g4k, kz_g4k, ky_g4k, kj_g4k, kh_g4k, kk_g4k
            #return teff_g4k, logg_g4k, feh_g4k, kg_g4k, kr_g4k, ki_g4k, kz_g4k, ky_g4k, kj_g4k, kh_g4k, kk_g4k
            return sam_params

    elif use_interpolated_kurucz_models_lesser_than_4000K is True:
            indp = np.where(teff<4000)[0]
            teff_l4k, logg_l4k, feh_l4k, kg_l4k, kr_l4k, ki_l4k, kz_l4k, ky_l4k, kj_l4k, kh_l4k, kk_l4k =  teff[indp], logg[indp], feh[indp], kg[indp], kr[indp], ki[indp], kz[indp], ky[indp], kj[indp], kh[indp], kk[indp]
            sam_params = teff_l4k, logg_l4k, feh_l4k, kg_l4k, kr_l4k, ki_l4k, kz_l4k, ky_l4k, kj_l4k, kh_l4k, kk_l4k
            #return teff_l4k, logg_l4k, feh_l4k, kg_l4k, kr_l4k, ki_l4k, kz_l4k, ky_l4k, kj_l4k, kh_l4k, kk_l4k
            return sam_params
    else:
            return sam_params

def select_phoenix_models(self):#, use_interpolated_phoenix_models_greater_than_4000K = False, use_interpolated_phoenix_models_greater_than_4000K = False):
    teff,logg,feh,kg,kr,ki,kz,ky,kj,kh,kk = self.read_sam_file(use_phoenix='Phoenix')

    if self.use_interpolated_phoenix_models_greater_than_4000K is True:
            indp = np.where(teff>4000)[0]
            teff_g4k, logg_g4k, feh_g4k, kg_g4k, kr_g4k, ki_g4k, kz_g4k, ky_g4k, kj_g4k, kh_g4k, kk_g4k =  teff[indp], logg[indp], feh[indp], kg[indp], kr[indp], ki[indp], kz[indp], ky[indp], kj[indp], kh[indp], kk[indp]
            return teff_g4k, logg_g4k, feh_g4k, kg_g4k, kr_g4k, ki_g4k, kz_g4k, ky_g4k, kj_g4k, kh_g4k, kk_g4k
        
    elif self.use_interpolated_phoenix_models_lesser_than_4000K is True:
            indp = np.where(teff<4000)[0]
            teff_l4k, logg_glk, feh_glk, kg_glk, kr_glk, ki_l4k, kz_l4k, ky_l4k, kj_l4k, kh_l4k, kk_l4k =  teff[indp], logg[indp], feh[indp], kg[indp], kr[indp], ki[indp], kz[indp], ky[indp], kj[indp], kh[indp], kk[indp]
            return teff_l4k, logg_l4k, feh_l4k, kg_l4k, kr_l4k, ki_l4k, kz_l4k, ky_l4k, kj_l4k, kh_l4k, kk_l4k
        
    else:
            return teff, logg, feh, kg, kr, ki, kz, ky, kj, kh, kk