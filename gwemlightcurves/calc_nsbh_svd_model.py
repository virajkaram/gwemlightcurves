# This creates the svdmodels from a set of binry model files. The model files need to be located in
# ../output/<model_dir>, where <model_dir> can be looked up in calc_svd_models function of svd_utils.py
import sys
sys.path.append('/scr2/viraj/gwemlightcurves')
from gwemlightcurves import lightcurve_utils, Global, svd_utils
import os, sys, glob, pickle

svd_model = svd_utils.calc_svd_mag(tini=0.01, tmax=20, dt=0.02, n_coeff=100, model='Bu2019nsbh')

with open('../output/svdmodels/Bu2019nsbh_mag.pkl', 'wb') as pklfile:
    pickle.dump(svd_model, pklfile)
