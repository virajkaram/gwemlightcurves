from __future__ import division
import matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('xtick', labelsize=16)
matplotlib.rc('ytick', labelsize=16)
matplotlib.rc('font', **font)
from argparse import ArgumentParser
import numpy as np
import astropy.units as u
from astropy.table import (Table, Column, vstack)
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from astropy.time import Time
from astropy.io import ascii
import pickle
import argparse
import sys
import os
import requests
import glob
import scipy.stats as ss
from scipy.stats import rv_continuous
from scipy.integrate import quad
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from pandas import read_csv, DataFrame
import healpy as hp
from ligo.skymap import postprocess, distance
from ligo.skymap.io import fits
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import pandas as pd



### non-standard libraries

from gwemlightcurves.KNModels import KNTable
from gwemlightcurves import __version__
from gwemlightcurves.EOS.EOS4ParameterPiecewisePolytrope import EOS4ParameterPiecewisePolytrope
#from twixie import kde
from gwemlightcurves import lightcurve_utils
from mass_grid import run_EOS
#

fig, ax = plt.subplots(figsize=(16, 12))

mej_theta_data=np.loadtxt('mej_theta_data_BNS_alsing.txt')


mej_data, thetas = mej_theta_data[:,0], mej_theta_data[:,1] 
l= len(mej_data)

phis = 30+30*np.random.rand(l)
samples = Table((mej_data, phis, thetas), names=('mej', 'phi', 'theta'))

tini = 0.1
tmax = 50.0
dt = 0.1

vmin = 0.02
th = 0.2
ph = 3.14
kappa = 10.0
eps = 1.58*(10**10)
alp = 1.2
eth = 0.5
flgbct = 1

beta = 3.0
kappa_r = 0.1
slope_r = -1.2
theta_r = 0.0
Ye = 0.3

samples['tini'] = tini
samples['tmax'] = tmax
samples['dt'] = dt
samples['vmin'] = vmin
samples['th'] = th
samples['ph'] = ph
samples['kappa'] = kappa
samples['eps'] = eps
samples['alp'] = alp
samples['eth'] = eth
samples['flgbct'] = flgbct
samples['beta'] = beta
samples['kappa_r'] = kappa_r
samples['slope_r'] = slope_r
samples['theta_r'] = theta_r
samples['Ye'] = Ye


ModelPath = "/home/cosmin.stachie/gwemlightcurves/output/svdmodels"
kwargs = {'SaveModel':False,'LoadModel':True,'ModelPath':ModelPath}
kwargs["doAB"] = True
kwargs["doSpec"] = False

#model = "Ka2017"
model = "Bu2019inc"
model_tables = {}
model_tables[model] = KNTable.model(model, samples, **kwargs)

idx = np.where(model_tables[model]['mej'] <= 1e-3)[0]
print("idx")
print(idx)
model_tables[model]['mag'][idx] = 10.
model_tables[model]['lbol'][idx] = 1e30
 
mags = model_tables[model]['mag'][0]
t = model_tables[model]['t'][0]
#print(np.shape(t))
for i, band in enumerate(mags):
    l = len(band)
    plt.plot(t, band)
plt.ylim([10,-20])

plt.savefig('mag_test.pdf')

#ax.plot([.03,.03],[1e-6,1e6], color = 'black', linestyle='--')
#ax.plot([.05,.05],[1e-6,1e6], color = 'black', linestyle='--')
#ax.set_yscale('log')
#ax.set_xscale('log')
#ax.set_xlim(-9.8, -20.2)
#plt.ylim(1e-3,1.1)
plt.xlabel('Magnitude')
plt.ylabel('Probability')
plt.legend()
plt.grid()
plt.savefig('mag_test.pdf')

plt.close(fig) 
