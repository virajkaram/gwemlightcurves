from gwemlightcurves.KNModels import KNTable
import pickle
import numpy as np
from gwemlightcurves.EjectaFits.BNS_fits import calc_meje_dyn, calc_meje_wind, calc_vej


t = KNTable.read_samples('input/bns_samples_realistic.dat',Nsamples=200)
t = t.calc_radius(EOS='ap4', TOV='Monica', polytrope=True)
t = t.calc_compactness()
t = t.calc_baryonic_mass(EOS='ap4',TOV='Monica')

mej_dyn = calc_meje_dyn(t['m1'],t['c1'],t['m2'],t['c2'])
mej_wind = calc_meje_wind(t['m1'],t['c1'],t['m2'],t['c2'])
t['mej'] = t['mej_dyn'] + t['mej_wind']
t['vej'] = calc_vej(t['m1'],t['c1'],t['m2'],t['c2'])

t['tini'] = 0.1; t['dt']=0.1;t['tmax']=15


Xlans = [1e-5,3e-5,6e-5,1e-4] #Similar to Andreoni et al. 2021
filt = 'g'
for Xlan in Xlans:
    print('Currently doing lanthanide fraction',Xlan)
    t['Xlan'] = Xlan
    model_table1 = KNTable.model('Ka2017',t,LoadModel=True,ModelPath='output/svdmodels/',filterlist=[filt])

    pickle.dump(model_table1,open('output/bns_samples_realistic_ejecta_masses_Ka2017_%sband_Xlan%s.dat'%(filt,theta),'wb'))
