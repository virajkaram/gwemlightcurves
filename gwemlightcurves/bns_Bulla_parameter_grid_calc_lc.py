from gwemlightcurves.KNModels import KNTable
import numpy as np
import pickle

t = KNTable.read_samples('input/bns_ejecta_parameters_Andreoni.csv',Nsamples=10000) 


t['tini'] = 0.1; t['dt']=0.1;t['tmax']=15

filt = 'J'
model_table = KNTable.model('Bu2019lm',t,LoadModel=True,ModelPath='output/svdmodels/',filterlist=[filt])

pickle.dump(model_table,open('output/bns_Bulla_parameter_grid_Andreoni_%sband.dat'%(filt),'wb'))
