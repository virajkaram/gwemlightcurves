from gwemlightcurves.KNModels import KNTable
import numpy as np
import pickle
import sys
from datetime import datetime

t_start = datetime.utcnow()
print('Started at',t_start)

filt = str(sys.argv[1])
num = str(sys.argv[2])
t = KNTable.read_samples('../input/bns_ejecta_parameters_Andreoni_large_ranges%s.csv'%(num),Nsamples=100000) 


t['tini'] = 0.01; t['dt']=0.1;t['tmax']=15

#filt = str(sys.argv[1])
print('Doing filter',filt)
model_table = KNTable.model('Bu2019lm',t,LoadModel=True,ModelPath='../output/svdmodels/',filterlist=[filt])

pickle.dump(model_table,open('../output/bns_Bulla_parameter_grid_Andreoni_large_ranges%s_%sband.dat'%(num,filt),'wb'))

t_end = datetime.utcnow()
print('Ended at',t_end)
print('Elapsed time',t_end-t_start)
