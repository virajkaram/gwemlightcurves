from gwemlightcurves.KNModels import KNTable
import numpy as np
import pickle
import sys
from datetime import datetime

t_start = datetime.utcnow()
print('Started at',t_start)

filt = str(sys.argv[1])
#num = str(sys.argv[2])
t = KNTable.read_samples('../input/nsbh_ejecta_parameters_Andreoni.csv',Nsamples=100000)

t['tini'] = 0.01; t['dt']=0.1;t['tmax']=15

#filt = str(sys.argv[1])
print('Doing filter',filt)
model_table = KNTable.model('Bu2019nsbh',t,LoadModel=True,ModelPath='../output/svdmodels/',filterlist=[filt])

pickle.dump(model_table,open('../output/nsbh_Bulla_parameter_grid_Andreoni_%sband.dat'%(filt),'wb'))

t_end = datetime.utcnow()
print('Ended at',t_end)
print('Elapsed time',t_end-t_start)
