from gwemlightcurves.KNModels import KNTable

#t = KNTable.read_samples('input/bns_samples_realistic_test_Bu2019lm.dat')

#t['tini'] = 1; t['dt']=0.1;t['tmax']=10 

model_table = KNTable.model('Bu2019lm',t,LoadModel=True,ModelPath='output/svdmodels/',filterlist=['g'])

t['tini'] = 1; t['dt']=0.1;t['tmax']=10
model_table = KNTable.model('Bu2019lm',t,LoadModel=True,ModelPath='output/svdmodels/',filterlist=['g'])

print(model_table)