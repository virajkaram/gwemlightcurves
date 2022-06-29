from gwemlightcurves.KNModels import KNTable
import pickle
import numpy as np

filt = 'r'
t = KNTable.read_samples('../input/bns_samples_realistic_ejecta_masses.dat',Nsamples=20)
t['tini'] = 0.01; t['dt']=0.02;t['tmax']=7
t['phi'] = 45.0
#thetas = [0,60,90]
'''
for theta in thetas:
    print('Currently doing angle',theta)
    t['theta'] = theta
    model_table1 = KNTable.model('Bu2019lm',t,LoadModel=True,ModelPath='output/svdmodels/',filterlist=[filt])

    pickle.dump(model_table1,open('output/bns_samples_realistic_ejecta_masses_%sband_theta%s.dat'%(filt,theta),'wb'))
'''
t['theta'] = t['theta_jn']*180/np.pi

model_table = KNTable.model('Bu2019lm',t,LoadModel=True,ModelPath='../output/svdmodels/',filterlist=[filt])
pickle.dump(model_table,open('../output/bns_samples_realistic_Bulla_%sband.pickle'%(filt),'wb'))

