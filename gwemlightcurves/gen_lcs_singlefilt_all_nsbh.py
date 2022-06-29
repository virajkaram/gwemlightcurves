from gwemlightcurves.KNModels import KNTable
import pickle
import numpy as np
import glob
import sys

filt = str(sys.argv[1])
phi = float(sys.argv[2])
print(filt,phi)

ls = ['bns_samples_optimistic.dat','bns_samples_optimistic_catalog.dat','bns_samples_pessimistic.dat','bns_samples_pessimistic_catalog.dat','bns_samples_realistic.dat','bns_samples_realistic_catalog.dat']

for l in ls:
	print('Now doing l')
	t = KNTable.read_samples('../input/%s_ejecta_masses.dat'%(l.split('.dat')[0]),Nsamples=10)
	t['tini'] = 0.01; t['dt']=0.02;t['tmax']=15
	t['phi'] = phi
	#thetas = [0,60,90]
	'''
	for theta in thetas:
    	print('Currently doing angle',theta)
    	t['theta'] = theta
    	model_table1 = KNTable.model('Bu2019lm',t,LoadModel=True,ModelPath='output/svdmodels/',filterlist=[filt])

    	pickle.dump(model_table1,open('output/bns_samples_realistic_ejecta_masses_%sband_theta%s.dat'%(filt,theta),'wb'))
	'''
	t['theta'] = np.arccos(np.abs(np.cos(t['theta_jn'])))*180/np.pi

	model_table = KNTable.model('Bu2019nsbh',t,LoadModel=True,ModelPath='/scr2/viraj/gwemlightcurves/bin/',filterlist=[filt])
	pickle.dump(model_table,open('../output/nsbh_%sband_phi%s_abscos.pickle'%(l.split('.dat')[0],filt,phi),'wb'))
