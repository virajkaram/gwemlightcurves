import os
os.environ['PATH'] = '/scr2/viraj/gwemlightcurves/input/Monica'
from gwemlightcurves.EjectaFits.BNS_fits import calc_meje_dyn, calc_meje_wind
from gwemlightcurves.KNModels import KNTable
import matplotlib.pyplot as plt
from astropy.table import Column


ls = ['bns_samples_optimistic.dat','bns_samples_optimistic_catalog.dat','bns_samples_pessimistic.dat','bns_samples_pessimistic_catalog.dat','bns_samples_realistic.dat','bns_samples_realistic_catalog.dat']

for l in ls:

	t = KNTable.read_samples('../input/%s'%(l),Nsamples=1000)
	t = t.calc_radius(EOS='ap4', TOV='Monica', polytrope=True)
	t = t.calc_compactness()
	t = t.calc_baryonic_mass(EOS='ap4',TOV='Monica')

	mej_dyn = calc_meje_dyn(t['m1'],t['c1'],t['m2'],t['c2'])
	mej_wind = calc_meje_wind(t['m1'],t['c1'],t['m2'],t['c2'])

	t.add_column(Column(name='mej_dyn',data=mej_dyn))
	t.add_column(Column(name='mej_wind',data=mej_wind))

	t.write('../input/%s_ejecta_masses.dat'%(l.split('.dat')[0]),format='ascii')
	
	plt.plot(t['q'],mej_dyn,'.',label='Mej_dyn')
	plt.plot(t['q'],mej_wind,'.',label='Mej_wind')
	plt.xlabel('q')
	plt.ylabel('M')
	plt.legend(fontsize=12)
	plt.savefig('../input/%s_ejecta_masses.pdf'%(l.split('.dat')[0]),bbox_inches='tight')
	plt.savefig('../input/%s_ejecta_masses.png'%(l.split('.dat')[0]),bbox_inches='tight')
