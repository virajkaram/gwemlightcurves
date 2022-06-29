import os
os.environ['PATH'] = '/scr2/viraj/gwemlightcurves/input/Monica'
from gwemlightcurves.EjectaFits.BNS_fits import calc_meje_dyn, calc_meje_wind
from gwemlightcurves.KNModels import KNTable
import matplotlib.pyplot as plt
from astropy.table import Column

ls = glob('../input/bns_samples_*.dat')
t = KNTable.read_samples('../input/bns_samples_realistic.dat',Nsamples=200)
t = t.calc_radius(EOS='ap4', TOV='Monica', polytrope=True)
t = t.calc_compactness()
t = t.calc_baryonic_mass(EOS='ap4',TOV='Monica')

mej_dyn = calc_meje_dyn(t['m1'],t['c1'],t['m2'],t['c2'])
mej_wind = calc_meje_wind(t['m1'],t['c1'],t['m2'],t['c2'])

t.add_column(Column(name='mej_dyn',data=mej_dyn))
t.add_column(Column(name='mej_wind',data=mej_wind))

t.write('../input/bns_samples_realistic_ejecta_masses.dat',format='ascii')
plt.plot(t['q'],mej_dyn,'.',label='Mej_dyn')
plt.plot(t['q'],mej_wind,'.',label='Mej_wind')
plt.xlabel('q')
plt.ylabel('M')
plt.legend(fontsize=12)
plt.savefig('../input/ejecta_masses.pdf',bbox_inches='tight')
plt.savefig('../input/ejecta_masses.png',bbox_inches='tight')
