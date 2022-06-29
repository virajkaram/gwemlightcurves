from astropy.table import Table
from gwemlightcurves.KNModels import KNTable
import numpy as np
from astropy.table import Column

#Mej_dyn = np.linspace(0.001,0.02,10)
#Mej_wind = np.linspace(0.01,0.13,10)
Mej_dyn = np.linspace(0.01,0.09,10)
Mej_wind = np.linspace(0.01,0.09,10)

phis = np.array([10,30,45,60,75,90])#np.linspace(30,60,5)
costhetas = np.linspace(0,1,5)
thetas = np.arccos(costhetas)*(180/np.pi)

#thetas = np.array([ 4.01400815,  5.66079142,  7.77481813,  8.06401156, 16.99208691,
#       17.00978436, 18.00516238, 21.24470016, 23.40579027, 24.0112735 ,
#       25.28370486, 25.83485555, 26.32700554, 26.80824358, 27.47679023,
#       29.15216991, 29.33288527, 29.87260256, 30.44672259, 31.40503272,
#       31.9187994 , 33.05273705, 34.1613424 , 34.62285375, 35.43257655,
#       35.89145317, 36.16977171, 39.2555548 , 42.2579595 , 44.04547949,
#       44.26885816, 44.39793328, 45.74930506, 45.89198028, 45.91909639,
#       47.61342572, 48.50119428, 50.04437586, 51.3506148 , 55.61290305,
#       55.82400947, 61.03718753, 61.86678952, 62.80751381, 65.87598187,
#       67.88600199, 68.70225547, 69.33249482, 69.36173921, 70.52121799,
#       70.89863086, 71.58184137, 75.0208062 , 81.03504812])

t = Table(names=['mej_dyn','mej_wind','phi','theta'])
for mdyn in Mej_dyn:
    for mwind in Mej_wind:
        for theta in thetas:
            for phi in phis:
                t.add_row([mdyn,mwind,phi,theta])

m1 = np.ones(len(t))
m2 = np.ones(len(t))
q = m1/m2

t.add_column(Column(name='m1',data=m1))
t.add_column(Column(name='m2',data=m2))
t.add_column(Column(name='q',data=q))

t.write('../output/nsbh_ejecta_parameters_Andreoni.csv',overwrite=True)
