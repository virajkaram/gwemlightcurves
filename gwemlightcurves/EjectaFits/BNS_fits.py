import numpy as np
import scipy
#Viraj's notes : These fitting formulae are the ones stated in Stachie et al. 2021. They are also derived in various other papers including Coughlin et al. 2019a, Dietrich et al. 2020

def calc_meje_dyn(m1,c1,m2,c2):
    a= -0.0719
    b= 0.2116
    d= -2.42
    n= -2.905

    log10_mej = a*(m1*(1-2*c1)/c1 + m2*(1-2*c2)/c2) + b*(m1*(m2/m1)**n + m2*(m1/m2)**n)+d
    meje_dynamical_fit = 10**log10_mej     

    return meje_dynamical_fit


def calc_meje_wind(m1,c1,m2,c2):
    lambda_coeff = np.array([374839, -1.06499e7, 1.27306e8, -8.14721e8, 2.93183e9, -5.60839e9, 4.44638e9])
    coeff = lambda_coeff[::-1]
    p = np.poly1d(coeff)
    lambda1 = p(c1)
    lambda2 = p(c2)
    lambda1[lambda1 < 0] = 0
    lambda2[lambda2 < 0] = 0
    q = m1/m2

    lambdatilde = (16.0/13.0)*(lambda2 + lambda1*(q**5) + 12*lambda1*(q**4) + 12*lambda2*q)/((q+1)**5)
    mc = ((m1*m2)**(3./5.)) * ((m1 + m2)**(-1./5.))

    mTOV = 2.17
    R16 = mc * (lambdatilde/0.0042)**(1.0/6.0)
    rat = mTOV/R16
    mth = (2.38 - 3.606*mTOV/R16)*mTOV

    a0 = -1.581
    da = -2.439
    b0 = -0.538
    db = -0.406
    c = 0.953
    d = 0.0417
    beta = 3.910
    qtrans = 0.900

    x = lambdatilde*1.0
    mtot = m1+m2

    eta = 0.5 * np.tanh(beta*(q-qtrans))
    a = a0 + da * eta
    b = b0 + db * eta

    mdisk = a*(1+b*np.tanh((c-mtot/mth)/d))

    mdisk[mdisk<-3] = -3.0
    mdisk[rat>0.32] = -3.0
    mdisk = 10**mdisk
 
    zeta = 0.15
    meje_wind_fit = zeta*mdisk
    meje_wind_fit[meje_wind_fit > 0.1] = 0.1

    return meje_wind_fit


def calc_vej(m1,c1,m2,c2):
    """
.. py:function:: calc_vrho(m1,c1,m2,c2)

    velocity mass ejecta

    https://arxiv.org/pdf/1612.03665.pdf#equation.3.5

    https://arxiv.org/pdf/1612.03665.pdf#equation.3.6

    a = −0.219479,  b= 0.444836,  c=−2.67385

   :param float m1: mass of larger ns (MSun)
   :param float c1: compactness of the larger neutron star
   :param float m2: mass of samller ns (MSun)
   :param float c2: compactness of the smaller neutron star
   :return: velocity of ejecta mass (Msun)
   :rtype: float
    """
    a=-0.3090
    b=0.657
    c=-1.879

    return a*(m1/m2)*(1+c*c1) + a*(m2/m1)*(1+c*c2)+b