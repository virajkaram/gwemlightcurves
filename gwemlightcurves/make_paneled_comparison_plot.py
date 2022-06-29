#!/usr/bin/env python
# coding: utf-8

# In[3]:


from astropy.io import ascii
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from astropy.table import vstack
from glob import glob
from astropy.table import Table
import pickle
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib import colors
import matplotlib.patches as mpatches
import pandas as pd
import matplotlib
import os

def init():
    matplotlib.rcParams['xtick.minor.size'] = 6
    matplotlib.rcParams['xtick.major.size'] = 6
    matplotlib.rcParams['ytick.major.size'] = 6
    matplotlib.rcParams['ytick.minor.size'] = 6
    matplotlib.rcParams['lines.linewidth'] = 1.5
    matplotlib.rcParams['axes.linewidth'] = 1.5
    matplotlib.rcParams['font.size']= 16
    matplotlib.rcParams['font.family']= 'sans-serif'
    matplotlib.rcParams['xtick.major.width']= 2.
    matplotlib.rcParams['ytick.major.width']= 2.
    matplotlib.rcParams['ytick.direction']='in'
    matplotlib.rcParams['xtick.direction']='in'

#init()
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['text.usetex'] = True

# In[2]:


filt = 'r'
with open('../output/bns_Bulla_parameter_grid_Andreoni_%sband.dat'%(filt),'rb') as f:
    r = pickle.load(f)
    
filt = 'J'    
with open('../output/bns_Bulla_parameter_grid_Andreoni_%sband.dat'%(filt),'rb') as f:
    J = pickle.load(f)
    
peakmags = []
for row in r:
    peakmags.append(np.min(row['mag'][0]))
r['peak_mag'] = peakmags

peakmags = []
for row in J:
    peakmags.append(np.min(row['mag'][0]))
J['peak_mag'] = peakmags


'''
# In[3]:


filt = 'r'
with open('output/bns_Bulla_parameter_grid_Andreoni_%sband.dat'%(filt),'rb') as f:
    r = pickle.load(f)
    
filt = 'J'    
with open('output/bns_Bulla_parameter_grid_Andreoni_%sband.dat'%(filt),'rb') as f:
    J = pickle.load(f)
''' 


# In[4]:


thetas = np.unique(J['theta'])
phis = np.unique(J['phi'])
mej_dyns = np.unique(J['mej_dyn'])
mej_winds = np.unique(J['mej_wind'])


# In[6]:


#ls = glob('/Users/viraj/winter/gwemopt_sims/output_parallel/realistic_isPE1_isGal0_*')
sched = Table.read('../input/realistic_isPE1_isGal0_isGalTiling0_3days.csv')

samples = pd.read_csv('../input/bns_samples_realistic.dat',delimiter='\t')
costhetas = []
distances = []
for row in sched:
    if row['distance'] in distances:
        continue
    costhetas.append(np.abs(np.cos(samples.iloc[int(row['event'])]['theta_jn'])))
    distances.append(round(row['distance'],3))
    print(row['distance'],np.abs(np.cos(samples.iloc[int(row['event'])]['theta_jn'])))

sinds = np.argsort(distances)

costhetas = np.array(costhetas)
distances = np.array(distances)

realistic_costhetas = costhetas[sinds]
realistic_distances = distances[sinds]
print(realistic_costhetas,realistic_distances)

mask = (realistic_distances<400)
realistic_distances = realistic_distances[mask]
realistic_costhetas = realistic_costhetas[mask]
#area90s = area90s[mask]
#realistic_costhetas = np.unique(costhetas)
#realistic_distances = np.unique(distances)

#ls = glob('/Users/viraj/winter/gwemopt_sims/output_parallel/realistic_isPE1_isGal0_*')
#sched = Table.read(ls[2])
for row in sched:
    print(row['area90']/(8*3600/450), row['distance'])


# In[9]:


cmap = colors.LinearSegmentedColormap.from_list("",['ghostwhite','cornflowerblue','darkblue'])
color_list = [plt.get_cmap('coolwarm')(i) for i in np.linspace(0, 1, 10)]
tablelist = [r,J]
cols = ['red','orange']
labs = ['r','J']

fig = plt.figure(figsize=(25,15))
gs = GridSpec(2,2,hspace=0.5,wspace=0.2)

axs = []
#is_dets = [[],[]]

#is_dets = [[],[]]
gs1 = gs[0].subgridspec(1,2,wspace=0.1)


phi_plots = [phis[0],phis[-1]]

delim_angle = 60
theta_plot_lims = [(0,delim_angle),(delim_angle,90)]

enum = -1
kn_types = [['Blue on-axis','Blue off-axis'],['Red on-axis','Red off-axis']]
for phi_ind,phi_val in enumerate(phi_plots):
    for theta_ind,theta_lim in enumerate(theta_plot_lims):
        enum = enum+1
        pmarrays = []
        
        max_dms = [[],[]]
        gs1 = gs[enum].subgridspec(1,2,wspace=0.1)
        
        for ind in range(len(tablelist)):
            t = tablelist[ind]
            phi = phi_val
            t1 = t[(t['theta']>=theta_lim[0]) & (t['theta']<theta_lim[1]) & (t['phi']==phi)]
            dmarray = []
            for md in mej_dyns:
                t2 = t1[t1['mej_dyn']==md]
                t2.sort('mej_wind')
                t2 = t2.group_by('mej_wind')
                dmbins = []
                for grp in t2.groups:
                    dmbins.append(21-np.median(grp['peak_mag']))
                dmarray.append(dmbins)
            #plt.plot(t1['theta'],t1['peak_mag'],'.',c=cols[ind])
            dmarray = np.array(dmarray)
            #ax = fig.add_subplot(gs[ind])
            #is_det = (pmarray<21)
            #is_dets[ind].append(is_det)
            max_dms[ind].append(dmarray)
        
        max_dms = np.array(max_dms)
        max_dists = 10*(10**(max_dms/5))/1e6
        
        if theta_lim[1]==delim_angle:
            dists_to_plot = realistic_distances[(realistic_costhetas>np.cos(delim_angle*np.pi/180))]
            print('Blue r',len(dists_to_plot),len(dists_to_plot[dists_to_plot>np.max(max_dists[0][0])]))
            print('Blue GW170817 r',len(dists_to_plot),len(dists_to_plot[dists_to_plot<max_dists[0][0][3][5]]))
            print('Blue J',len(dists_to_plot),len(dists_to_plot[dists_to_plot>np.max(max_dists[1][0])]))
            print('Blue GW170817 J',len(dists_to_plot),len(dists_to_plot[dists_to_plot<max_dists[1][0][3][5]]))
        else:
            dists_to_plot = realistic_distances[(realistic_costhetas<np.cos(delim_angle*np.pi/180))]
            print('Red r',len(dists_to_plot),len(dists_to_plot[dists_to_plot>np.max(max_dists[0][0])]))
            print('Red J',len(dists_to_plot),len(dists_to_plot[dists_to_plot>np.max(max_dists[1][0])]))
            print('Red GW170817 r',len(dists_to_plot),len(dists_to_plot[dists_to_plot<max_dists[0][0][3][5]]))
            print('Red GW170817 J',len(dists_to_plot),len(dists_to_plot[dists_to_plot<max_dists[1][0][3][5]]))
        ax = fig.add_subplot(gs1[0])
        im = ax.imshow(max_dists[0][0],extent=[mej_winds.min(),mej_winds.max(),mej_dyns.min(),mej_dyns.max()],origin='lower',aspect=7,vmin=150,vmax=500,cmap=cmap,interpolation='bicubic')
        cs = ax.contour(max_dists[0][0],[max_dists[0][0][3][5]],colors=[color_list[-1]],extent=[mej_winds.min(),mej_winds.max(),mej_dyns.min(),mej_dyns.max()],origin='lower',linestyles='--',linewidth=3)
        cs1 = ax.contour(max_dists[0][0],levels=dists_to_plot,extent=[mej_winds.min(),mej_winds.max(),mej_dyns.min(),mej_dyns.max()],origin='lower',linestyles='-',alpha=0.6,vmin=150,vmax=500,cmap='magma',linewidth=0.7)
        plt.title('%s-band'%(labs[0]),size=25)
        ax.set_xticks([0.03,0.06,0.09,0.12])
        ax.set_xticklabels([0.03,0.06,0.09,0.12],size=20)
        ax.set_xlabel(r'M$_{\rm{ej, wind}}$ (M$_\odot$)',size=28)
        ax.set_ylabel(r'M$_{\rm{ej, dyn}}$ (M$_\odot$)',size=28)
        ax.set_yticks([0.005,0.01,0.015,0.02])
        ax.set_yticklabels([0.005,0.01,0.015,0.02],size=20)
        ax.tick_params(size=10)
        #cbar = fig.colorbar(im)
        
        ax.plot(mej_winds[5],mej_dyns[3],'*',c=color_list[-1],markersize=20)
        ax.text(0.09,0.0035,'%i Mpc'%(max_dists[0][0][3][5]),size=20,color=color_list[-1])
        
        #fmt = {}
        #strs = ['%i Mpc'%(max_dists[0][0][3][5])]
        #for l, s in zip(cs.levels, strs):
        #    fmt[l] = s
        #ax.clabel(cs,cs.levels, fmt=fmt,  inline=True,fontsize=20)
        
        ax = fig.add_subplot(gs1[1])
        im = ax.imshow(max_dists[1][0],extent=[mej_winds.min(),mej_winds.max(),mej_dyns.min(),mej_dyns.max()],origin='lower',aspect=7,vmin=150,vmax=500,cmap=cmap,interpolation='bicubic')
        cs = ax.contour(max_dists[1][0],[max_dists[1][0][3][5]],colors=[color_list[-1]],extent=[mej_winds.min(),mej_winds.max(),mej_dyns.min(),mej_dyns.max()],origin='lower',linestyles='--',linewidth=3)
        cs1 = ax.contour(max_dists[1][0],levels=dists_to_plot,extent=[mej_winds.min(),mej_winds.max(),mej_dyns.min(),mej_dyns.max()],origin='lower',linestyles='-',alpha=0.6,vmin=150,vmax=500,cmap='magma',linewidth=0.7)
        plt.title('%s-band'%(labs[1]),size=25)
        #ax.set_yticks([0.005,0.01,0.015,0.02])
        ax.set_xticks([0.03,0.06,0.09,0.12])
        ax.set_xticklabels([0.03,0.06,0.09,0.12],size=20)
        ax.set_xlabel(r'M$_{\rm{ej, wind}}$ (M$_\odot$)',size=25)
        ax.set_yticks([0.005,0.01,0.015,0.02])
        ax.set_yticklabels([],size=12)
        ax.tick_params(size=10)
        sign = '>'
        if theta_lim[1] == delim_angle:
            sign = '<'
        ax.text(-0.1,0.022,r'%s :   $\Phi = %i^{\rm{o}}$, $\Theta_{\rm{obs}} %s %s^{\rm{o}}$'%(kn_types[phi_ind][theta_ind],phi,sign,delim_angle),size=25)
        
        ax.plot(mej_winds[5],mej_dyns[3],'*',c=color_list[-1],markersize=20)
        ax.text(0.09,0.0025,'%i Mpc'%(max_dists[1][0][3][5]),size=20,color=color_list[-1])
#patch100 = mpatches.Patch(facecolor='#95b6f3',edgecolor='black',label=r'100 Mpc')
#patch200 = mpatches.Patch(facecolor='#4263cc',edgecolor='black',label=r'200 Mpc')
#patch300 = mpatches.Patch(facecolor='#00008b',edgecolor='black',label=r'300 Mpc')
#fig.subplots_adjust(right=0.9,left=0)
cbar_ax = fig.add_axes([0.93,0.15,0.02,0.7])
cbar = fig.colorbar(im,cax=cbar_ax)
cbar.set_label('Max. detectable distance [Mpc]',size=25)
cbar.set_ticks([150,200,250,300,350,400,450,500])
cbar.set_ticklabels([150,200,250,300,350,400,450,500])
cbar.ax.tick_params(labelsize=20,pad=5)
#plt.legend(handles=[patch100,patch200,patch300],bbox_to_anchor=(1.5,1.4),fontsize=20)
plt.savefig('rJ_comparison_ejecta_masses_Bulla_models.pdf',bbox_inches='tight')


# In[195]:



# In[175]:


a = ascii.read('realistic_localisation_areas.csv')
area90s = []
area50s = []
distances = []
for row in a:
    ind = np.argmin(np.abs(row['Area90'] - sched['area90']))
    area90s.append(row['Area90'])
    area50s.append(row['Area50'])
    distances.append(sched[ind]['distance'])
    #print(row['Area90'],sched[ind]['area90'],sched[ind]['distance'])
area90s = np.array(area90s)
area50s = np.array(area50s)
distances = np.array(distances)


# In[192]:


cmap = colors.LinearSegmentedColormap.from_list("",['ghostwhite','cornflowerblue','darkblue'])
labs = ['r','J']
Ds = np.linspace(20,500,20)
tablelist = [r,J]
DMs = 5*np.log10(Ds*1e5)
t_covs = np.arange(0,15,0.3)

fig = plt.figure(figsize=(14,6))
gs = GridSpec(1,2,hspace=0.4,wspace=0.3)

axs = []
delim_angle = 60
theta_plot_lims = [(0,delim_angle),(delim_angle,90)]
enum =-1
phi_plots = [phis[0],phis[-1]]

enum = enum+1
pmarrays = []
        
#gs1 = gs[enum].subgridspec(1,2,wspace=0.1)
        
det_fracs = [[],[]]
#fdets = [[],[]]
        

for ind in range(len(tablelist)):
    t = tablelist[ind]
    t1 = t


    for DM in DMs:
        fdetarray = []
        for t_cov in t_covs:
            mapps = t1['mag'][:,0]+DM
            t_ind = np.argmin(np.abs(t1[0]['t']-t_cov))
            m_at_cov = mapps[:,t_ind]
            ndet = np.sum(m_at_cov<21)
            fdet = ndet/len(m_at_cov)
                    
            fdetarray.append(fdet)
        det_fracs[ind].append(fdetarray)    
        
    det_fracs = np.array(det_fracs)
        
ax = fig.add_subplot(gs[0])
im = ax.imshow(det_fracs[0],extent=[t_covs.min(),t_covs.max(),Ds.min(),Ds.max()],origin='lower',aspect=0.026,vmin=0,vmax=1,cmap=cmap,interpolation='bicubic')
cs = ax.contour(det_fracs[0],[0.1,0.5,0.9],colors='black',extent=[t_covs.min(),t_covs.max(),Ds.min(),Ds.max()],origin='lower',linestyles=['--','-',':'],linewidth=3)
#ax.set_title('%s-band'%(labs[0]),size=15)
ax.set_ylim(80,450)
ax.set_xlim(0,10)
ax.set_yticks([100,200,300,400])
ax.set_xticks([0,2,4,6,8,10])
#ax.set_xticklabels([0.03,0.06,0.09,0.12],size=20)
ax.set_ylabel(r'Distance [Mpc]',size=16)
ax.set_xlabel(r'Time searched [days]',size=16)
#ax.set_yticks([0.005,0.01,0.015,0.02])
#ax.set_yticklabels([0.005,0.01,0.015,0.02],size=20)

ax.text(4,400,'r-band',size=15)
ax.text(4,380,r'm$_{\rm{lim}} = 21$',size=12)
ax2 = ax.twiny()
ax2.set_xlim(np.array([0,10])*9*3600/450 * 1)
ax2.set_xlabel(r'Max. area searchable with WINTER [deg$^{2}$]',size=14,labelpad=10)

ax.plot(sched['area90']/(9*3600/450),sched['distance'],'x',color=color_list[-1],markersize=7)
#ax.plot(area50s/(9*3600/(450)),distances,'x',color='red',markersize=7)

ax.tick_params(labelsize=12)
ax2.tick_params(labelsize=12)

ax.plot(0,0, label=r'f$_{\rm{det}} = 0.1$', color='black', ls='--')
ax.plot(0,0, label=r'f$_{\rm{det}} = 0.5$', color='black', ls='-')
ax.plot(0,0, label=r'f$_{\rm{det}} = 0.9$', color='black', ls=':')
ax.legend(fontsize=12)
#ax2.set_xticks(xticks*8*3600/450 * 1)
#cbar = fig.colorbar(im)
        
#ax.plot(mej_winds[5],mej_dyns[3],'*',c='red',markersize=20)
#ax.text(0.09,0.0035,'%i Mpc'%(max_dists[0][0][3][5]),size=20,color='red')
        
#fmt = {}
#strs = ['%i Mpc'%(max_dists[0][0][3][5])]
#for l, s in zip(cs.levels, strs):
#    fmt[l] = s
#ax.clabel(cs,cs.levels, fmt=fmt,  inline=True,fontsize=20)
        
ax = fig.add_subplot(gs[1])
im = ax.imshow(det_fracs[1],extent=[t_covs.min(),t_covs.max(),Ds.min(),Ds.max()],origin='lower',aspect=0.026,vmin=0,vmax=1,cmap=cmap,interpolation='bicubic')
cs = ax.contour(det_fracs[1],[0.1,0.5,0.9],colors='black',extent=[t_covs.min(),t_covs.max(),Ds.min(),Ds.max()],origin='lower',linestyles=['--','-',':'],linewidth=3)
#cs1 = ax.contour(max_dists[1][0],levels=dists_to_plot,extent=[mej_winds.min(),mej_winds.max(),mej_dyns.min(),mej_dyns.max()],origin='lower',linestyles='-',alpha=0.6,vmin=150,vmax=500,cmap='magma',linewidth=0.7)
#ax.set_title('%s-band'%(labs[1]),size=15)
ax.set_ylim(80,450)


ax.set_xlim(0,10)
ax.set_yticks([100,200,300,400])
ax.set_xticks([0,2,4,6,8,10])
xticks = np.array([0,2,4,6,8,10])
ax2 = ax.twiny()
ax2.set_xlim(np.array([0,10])*9*3600/(450) * 1)
ax2.set_xlabel(r'Max. area searchable with WINTER [deg$^{2}$]',size=14,labelpad=10)
#ax2.set_xticks(xticks*8*3600/450 * 1)

ax.plot(area90s/(9*3600/(450)),distances,'x',color=color_list[-1],markersize=7)
#ax.plot(area50s/(9*3600/(450)),distances,'x',color='red',markersize=7)

ax.set_ylabel(r'Distance [Mpc]',size=16)
ax.set_xlabel(r'Time searched [days]',size=16)


ax.tick_params(labelsize=12)
ax2.tick_params(labelsize=12)
ax.text(4,400,'J-band',size=15)
ax.text(4,380,r'm$_{\rm{lim}} = 21$',size=12)
#ax.text(7.5,350,'< 21 mag',size=15)

cbar_ax = fig.add_axes([0.93,0.15,0.02,0.75])
cbar = fig.colorbar(im,cax=cbar_ax)
cbar.set_label(r'Detected fraction (f$_{\rm{det}}$)',size=15)
#cbar.set_ticks([150,200,250,300,350,400,450,500])
#cbar.set_ticklabels([150,200,250,300,350,400,450,500])
cbar.ax.tick_params(labelsize=15,pad=5)

ax.plot(0,0, label=r'f$_{\rm{det}} = 0.1$', color='black', ls='--')
ax.plot(0,0, label=r'f$_{\rm{det}} = 0.5$', color='black', ls='-')
ax.plot(0,0, label=r'f$_{\rm{det}} = 0.9$', color='black', ls=':')

ax.legend(fontsize=12,loc=1)

plt.savefig(r'tloc_dist_comparisons.pdf',bbox_inches='tight')


# In[5]:


def gen_det_fracs(mlim,tablelist,DMs,t_covs):
    for ind in range(len(tablelist)):
        t = tablelist[ind]
        t1 = t

        det_fracs = [[],[]]
        for DM in DMs:
            fdetarray = []
            for t_cov in t_covs:
                mapps = t1['mag'][:,0]+DM
                t_ind = np.argmin(np.abs(t1[0]['t']-t_cov))
                m_at_cov = mapps[:,t_ind]
                ndet = np.sum(m_at_cov<mlim)
                fdet = ndet/len(m_at_cov)
                    
                fdetarray.append(fdet)
            det_fracs[ind].append(fdetarray)    
        
        det_fracs = np.array(det_fracs)
    return det_fracs


# In[163]:


cmap = colors.LinearSegmentedColormap.from_list("",['ghostwhite','cornflowerblue','darkblue'])
labs = ['r','J']
Ds = np.linspace(20,500,20)
tablelist = [r,J]
DMs = 5*np.log10(Ds*1e5)
t_covs = np.arange(0,15,0.3)

fig = plt.figure(figsize=(7,6))
#gs = GridSpec(1,1,hspace=0.4,wspace=0.3)
        
#ax = fig.add_subplot(gs[0])
mlim = {}
mlim['450'] = 21
texp = 450

areas = 1*t_covs*9*3600/(texp) #sq. deg
det_fracs = gen_det_fracs(21,tablelist,DMs,t_covs)
#ax = fig.add_subplot(gs[0])
cs = plt.contour(areas,Ds,det_fracs[1],[0.1],colors='black',origin='lower',linestyles=['--','-'],linewidth=3,label='J = 21 mag')
#ax.set_ylim(20,450)

texp = 180
areas = 1*t_covs*9*3600/(texp)
det_fracs = gen_det_fracs(20.5,tablelist,DMs,t_covs)
cs = plt.contour(areas,Ds,det_fracs[1],[0.1],colors='blue',origin='lower',linestyles=['--','-'],linewidth=3,label='J = 20 mag')

texp = 40
areas = 1*t_covs*9*3600/(texp)
det_fracs = gen_det_fracs(19.5,tablelist,DMs,t_covs)
cs = plt.contour(areas,Ds,det_fracs[1],[0.1],colors=[color_list[-1]],origin='lower',linestyles=['--','-'],linewidth=3,label='J = 19 mag')
plt.xlabel('Area searchable with WINTER [deg$^{2}$]',size=14)
plt.ylabel('Distance [Mpc]',size=14)

plt.xscale('log')
plt.xlim(10)
plt.ylim(40,400)
plt.tick_params(labelsize=12)
#plt.text(1000,350,'J = 21 mag',c='black')
#plt.text(1000,320,'J = 20 mag',c='blue')
#plt.text(1000,290,'J = 19 mag',c='red')
#ax.set_ylim(20,450)
#ax.set_xlim(0,10)
#ax.set_yticks([100,200,300,400])
#ax.set_xticks([0,2,4,6,8,10])
#xticks = np.array([0,2,4,6,8,10])

#ax2 = ax.twiny()
#ax2.set_xlim(np.array([0,10])*9*3600/(450) * 1)
#ax2.set_xlabel(r'Area searched with WINTER [sq. deg.]',size=12)
#ax2.set_xticks(xticks*8*3600/450 * 1)

#ax.plot(sched['area90']/(8*3600/(450)),sched['distance'],'x',color='red',markersize=7)
#ax.set_ylabel(r'Distance [Mpc]',size=12)
#ax.set_xlabel(r'Time allotted for search',size=12)


#ax.tick_params(size=6)
#ax.text(8,370,'J-band',size=15)
#ax.text(7.5,350,'< 21 mag',size=15)

#cbar_ax = fig.add_axes([0.93,0.15,0.02,0.75])
#cbar = fig.colorbar(im,cax=cbar_ax)
#cbar.set_label('Detected fraction',size=15)
#cbar.set_ticks([150,200,250,300,350,400,450,500])
#cbar.set_ticklabels([150,200,250,300,350,400,450,500])
#cbar.ax.tick_params(labelsize=15,pad=5)

plt.plot(0,0, label=r'J = 21  , f$_{\rm{det}} = 0.1$', color='black', ls='--')
#plt.plot(0,0, label=r'J = 21, f$_{\rm{det}} = 0.5$', color='black', ls='-')

plt.plot(0,0, label=r'J = 20.5, f$_{\rm{det}} = 0.1$', color='blue', ls='--')
#plt.plot(0,0, label=r'J = 20.5, f$_{\rm{det}} = 0.5$', color='blue', ls='-')

plt.plot(0,0, label=r'J = 19.5, f$_{\rm{det}} = 0.1$', color=color_list[-1], ls='--')
#plt.plot(0,0, label=r'J = 19.5, f$_{\rm{det}} = 0.5$', color='red', ls='-')

plt.legend(fontsize=12)
plt.savefig(r'tloc_dist_comparisons_contours_Jband.pdf',bbox_inches='tight')


# In[153]:


cmap = colors.LinearSegmentedColormap.from_list("",['ghostwhite','cornflowerblue','darkblue'])
labs = ['r','J']
Ds = np.linspace(20,500,20)
tablelist = [r,J]
DMs = 5*np.log10(Ds*1e5)
t_covs = np.arange(0,15,0.3)

fig = plt.figure(figsize=(7,6))
#gs = GridSpec(1,1,hspace=0.4,wspace=0.3)
        
#ax = fig.add_subplot(gs[0])
mlim = {}
mlim['450'] = 21
texp = 450

areas_plot = np.array([100,200,500,1000,2000])
area_labs = ['A','B','C','D','E']
t_plot_21 = areas_plot/(1*9*3600/450)
t_plot_20 = areas_plot/(1*9*3600/120)
t_plot_19 = areas_plot/(1*9*3600/60)

areas = 1*t_covs*9*3600/(texp) #sq. deg
det_fracs = gen_det_fracs(21,tablelist,DMs,t_covs)
#ax = fig.add_subplot(gs[0])
cs = plt.contour(t_covs,Ds,det_fracs[1],[0.1],colors='black',origin='lower',linestyles=['--','-'],linewidth=3,label='J = 21 mag')
#ax.set_ylim(20,450)

X, Y = cs.collections[0].get_paths()[0].vertices.T
for idx,t_plot in enumerate(t_plot_21):
    if t_plot>10 or t_plot < 0:
        continue
    ind = np.argmin(np.abs(X-t_plot))
    plt.plot(X[ind],Y[ind],'.',color='black',markersize=10)
    plt.text(X[ind]-0.1,Y[ind]+10,area_labs[idx],color='black',size=10)
    
texp = 180
areas = 1*t_covs*9*3600/(texp)
det_fracs = gen_det_fracs(20,tablelist,DMs,t_covs)
cs = plt.contour(t_covs,Ds,det_fracs[1],[0.1],colors='blue',origin='lower',linestyles=['--','-'],linewidth=3,label='J = 20 mag')

X, Y = cs.collections[0].get_paths()[0].vertices.T
for idx,t_plot in enumerate(t_plot_20):
    if t_plot>10 or t_plot < 0:
        continue
    ind = np.argmin(np.abs(X-t_plot))
    plt.plot(X[ind],Y[ind],'.',color='blue',markersize=10)
    plt.text(X[ind]-0.1,Y[ind]+10,area_labs[idx],color='blue',size=10)
    
texp = 40
areas = 1*t_covs*9*3600/(texp)
det_fracs = gen_det_fracs(19,tablelist,DMs,t_covs)
cs = plt.contour(t_covs,Ds,det_fracs[1],[0.1],colors=[color_list[-1]],origin='lower',linestyles=['--','-'],linewidth=3,label='J = 19 mag')

X, Y = cs.collections[0].get_paths()[0].vertices.T
for idx,t_plot in enumerate(t_plot_19):
    if t_plot>10 or t_plot < 0:
        continue
    ind = np.argmin(np.abs(X-t_plot))
    plt.plot(X[ind],Y[ind],'.',color='red',markersize=10)
    plt.text(X[ind]-0.1,Y[ind]+10,area_labs[idx],color=color_list[-1],size=10)
    
plt.xlabel('Time searched',size=12)
plt.ylabel('Distance [Mpc]',size=12)



#plt.vlines(t_plot_21,ymin=350,ymax=400,color='black',linestyle='--')
#plt.vlines(t_plot_20,ymin=350,ymax=400,color='blue',linestyle='--')
#plt.vlines(t_plot_19,ymin=350,ymax=400,color='red',linestyle='--')

#plt.xscale('log')
plt.xlim(0,10)
plt.ylim(40,400)
#plt.text(1000,350,'J = 21 mag',c='black')
#plt.text(1000,320,'J = 20 mag',c='blue')
#plt.text(1000,290,'J = 19 mag',c='red')

#ax.set_ylim(20,450)
#ax.set_xlim(0,10)
#ax.set_yticks([100,200,300,400])
#ax.set_xticks([0,2,4,6,8,10])
#xticks = np.array([0,2,4,6,8,10])

#ax2 = ax.twiny()
#ax2.set_xlim(np.array([0,10])*9*3600/(450) * 1)
#ax2.set_xlabel(r'Area searched with WINTER [sq. deg.]',size=12)
#ax2.set_xticks(xticks*8*3600/450 * 1)

#ax.plot(sched['area90']/(8*3600/(450)),sched['distance'],'x',color='red',markersize=7)
#ax.set_ylabel(r'Distance [Mpc]',size=12)
#ax.set_xlabel(r'Time allotted for search',size=12)


#ax.tick_params(size=6)
#ax.text(8,370,'J-band',size=15)
#ax.text(7.5,350,'< 21 mag',size=15)

#cbar_ax = fig.add_axes([0.93,0.15,0.02,0.75])
#cbar = fig.colorbar(im,cax=cbar_ax)
#cbar.set_label('Detected fraction',size=15)
#cbar.set_ticks([150,200,250,300,350,400,450,500])
#cbar.set_ticklabels([150,200,250,300,350,400,450,500])
#cbar.ax.tick_params(labelsize=15,pad=5)

plt.plot(0,0, label=r'J = 21, f$_{\rm{det}} = 0.1$', color='black', ls='--')
#plt.plot(0,0, label=r'J = 21, f$_{\rm{det}} = 0.5$', color='black', ls='-')

plt.plot(0,0, label=r'J = 20.5, f$_{\rm{det}} = 0.1$', color='blue', ls='--')
#plt.plot(0,0, label=r'J = 20.5, f$_{\rm{det}} = 0.5$', color='blue', ls='-')

plt.plot(0,0, label=r'J = 19.5, f$_{\rm{det}} = 0.1$', color=color_list[-1], ls='--')
#plt.plot(0,0, label=r'J = 19.5, f$_{\rm{det}} = 0.5$', color='red', ls='-')

plt.text(7.2,250,r'A :  100 sq. deg')
plt.text(7.2,230,r'B :  200 sq. deg')
plt.text(7.2,210,r'C :  500 sq. deg')
plt.text(7.2,190,r'D : 1000 sq. deg')
plt.text(7.2,170,r'E : 2000 sq. deg')
plt.legend(fontsize=12)
plt.savefig(r'../output/tloc_dist_comparisons_contours_Jband_loc_areas.pdf',bbox_inches='tight')

# In[4]:



# In[6]:
