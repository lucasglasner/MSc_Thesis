#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 12:02:56 2022

@author: lucas
# =============================================================================
# 
# =============================================================================
"""

import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import geopandas as gpd
import numpy as np
import re
from scipy.interpolate import interp1d
import matplotlib as mpl
import scipy.stats as st
import cartopy.crs as ccrs
import re
from functions import local_minimum_filter
from sklearn.decomposition import PCA
plt.rc('font',size=18)

interval = slice("2001-09-08","2015-12-31")
#%%
# =============================================================================
# CATCHMENT ATTRIBUTES
# =============================================================================
basins = pd.read_csv('datos/main_basins.txt',header=None)
basins = list(np.array(basins).squeeze())
basin_attrs = pd.read_csv('datos/basins_attributes.csv')
basin_attrs.index = basin_attrs.gauge_name
index = [" ".join(re.findall('[A-Z][a-z]*', b)) for b in basins]
index[-1] = 'Rio Uble En San Fabian N 2'
basin_attrs = basin_attrs.loc[index]

basins.pop(-2)

polygons = [gpd.read_file('datos/vector/basins/'+b+'.shp') for b
            in basins]
polygons = pd.concat(polygons)


hypso = []
areas = []
for b in basins:
    path = 'datos/topography/basins/hypso/'+b+'_hypso.csv'
    h = pd.read_csv(path)
    h.index = h.height
    h.drop('height',axis=1,inplace=True)
    areas.append(h.iloc[:,1].max())
    hypso.append(h.fArea)
    
areas = pd.Series(areas,index=basins)
hypso = pd.concat(hypso,axis=1)
hypso.columns=basins

int_func = [interp1d(hypso.index,hypso[b]) for b in basins]
int_func = pd.Series(int_func,index=basins)

int_func_inv = [interp1d(hypso[b],hypso.index) for b in basins]
int_func_inv = pd.Series(int_func_inv,index=basins)
del index
#%%
# =============================================================================
# SNOW
# =============================================================================
SCA = []
for b in basins:
    path = 'datos/ianigla/'+b+'_SCA_s_comp.filtro_MA.3días.csv'
    sca = pd.read_csv(path)
    sca.index = pd.to_datetime(sca.fecha)
    SCA.append(sca.iloc[:,1])
SCA = pd.concat(SCA,axis=1).dropna()/100
SCA.columns = basins

SCA_trend = SCA.shift(-1)-SCA.shift(1)

SL = [int_func_inv[b](1-SCA[b]) for b in basins]
SL = pd.DataFrame(SL,index=basins,columns=SCA.index).T
#%%
# =============================================================================
# RAIN
# =============================================================================
pr_cr2met = []
for b in basins:
    path='datos/estaciones/cr2met/pr_'+b+'_2000-2020.csv'
    pr = pd.read_csv(path)
    pr.index = pd.to_datetime(pr.time)
    pr_cr2met.append(pr.iloc[:,1])

pr_cr2met = pd.concat(pr_cr2met,axis=1)
pr_cr2met.columns = basins

pr_sanjose = pd.read_csv('datos/estaciones/pr_SanJoseMaipo.csv',dtype=str)
pr_sanjose.index = pr_sanjose.iloc[:,0]+"-"+pr_sanjose.iloc[:,1]+"-"+pr_sanjose.iloc[:,2]
# pr_sanjose.index = pd.to
pr_sanjose.index = pd.to_datetime(pr_sanjose.index)
pr_sanjose = pd.to_numeric(pr_sanjose.iloc[:,3])


pr_stations_h = []
for b in basins:
    path = 'datos/estaciones/vismet/pr_'+b+'.csv'
    pr = pd.read_csv(path)
    pr.index = pd.to_datetime(pr.Fecha)
    pr_stations_h.append(pr['Valor'])
pr_stations_h = pd.concat(pr_stations_h,axis=1)
pr_stations_h.columns = basins

pr_stations = pr_stations_h.resample('d').max()

pr_30cum = pd.DataFrame(np.empty(pr_cr2met.shape),
                        index=pr_cr2met.index,
                        columns=pr_cr2met.columns)


x = pr_cr2met.cumsum()
for i in range(len(pr_30cum)):
    if i>30:
        pr_30cum.iloc[i,:] = x.iloc[i,:]-x.iloc[i-30,:]
    else:
        pr_30cum.iloc[i,:] = x.iloc[i,:]

#%%
# =============================================================================
# FREEZING LEVEL
# =============================================================================
FL = []
for b in basins:
    path='datos/era5/horarias/H0_'+b+'_2000-2020.csv'
    fl = pd.read_csv(path)
    fl.index = pd.to_datetime(fl.time)
    FL.append(fl.iloc[:,1])
FL = pd.concat(FL, axis=1)-300
FL.columns = basins
FL = FL.resample('d').mean()

pluv_area = [int_func[b](FL[b]) for b in basins]
pluv_area = pd.DataFrame(pluv_area,columns=FL.index,index=basins).T
pluv_area = pluv_area.where(pr_cr2met>3)

#%%
# =============================================================================
# RUNOFF
# =============================================================================

Q = []
for b in basins:
    path='datos/estaciones/dga/qinst_'+b+'.csv'
    try:
        q = pd.read_csv(path)
        q.index = pd.to_datetime(q.iloc[:,0])
        Q.append(q.iloc[:,1])
    except:
        Q.append(pd.Series([np.nan]))
        
Q = pd.concat(Q,axis=1)
Q.columns = basins

Qmaxd  = Q.resample('d').max()
Qmeand = Q.resample('d').mean()



baseflow = []
quickflow = []
for b in basins:
    b,q = local_minimum_filter(Q[b].dropna(),25)
    baseflow.append(b)
    quickflow.append(q)
baseflow = pd.concat(baseflow,axis=1)
quickflow= pd.concat(quickflow,axis=1)

baseflow.columns = basins
quickflow.columns = basins


mean_quickflow = quickflow.resample('d').mean()

#%%
# =============================================================================
# MERGE DAILY DATA    
# =============================================================================
    
data_daily = []
for i,b in enumerate(basins):
    d = pd.concat([pr_cr2met[b],SCA[b],SCA_trend[b],SL[b],FL[b],
                   pluv_area[b],pr_30cum[b],pr_stations[b],
                   Qmaxd[b],Qmeand[b],mean_quickflow[b]],
                  axis=1)
    ros_area = (SCA[b]-(1-pluv_area[b]))
    ros_area = ros_area.where(ros_area>0).fillna(0)
    d = pd.concat([d,ros_area],axis=1)
    d.columns = ['pr_cr2met','SCA','SCA_trend','SL','FL',
                 'pluv_area','pr_30cum','pr_max',
                 'Qmaxd','Qmeand','quickflow','ros_area']
    data_daily.append(d[interval])

data_daily = pd.concat(data_daily,keys=basins)
data_daily['ros_area'] = data_daily['ros_area'].where(data_daily['pr_cr2met']>3)
data_daily['ros_area'] = data_daily['ros_area'].where(data_daily['SL']<data_daily['FL'])

data_daily['delta_ros'] = data_daily['FL']-data_daily['SL']

data_daily['date'] = data_daily.index.get_level_values(1)
data_daily['new_events'] = data_daily.groupby(
    'pr_cr2met')['date'].apply(lambda s: s.diff().dt.days < 2)
data_daily['new_events'] = data_daily['new_events'].rolling(
    2, center=True).min().astype(bool)
data_daily['events'] = data_daily['new_events'].cumsum()
data_daily = data_daily.groupby([data_daily.index.get_level_values(0),
                                 data_daily.index.get_level_values(1),
                                 data_daily.events]).mean()
# data_daily = data_daily.where(data_daily.quickflow>0)

#%%
data_events = data_daily.quickflow.unstack(level=1).max(axis=1)
data_events = pd.concat([data_events], axis=1)
data_events.columns = ['quickflow']
data_events = data_events[data_events['quickflow']>0]

data_events = [data_events.loc[b] for b in basins]


#%%
for i,b in enumerate(basins):
    #precipitation acumulated in the event
    data_events[i]['pr_cum'] = data_daily.loc[b].pr_cr2met.groupby(
        'events').apply(lambda x: x.sum())
    #max precipitation intensity if there is in a nearby station
    data_events[i]['pr_max'] = data_daily.loc[b].pr_max.groupby(
        'events').apply(lambda x: x.max())
    #event start date
    data_events[i]['start'] = data_daily.loc[b].SCA.groupby(
        'events').apply(lambda x: x.index[0][0])
    #event last date
    data_events[i]['end'] = data_daily.loc[b].SCA.groupby(
        'events').apply(lambda x: x.index[-1][0])
    #event duration in days
    data_events[i]['duration'] = data_events[i]['end']-data_events[i]['start']
    data_events[i]['duration'] = data_events[i].duration.apply(
        lambda x: x.total_seconds()/3600)
    #initial snow cover
    data_events[i]['start_sca'] = data_daily.loc[b].SCA.loc[data_events[i].start].values
    #final snow cover
    data_events[i]['end_sca'] = data_daily.loc[b].SCA.loc[data_events[i].end].values
    #maximum pluvial area on event
    data_events[i]['max_pluv_area'] = data_daily.loc[b].pluv_area.groupby(
        'events').apply(lambda x: x.max())
    #minimum pluvial area on event
    data_events[i]['min_pluv_area'] = data_daily.loc[b].pluv_area.groupby(
        'events').apply(lambda x: x.min())
    #maximum ros area on event
    data_events[i]['max_ros_area'] = data_daily.loc[b].ros_area.groupby(
        'events').apply(lambda x: x.max())
    #mean ros area on event
    data_events[i]['mean_ros_area'] = data_daily.loc[b].ros_area.groupby(
        'events').apply(lambda x: x.mean())
    #maximum runoff in the event
    data_events[i]['Qmaxd'] = data_daily.loc[b]['Qmaxd'].groupby(
        'events').apply(lambda x: x.max())
    #return period of the maxmimum runoff
    data_events[i]['exc_prob'] = data_events[i]['Qmaxd'].apply(
        lambda x: (st.percentileofscore(Qmaxd[b].dropna(),x)))
    data_events[i]['return_period'] = 1/(1-data_events[i]['exc_prob']/100)
    #direct runoff volume
    data_events[i]['direct_runoff_volume'] = data_daily.loc[b].groupby(
        'events').apply(lambda x: np.sum(x.quickflow*3600*24/1e6))
    #aviable rain volume
    data_events[i]['precipitated_volume'] = data_daily.loc[b].groupby(
        'events').apply(lambda x: areas.loc[b]*np.sum(x['pr_cr2met']*x['pluv_area'])/1e3)
    #soil moisture proxy
    data_events[i]['mean_pr30cum'] = data_daily.loc[b].pr_30cum.groupby(
        'events').apply(lambda x: x.mean())
    #freezing level
    data_events[i]['FL'] = data_daily.loc[b].FL.groupby(
        'events').apply(lambda x: x.mean())

data_events = pd.concat(data_events,keys=basins)
data_events['basin'] = data_events.index.get_level_values(0)
# data_events = data_events[data_events['return_period']>5]
data_events['delta_sca'] = data_events['end_sca']-data_events['start_sca']

data_events = data_events.iloc[np.where(data_events.direct_runoff_volume>0)[0]]
# data_events = data_events[data_events['pr_cum']>10]

# %%
# =============================================================================
# PCA
# =============================================================================
from scipy.cluster.vq import whiten
from sklearn.preprocessing import StandardScaler
from pca import pca



# PCA = pca(len(names))
# p = PCA.fit_transform(whiten(matrix[var]), col_labels=names)
# pc = p['PC']
# load = p['loadings']
# r2 = pd.Series(p['variance_ratio'])


# markers = ['o','^','s',"v",'*','P','X','D','1']

titles = [re.findall('[A-Z][^A-Z]*', b)[1] for b in basins]
titles[-1] = 'Ñuble'

# fig,ax = plt.subplots(1,1)
# # ax = [ax]
# ax.set_xlabel('PC1 ('+'{:.0%}'.format(r2[0])+')', fontsize=14)
# ax.set_ylabel('PC2 ('+'{:.0%}'.format(r2[1])+')', fontsize=14)

# mask = matrix[matrix.basin == b].dropna().max_ros_area>0.2
# comp = pc.iloc[:,:2]
# loadings = load.iloc[:2,:]
# sx = 1/(comp.iloc[:,0].max()-comp.iloc[:,0].min())
# sy = 1/(comp.iloc[:,1].max()-comp.iloc[:,0].min())

# for m,b in zip(markers,matrix.basin.unique()):
#     mask = matrix.basin==b
#     mask = mask.values
#     ax.scatter(comp.iloc[mask,0]*sx,comp.iloc[mask,1]*sy,
#                 s=12,marker=m,
#                 alpha=0.4, color='k')
# # ax[i].scatter(comp[mask.values].iloc[:,0]*sx,comp[mask.values].iloc[:,1]*sy,
# #               s=8, alpha=0.5, color='tab:blue')

# for j in range(len(loadings.columns)):
#     # ax[i].annotate("",xy=(0,0),
#     #                 xytext=(loadings.iloc[0,j]*.9,
#     #                         loadings.iloc[1,j]*.9),
#     #                 arrowprops=dict(headwidth=0, width=0.05,
#     #                                 color=colors[j]))
#     ax.annotate(names[j],xy=(0,0),xytext=(loadings.iloc[0,j],
#                                     loadings.iloc[1,j]),
#                 arrowprops=dict(arrowstyle="-",
#                                 color=colors[j]),
#                 fontsize=10)
# # yabs_max = abs(max(ax[i].get_ylim(), key=abs))
# # xabs_max = abs(max(ax[i].get_xlim(), key=abs))

# xmin = min([loadings.T.min().values[0]*1.1, comp.iloc[:,0].min()*sx*1.1])
# xmax = max([loadings.T.max().values[0]*1.1, comp.iloc[:,0].max()*sx*1.1])


# ymin = min([loadings.T.min().values[1]*1.1, comp.iloc[:,1].min()*sy*1.1])
# ymax = max([loadings.T.max().values[1]*1.1, comp.iloc[:,1].max()*sy*1.1])

# ax.set_ylim(ymin=ymin,
#                ymax=ymax)
# ax.set_xlim(xmin=xmin,
#                xmax=xmax)
# # ax.set_title(titles[i],loc='left')
# # ax[i].set_xlim(-0.8,0.8)
# # ax[i].set_ylim(-0.8,0.8)

var = ['direct_runoff_volume','precipitated_volume','mean_pr30cum',
        'delta_sca','duration','max_ros_area']

names = [r"$\int Q' dt$",r'$\int PR \cdot A_{p}dt$',r'$PR_{30d}$',
         r'$\Delta SCA$','Duration',r'$A_{ROS}$']

matrix = data_events.droplevel(0)
matrix = matrix[matrix.return_period>100]


PCA = []
pc = [] 
load = []
r2 = []

for b in matrix.dropna().basin.unique():
    v = matrix[matrix['basin']==b][var].dropna()
    # scaler = StandardScaler().fit(v.values)
    # scaler = scaler.transform(v.values)
    scaler = whiten(v.values)
    pca_x = pca(n_components = len(names))
    p = pca_x.fit_transform(scaler, col_labels=names)
    p['PC'].index = v.index
    PCA.append(pca_x)
    load.append(p['loadings'])
    pc.append(p['PC'])
    r2.append(pd.Series(p['variance_ratio']))
#     pc.append(pd.DataFrame(p))
#     load.append(pd.DataFrame(pca_x.components_))
#     r2.append(pd.Series(pca_x.explained_variance_ratio_))
    
   
pc = pd.concat(pc,keys=matrix.dropna().basin.unique())    
load = pd.concat(load, keys=matrix.dropna().basin.unique())
r2 = pd.concat(r2, keys=matrix.dropna().basin.unique())    
# #%%

# fig, ax = plt.subplots(2,4, sharex=True,sharey=True, figsize=(12,4))
# ax = ax.ravel()

# for i,b in enumerate(matrix.dropna().basin.unique()):
#     ax[i].bar(np.arange(len(names)),r2.loc[b].cumsum(),ec='k',width=0.5)
#     ax[i].set_xticks(np.arange(len(names)))
#     ax[i].set_yticks(np.arange(0,1+0.25,0.25))
#     ax[i].grid(axis='y', ls=":")
    
#%%
colors = mpl.cm.nipy_spectral(np.linspace(0,.9,len(names)))

# colors = mpl.colors.BASE_COLORS.values()
colors = list(colors)

fig,ax = plt.subplots(2,4, figsize=(14,6))
fig.tight_layout(pad=2)

plt.rc('font',size=18)
ax = ax.ravel()
for i,b in enumerate(matrix.dropna().basin.unique()):
    ax[i].set_xlabel('PC1 ('+'{:.0%}'.format(r2.loc[b][0])+')', fontsize=14)
    ax[i].set_ylabel('PC2 ('+'{:.0%}'.format(r2.loc[b][1])+')', fontsize=14)
    
    ax[i].text(0,0,'{:.0%}'.format(r2.loc[b][0]+r2.loc[b][1]),
               transform=ax[i].transAxes, fontsize=12)
    
    comp = pc.loc[b].iloc[:,[0,1]]
    loadings = load.loc[b].iloc[[0,1],:]
    sx = 1/(comp.iloc[:,0].max()-comp.iloc[:,0].min())
    sy = 1/(comp.iloc[:,1].max()-comp.iloc[:,0].min())
    
    
    mask = (matrix[matrix.basin==b].max_ros_area>0.1)
    mask = (matrix[matrix.basin==b].delta_sca<0) & mask
    mask = mask.reindex(comp.index)
    for j in range(len(loadings.columns)):
        # ax[i].annotate("",xy=(0,0),
        #                 xytext=(loadings.iloc[0,j]*.9,
        #                         loadings.iloc[1,j]*.9),
        #                 arrowprops=dict(headwidth=0, width=0.05,
        #                                 color=colors[j]))
        ax[i].annotate("",xy=(0,0),xytext=(loadings.iloc[0,j],
                                           loadings.iloc[1,j]),
                       arrowprops=dict(arrowstyle="<-",
                                       color=colors[j],
                                       linewidth=2,
                                       mutation_scale=9))
    
    ax[i].scatter(comp[~mask].iloc[:,0]*sx,comp[~mask].iloc[:,1]*sy,
                  s=5,
                  alpha=0.5, color='grey')
    ax[i].scatter(comp[mask].iloc[:,0]*sx,comp[mask].iloc[:,1]*sy,
                   s=5, c = 'darkred')
                  # s=20, c = matrix.loc[mask[mask].index].max_ros_area)
    # ax[i].scatter(comp[mask.values].iloc[:,0]*sx,comp[mask.values].iloc[:,1]*sy,
    #               s=8, alpha=0.5, color='tab:blue')
    

    yabs_max = abs(max(ax[i].get_ylim(), key=abs))
    xabs_max = abs(max(ax[i].get_xlim(), key=abs))
    
    xmin = min([loadings.T.min().values[0]*1.5, comp.iloc[:,0].min()*sx*1.5])
    xmax = max([loadings.T.max().values[0]*1.1, comp.iloc[:,0].max()*sx*1.1])
    
    
    ymin = min([loadings.T.min().values[1]*1.1, comp.iloc[:,1].min()*sy*1.1])
    ymax = max([loadings.T.max().values[1]*1.1, comp.iloc[:,1].max()*sy*1.1])
    
    ax[i].set_ylim(ymin=ymin,
                   ymax=ymax)
    ax[i].set_xlim(xmin=xmin,
                   xmax=xmax)
    ax[i].set_title(titles[i],loc='left')
    
    ax[i].axvline(0,ls=":", color='grey')
    ax[i].axhline(0,ls=":",color='grey')
    # ax[i].set_xticks([0])
    # ax[i].set_yticks([0])
    # ax[i].grid(True, ls=":")
    # ax[i].set_xlim(-0.8,0.8)
    # ax[i].set_ylim(-0.8,0.8)

for j in range(len(names)):
    ax[0].plot([],[], color=colors[j], label=names[j])
ax[0].legend(frameon=False, loc=(0,-2.5), ncol=3)
    

# #plt.savefig('plots/BIPLOTS_PC1PC2.pdf',dpi=150,bbox_inches='tight')

#%%
fig,axis = plt.subplots(2,4,figsize=(12,6),sharex=True,sharey=True)
fig.tight_layout(pad=1)
plt.rc('font', size=18)

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, LogLocator, NullFormatter)
ax = axis.ravel()
v = data_events
c = 'max_ros_area'



for i,b in enumerate(basins):
    mask = v.loc[b]['max_ros_area']>0.
    mask = mask & (v.loc[b].delta_sca<0)
    ax[i].set_title(titles[i], loc='left')
    ax[i].set_yscale('log')
    ax[i].set_xscale('log')
    ax[i].set_xlim(0.25,2000)
    ax[i].set_ylim(0.1,1000)
    ax[i].set_xticks([1,10,100,1000])
    # ax[i].xaxis.set_minor_locator(LogLocator(base=10,
                                             # subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 )))
    # ax[i].xaxis.set_minor_formatter(NullFormatter())


    ax[i].scatter(v.loc[b][~mask]['precipitated_volume'],
                  v.loc[b][~mask]['direct_runoff_volume'],
                  c='grey', alpha=0.5,
                  s=v.loc[b][~mask]['mean_pr30cum']*1)
    
    m = ax[i].scatter(v.loc[b][mask]['precipitated_volume'],
                      v.loc[b][mask]['direct_runoff_volume'],
                      c=v.loc[b][mask][c],
                      ec='k', zorder=1,
                      s=v.loc[b][mask]['mean_pr30cum']*1,
                      alpha=0.8,
                      cmap='Blues')
    ax[i].plot([1e-2,1e4],[1e-2,1e4], color='r', ls="--")

# ax[-1].axis('off')
box = axis[0,-1].get_position()
box1 = axis[-1,-1].get_position()
cax = fig.add_axes([box.xmax*1.05,box1.ymin,0.02,box.ymax-box1.ymin])

fig.colorbar(m,cax=cax, label='Maximum ROS Area (-)')

fig.text(0.5,-0.05,r"Precipitated Volume $(hm^3)$"+"\n $\int PR \cdot A_p dt$", ha='center',va='center')
fig.text(0,0.5,r"Direct Runoff Volume $(hm^3)$"+"\n $\int Q' dt$", ha='center',va='center',
         rotation='vertical')


# produce a legend with a cross section of sizes from the scatter
handles, labels = m.legend_elements(prop="sizes", alpha=0.5)
legend2 = ax[3].legend(handles, labels, loc=(1.8,-1),
                       title="30 days\nacumulated\nprecipitation\n(mm)",
                       frameon=False)

# #plt.savefig('plots/QvsPRvolume_plots.pdf',dpi=150,bbox_inches='tight')

#%%

ros_areas = (data_daily.ros_area>0.1).droplevel(2).unstack().T
# ros_areas = (ros_areas) & (data_daily.SCA_trend<0).droplevel(2).unstack().T

ac = ros_areas.groupby([ros_areas.index.month,ros_areas.index.year]).sum()
ac = ac.mean(level=0)

fig,ax = plt.subplots(2,4,sharex=True,sharey=True,figsize=(14,5))

ax = ax.ravel()

colors=mpl.cm.nipy_spectral(np.linspace(0.1,0.9,len(ac.columns)))
for i,b in enumerate(ac.columns):
    ax[i].bar(np.arange(12),ac[b],edgecolor='k',
              color=colors[i], label=titles[i])
    # ax[i].set_title(titles[i],loc='left')
    ax[i].legend(frameon=False, fontsize=12,loc='upper left')
    ax[i].set_xticks(np.arange(12))
    ax[i].set_xticklabels(['JAN','FEB','MAR','APR','MAY',
                           'JUN','JUL','AUG','SEP','OCT',
                           'NOV','DIC'],fontsize=13)
    ax[i].tick_params(axis='x',rotation=45)
    
fig.text(0.08,0.5,'Mean anual cycle of ROS days', rotation=90, ha='center',va='center')
    
# #plt.savefig('plots/ros_anualcycle.pdf',dpi=150,bbox_inches='tight')


#%%

ycum = ros_areas.groupby(ros_areas.index.year).sum()
fig,ax = plt.subplots(2,4,sharex=True,sharey=True,figsize=(14,5))
ax = ax.ravel()
for i,b in enumerate(ac.columns):
    m = st.linregress(ycum[b].index,ycum[b])
    ax[i].plot(ycum[b].index,ycum[b].index*m.slope+m.intercept,
               color='k', ls="--")
    ax[i].plot(ycum[b].index,ycum[b], color=colors[i], label=titles[i])
    ax[i].legend(frameon=False, fontsize=12)
    ax[i].set_yticks(np.arange(0,50+10,10))
    ax[i].set_xticks(np.arange(2001,2021,3))
    ax[i].tick_params(axis='x',rotation=45)
    ax[i].text(0,1,'pvalue: '+'{:.0%}'.format(m.pvalue),transform=ax[i].transAxes,
               fontsize=14)
fig.text(0.08,0.5,'#ROS days per year', rotation=90, ha='center',va='center')
# #plt.savefig('plots/ROS_days_per_year.pdf',dpi=150,bbox_inches='tight')

#%%
dd = data_daily.loc['RioMaipoEnElManzano'].droplevel(1)
ros_days = dd.ros_area > 0.1
ros_days = ros_days.groupby([ros_days.index.year,ros_days.index.month]).sum()
ros_days = ros_days.unstack().T

de = data_events.loc['RioMaipoEnElManzano']

plt.rc('font',size=18)
fig,ax = plt.subplots(1,4,figsize=(18,3))
ax[0].bar(np.arange(12),ros_days.mean(axis=1),yerr=1.96*ros_days.std(axis=1)/np.sqrt(14),
          capsize=3, edgecolor='k')
ax[0].set_ylim(0,3.5)
ax[0].text(0,1.06,"Mean #Events: "+"{:.01f}".format(de[de.max_ros_area>0.1].start.count()/14),
           transform=ax[0].transAxes)
ax[0].set_xticks(np.arange(12))
ax[0].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
ax[0].set_ylabel('Mean #ROS days')

ax[1].boxplot(de[(de.return_period>100)].Qmaxd, positions=[0],patch_artist=True,
              boxprops=dict(facecolor='tab:blue'),medianprops=dict(color='k'))
ax[1].boxplot(de[(de.return_period>100)&((de.delta_sca<0))].Qmaxd,
              positions=[0.5],patch_artist=True,
              boxprops=dict(facecolor='tab:purple'),medianprops=dict(color='k'))
ax[1].set_xticklabels(['100yr\nfloods','100yr\nROS\nfloods'])

# ax[2].hist(dd.SL[(dd.pr_cr2met>3)&(dd.ros_area>0)],alpha=0.5,density=True)
# ax[2].hist(dd.SL[dd.pr_cr2met>3],alpha=0.5,density=True)

ax[2].boxplot(dd.pluv_area[dd.pr_cr2met>3].dropna(),positions=[0],
              patch_artist=True,medianprops=dict(color='k'),boxprops=dict(facecolor='tab:red'),sym="")
ax[2].boxplot(dd.pluv_area[(dd.pr_cr2met>3) & (dd.ros_area>0.1) &(dd.SCA_trend<0)].dropna(),positions=[0.25],
              patch_artist=True,medianprops=dict(color='k'),boxprops=dict(facecolor='tab:purple'),sym="")


ax[2].boxplot(dd.SCA[dd.pr_cr2met>3].dropna(),positions=[1],
              patch_artist=True,medianprops=dict(color='k'),boxprops=dict(facecolor='cornflowerblue'),sym="")
ax[2].boxplot(dd.SCA[(dd.pr_cr2met>3) & (dd.ros_area>0) &(dd.SCA_trend<0)].dropna(),positions=[1.25],
              patch_artist=True,medianprops=dict(color='k'),boxprops=dict(facecolor='tab:purple'),sym="")

ax[2].set_xticks([0.1,1.1])
ax[2].set_xticklabels(['Pluvial\nArea','Snow\nCover\nArea'])

ax[2].plot([],[],color='tab:purple',label='ROS days')
ax[2].legend(frameon=False,loc=(0,1))

# ax[3].set_yscale('log')
pr_sanjose = pr_sanjose.reindex(dd.index)
ax[3].boxplot(pr_sanjose[pr_sanjose>3].dropna(),positions=[0], sym="",widths=0.3,
              patch_artist=True,medianprops=dict(color='k'),boxprops=dict(facecolor='tab:green'))
ax[3].boxplot(pr_sanjose[(pr_sanjose>3) & (dd.ros_area>0.1) & (dd.SCA_trend<0)].dropna(),positions=[2], sym="",widths=0.3,
              patch_artist=True,medianprops=dict(color='k'),boxprops=dict(facecolor='tab:purple'))
ax[3].set_xticklabels(['Daily\nprecipitation','Daily\nPrecipitation\non ROS\ndays'])
ax[3].set_yticks(np.arange(0,60,10))

ax[0].text(0.9,1.05,'(a)',transform=ax[0].transAxes)
ax[1].text(0.9,1.05,'(b)',transform=ax[1].transAxes)
ax[2].text(0.9,1.05,'(c)',transform=ax[2].transAxes)
ax[3].text(0.9,1.05,'(d)',transform=ax[3].transAxes)

# ax[3].set_ylim(1e-1,1e3)
# ax[3].boxplot(dd.pr_max[dd.pr_cr2met>3].dropna(),positions=[1])
# ax[3].boxplot(dd.pr_max[dd.pr_cr2met>3].dropna(),positions=[1.25])

#plt.savefig('plots/maipomanzano/ROS_FINAL.pdf',dpi=150,bbox_inches='tight')
#%%

polygons.index = basins
polygons['mean_ros_days'] = ycum.mean()
polygons['mean_rain_days'] = [(data_daily.pr_cr2met>3).loc[b].groupby(data_daily.loc[b].index.get_level_values(0).year).sum().mean() for b in basins]
polygons['ros_rain_ratio'] = polygons['mean_ros_days']/polygons['mean_rain_days']
polygons['max_runoff_ros_day'] = data_events[(data_events.delta_sca<0) & (data_events.max_ros_area>0.1)].Qmaxd.unstack().T.max()
polygons['mean_runoff_ros_day'] = data_events[(data_events.delta_sca<0) & (data_events.max_ros_area>0.1)].Qmaxd.unstack().T.median()
polygons['extreme_runoff'] = [np.percentile(data_events.Qmaxd[b],90) for b in basins]
polygons['qros_over_qextreme'] = polygons['mean_runoff_ros_day']/polygons['extreme_runoff']

polygons['max_deltasca_ros_day'] = -100*data_events[(data_events.delta_sca<0) & (data_events.max_ros_area>0.1)].delta_sca.unstack().T.min()
polygons['mean_deltasca_ros_day'] = -100*data_events[(data_events.delta_sca<0) & (data_events.max_ros_area>0.1)].delta_sca.unstack().T.mean()

polygons['p90'] = [np.percentile(data_events[data_events.pr_cum>1].pr_cum[b],90) for b in basins]
polygons['pmean_ROS'] = [data_events[data_events.max_ros_area>0.1].pr_cum[b].median() for b in basins]

polygons['prros_over_prextreme'] = polygons['pmean_ROS']/polygons['p90']

#%%
import cartopy.feature as cf
from functions import add_labels
import cmocean as cm

fig,ax = plt.subplots(1,7,sharex=True,sharey=True,figsize=(14,10),
                      subplot_kw={'projection':ccrs.PlateCarree()})

titles=['Mean ROS\ndays per\nyear','ROS/Rain\ndays\nratio',
        'Median runoff\non ROS days\nover extreme\nrunoff (Q90)','Maximum\nrunoff on\nROS days',
        'Mean SCA\nloss on\nROS days',
        'Maximum SCA\nloss on\nROS days',
        'Median\nprecipitation\non ROS days\nover extreme\nprecipitation\n(P90)'
        ]
cax = []
for axis,c in zip(ax,range(len(titles))):
    axis.coastlines()
    axis.set_extent([-72.1, -69.5, -32.4, -37.3])
    axis.add_feature(cf.OCEAN,rasterized=True)
    axis.add_feature(cf.LAND,rasterized=True)
    axis.add_feature(cf.BORDERS)
    
    box = axis.get_position()
    cax.append(fig.add_axes([box.xmin,box.ymin-0.12,box.xmax-box.xmin,
                             0.01]))
    # for b,color in zip(basins,colors):
    #     gpd.GeoSeries(polygons.boundary.loc[b]).plot(ax=axis,
    #                                                   color=color,
    #                                                   linewidth=0.75)
    polygons.boundary.plot(ax=axis,color='k',linewidth=0.5)
    axis.set_xticks([-72,-71,-70])
    axis.set_xticklabels(["72°W","71°W","70°W"])
    axis.tick_params(axis="x",labelsize=14,rotation=45)
    axis.set_title(titles[c],loc='left',fontsize=15)

polygons.plot(ax=ax[0],column='mean_ros_days',
              legend=True,
              cax=cax[0],
              legend_kwds={'orientation': "horizontal",
                           'ticks':[5,15,25],
                           'label':'(days)'})
polygons.plot(ax=ax[1],column='ros_rain_ratio',cmap='BuPu',
              legend=True,
              cax=cax[1],
              legend_kwds={'orientation': "horizontal",
                           'ticks':[0.2,0.3,0.4],
                           'label':'(-)'})
polygons.plot(ax=ax[2],column='qros_over_qextreme',cmap=cm.cm.matter,
              legend=True,
              cax=cax[2],vmax=1,
              legend_kwds={'orientation': "horizontal",
                           'ticks':[0.5,1],
                           'label':'(-)'})
polygons.plot(ax=ax[3],column='max_runoff_ros_day', cmap=cm.cm.turbid,
              legend=True,vmax=1.5e3,
              cax=cax[3],
              legend_kwds={'orientation': "horizontal",
                           'ticks':[500,1200],
                           'label':r'$(m3/s)$'})
polygons.plot(ax=ax[4],column='mean_deltasca_ros_day', cmap='Blues_r',
              legend=True,
              cax=cax[4],
              legend_kwds={'orientation': "horizontal",
                           'label':'(%)'})
polygons.plot(ax=ax[5],column='max_deltasca_ros_day', cmap=cm.cm.ice,
              legend=True,
              cax=cax[5],
              legend_kwds={'orientation': "horizontal",
                           'label':'(%)'})
polygons.plot(ax=ax[6],column='prros_over_prextreme', cmap=cm.cm.dense_r,
              legend=True,
              cax=cax[6],
              legend_kwds={'orientation': "horizontal",
                           'label':'(%)'})

ax[0].set_yticks(np.arange(-33,-38,-1))
ax[0].set_yticklabels(["33°S","34°S","35°S","36°S","37°S"])
ax[0].tick_params(axis="y",labelsize=14)

#plt.savefig('plots/ROS_STATS_final_basins.pdf',dpi=150,bbox_inches='tight')