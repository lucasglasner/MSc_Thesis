#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 14:34:49 2022

@author: lucas
"""

import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
#%%
interval = slice("2000-01-01","2022-01-01")
pr_sanjose = pd.read_csv('datos/estaciones/pr_SanJoseMaipo.csv',dtype=str)
pr_sanjose.index = pr_sanjose.iloc[:,0]+"-"+pr_sanjose.iloc[:,1]+"-"+pr_sanjose.iloc[:,2]
# pr_sanjose.index = pd.to
pr_sanjose.index = pd.to_datetime(pr_sanjose.index)
pr_sanjose = pd.to_numeric(pr_sanjose.iloc[:,3])
pr_sanjose = pr_sanjose[interval]


pr_cr2met = xr.open_mfdataset('datos/cr2met/complete_cr2met/CR2MET_pr*')

pr_cr2met_sanjose = pr_cr2met.sel(lon=-70.35,
                                  lat=-33.64,
                                  method='nearest').pr.to_series()
pr_cr2met_sanjose = pr_cr2met_sanjose[interval]
pr_cr2met_qn = pr_cr2met.sel(lon=-70.68,
                              lat=-33.44,
                              method='nearest').pr.to_series()
pr_cr2met_qn = pr_cr2met_qn[interval]

del pr_cr2met

# pr_qn = pd.read_csv('datos/estaciones/qn/DATOSUTC_2004-2019.csv',index_col=0)
# pr_qn.index = pd.to_datetime(pr_qn.index)
# pr_qn = pr_qn[interval].iloc[:,9].resample('d').sum()

pr_qn = pd.read_csv('datos/estaciones/pr_quintanormal.csv',dtype=str)
pr_qn.index = pr_qn.iloc[:,0]+"-"+pr_qn.iloc[:,1]+"-"+pr_qn.iloc[:,2]
pr_qn.index = pd.to_datetime(pr_qn.index)
pr_qn = pd.to_numeric(pr_qn.iloc[:,3])[interval]
#%%

pr = pd.concat([pr_sanjose,pr_qn,pr_cr2met_sanjose,pr_cr2met_qn],axis=1)
pr.columns = ['sanjose','qn','cr2met_sanjose','cr2met_qn']
#%%
def tick_density(ax,which="x",every=1):
    if which=="x":
        ticks = ax.get_xlim()
        mn,mx = np.min(ticks),np.max(ticks)
        ax.set_xticks(np.arange(mn,mx,every))
    elif which=="y":
        ticks=ax.get_yticks()
        mn,mx = np.min(ticks),np.max(ticks)
        ax.set_yticks(np.arange(mn,mx,every))
    else:
        assert False    

plt.rc('font',size=18)
pr_month = pr.resample('m').sum()
pr_month = pr_month.groupby(pr_month.index.month).mean()

pr_year = pr.resample('y').sum()

fig,ax = plt.subplots(2,3,sharex='col',figsize=(12,4))

fig.tight_layout(pad=1)
ax[0,0].bar(np.arange(0,12)-0.25,pr_month.sanjose,width=0.25,ec='k')
ax[0,0].bar(np.arange(0,12),pr_month.cr2met_sanjose,width=0.25,ec='k')

ax[1,0].bar(np.arange(0,12)-0.25,pr_month.qn,width=0.25,ec='k')
ax[1,0].bar(np.arange(0,12),pr_month.cr2met_qn,width=0.25,ec='k')


ax[1,0].scatter([],[],marker='s',label='Rain Gauge',color='tab:blue')
ax[1,0].scatter([],[],marker='s',label='Reanalysis',color='tab:orange')
ax[1,0].legend(frameon=False,ncol=2,loc=(-0.2,-0.6),fontsize=15)

ax[0,1].plot(pr_year.sanjose)
ax[0,1].plot(pr_year.cr2met_sanjose)

ax[1,1].plot(pr_year.qn)
ax[1,1].plot(pr_year.cr2met_qn)


ax[0,2].scatter(pr.cr2met_sanjose.where(pr.cr2met_sanjose>3),
                pr.sanjose.where(pr_sanjose>3))
ax[1,2].scatter(pr.cr2met_qn.where(pr.cr2met_qn>3),
                pr.qn.where(pr_qn>3))

ax[0,2].plot([0,120],[0,120],ls="--",color='k')
ax[1,2].plot([0,120],[0,120],ls="--",color='k')

ax[0,2].set_xlim(0,120)
ax[0,2].set_ylim(0,120)

ax[1,2].set_xlim(0,120)
ax[1,2].set_ylim(0,120)

tick_density(ax[1,2],"x",50)
tick_density(ax[1,2],"y",50)


tick_density(ax[0,2],"y",50)


tick_density(ax[0,0],"y",50)

tick_density(ax[1,0],"y",25)


# tick_density(ax[0,0],"x",1)
ax[0,0].set_xticks(np.arange(12))
ax[0,0].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
# tick_density(ax[1,0],"x",1)

# tick_density(ax[1,1],"y",250)
ax[1,1].set_yticks([0,250,500])

ax[0,1].set_yticks([0,5e2,1e3])
ax[1,1].tick_params(axis="x",rotation=45)

tick_density(ax[0,-1],which='x',every=50)
ax[0,-1].set_xlim(0,120)
ax[0,-1].set_ylim(0,120)

ax[1,-1].set_ylim(0,120)

ax[0,0].set_ylabel('San Jose\nde Maipo')
ax[1,0].set_ylabel('Quinta\nNormal')
ax[0,0].set_title('Precipitation\n(mm/month)',loc='left')


ax[0,1].set_title('Precipitation\n(mm/year)',loc='left')


ax[0,2].set_title('Precipitation\n(mm/day)',loc='left')

ax[1,2].set_xlabel('CR2MET')

plt.savefig('plots/maipomanzano/PR_CR2MET_STUDY.pdf',dpi=150,bbox_inches='tight')
