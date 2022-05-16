
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 12:38:15 2022

@author: lucas
"""

import xarray as xr
import pandas as pd
from scipy.interpolate import interp1d
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import matplotlib as mpl
from functions import local_minimum_filter
from scipy.signal import find_peaks, butter, filtfilt
def butter_lowpass_filter(data, cutoff, fs, order):
    nyq = 0.5*fs
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

#%%
b = 'RioTenoDespuesDeJuntaConClaro'
data_daily = pd.read_csv('datos/TABLAS_DIARIAS_CRECIDAS_CC.csv',index_col=[0,1])
# data_daily = data_daily.loc['RioMaipoEnElManzano']
data_daily = data_daily.loc[b]
data_daily.index = pd.to_datetime(data_daily.index)
data_daily['ros_day'] = (data_daily['ros_area']>0.1) & (data_daily['FL']>data_daily['SL'])


data_events = pd.read_csv('datos/TABLAS_EVENTOS_CRECIDAS.csv',index_col=[0,1])
data_events = data_events.loc[b]
data_events['ros_event'] = data_events['max_ros_area']>0.1


hypso = pd.read_csv('datos/topography/basins/hypso/'+b+'_hypso.csv',
                    index_col=0)

pr = pd.read_csv('datos/estaciones/vismet/pr_'+b+'.csv',
                 index_col=0)['Valor']
pr.index = pd.to_datetime(pr.index)

q = pd.read_csv('datos/estaciones/dga/qinst_'+b+'.csv',
                index_col=0).squeeze()
q.index = pd.to_datetime(q.index)

SCA = pd.read_csv('datos/ianigla/'+b+'_SCA_s_comp.filtro_MA.3días.csv',
                  index_col=0)['SCA(%)']
SCA.index = pd.to_datetime(SCA.index)
SCA_trend = (SCA.shift(-3)-SCA.shift(1))/2
SCA = SCA.resample('h').interpolate('cubicspline')
SCA = SCA.where(SCA>0).fillna(0)
SCA = SCA.where(SCA<100).fillna(100)

H0 = pd.read_csv('datos/era5/horarias/H0_'+b+'_2000-2020.csv',
                 index_col=0).deg0l
H0.index = pd.to_datetime(H0.index)

int_func = interp1d(hypso.index,hypso.fArea)
int_func2 = interp1d(hypso.fArea,hypso.index)
SL = pd.Series(int_func2(SCA.values/100),index=SCA.index)
ap = pd.Series(int_func(H0-300).squeeze(),index=H0.index)


basin_data = pd.concat([pr,q,SCA/100,H0,SL,ap],axis=1).dropna()
basin_data.columns = ['pr','q','sca','h0','sl','pluv_area']
basin_data['ros_area'] = basin_data['sca']-(1-basin_data['pluv_area'])
basin_data['ros_area'] = basin_data['ros_area'].where(basin_data['pr']>0)
basin_data['ros_area'] = basin_data['ros_area'].where(basin_data['ros_area']>0.1)
basin_data['SCA_trend'] = SCA_trend.resample('h').fillna('ffill')
basin_data['quickflow'] = local_minimum_filter(basin_data.q,40)[1]
basin_data['quickflow'] = basin_data.quickflow.where(basin_data.quickflow>0)

#%%
plt.rc('font',size=14)
fig,ax = plt.subplots(1,4,figsize=(12,3),sharex=False)
fig.tight_layout(pad=1.5)

# ax[0].boxplot([basin_data[basin_data.pr>0.1].h0.dropna(),
#                basin_data[(basin_data.ros_area>0.1)&(basin_data.pr>0.1)].h0.dropna(),
#                basin_data[(basin_data.ros_area>0.1)&(basin_data.pr>0.1)&(basin_data.SCA_trend<0)].h0.dropna(),
#                ],
#               positions=[1,1.5,2],sym='',
#               patch_artist=False)

m1,m2,m3 = basin_data.pr>3/24,(basin_data.pr>3/24)&(basin_data.ros_area>0.1),(basin_data.pr>3/24)&(basin_data.ros_area>0.1)&(basin_data.SCA_trend<0)
c = ["royalblue","tab:green","tab:red"]
for i in range(3):
    ax[0].boxplot(basin_data[eval("m"+str(i+1))].h0.dropna()-300,
                  positions=[i/2], sym="", patch_artist=True,
                  medianprops={'color':'k'},
                  boxprops={'facecolor':c[i]})
    ax[1].boxplot(basin_data[eval("m"+str(i+1))].pr.dropna(),
                  positions=[i/2], sym="", patch_artist=True,
                  medianprops={'color':'k'},
                  boxprops={'facecolor':c[i]})
    ax[2].boxplot(basin_data[eval("m"+str(i+1))].quickflow.dropna(),
                  positions=[i/2], sym="", patch_artist=True,
                  medianprops={'color':'k'},
                  boxprops={'facecolor':c[i]})
    ax[3].boxplot(basin_data[eval("m"+str(i+1))].sca.dropna()*100,
                  positions=[i/2], sym="", patch_artist=True,
                  medianprops={'color':'k'},
                  boxprops={'facecolor':c[i]})

ax[0].set_ylabel('Freezing level (m)')
ax[1].set_ylabel('Precipitation (mm/h)')
ax[2].set_ylabel('Direct Runoff ($m^3/s$)')
ax[3].set_ylabel('Snow Cover (%)')
for axis in ax.ravel():
    axis.set_xticklabels(['ALL','ROS','MROS'])
    
ax[0].scatter([],[],marker="s",color=c[0], label='ALL\nN=300 days')
ax[0].scatter([],[],marker="s",color=c[1], label='ROS\nN=83 days')
ax[0].scatter([],[],marker="s",color=c[2], label='MROS\nN=41 days')
ax[0].legend(frameon=False,loc=(0,1),ncol=3)

plt.savefig('plots/teno_stats.pdf',dpi=150,bbox_inches='tight')

# ax[0].boxplot(basin_data[m1].h0.dropna(),
#               positions=[0.5],
#               patch_artist=True, sym="",
#               boxprops={'facecolor':'tab:blue'})
# ax[0].boxplot(basin_data[m2].h0.dropna(),
#               positions=[1],
#               patch_artist=True, sym="",
#               boxprops={'facecolor':'red'})
# ax[0].boxplot(basin_data[m3].h0.dropna(),
#               positions=[1.5],
#               patch_artist=True, sym="",
#               boxprops={'facecolor':'red'})

# ax[1].boxplot([basin_data[basin_data.pr>1].q,
#                basin_data[(basin_data.ros_area>0.1)&(basin_data.pr>1)].q],
#               positions=[1,1.5],sym='')



# ax[0].hist(basin_data[basin_data.pr>0.1].h0,bins=30, density=False, alpha=0.6)
# ax[0].hist(basin_data[(basin_data.pr>0.1)&(basin_data.ros_area>0.1)].h0,
#            bins=15, density=False, alpha=0.6)

# # ax[1].set_xscale('log')
# ax[1].set_yscale('log')
# ax[1].hist(basin_data[basin_data.pr>1].q.dropna(),bins='auto',alpha=0.6, density=False)
# ax[1].hist(basin_data[(basin_data.pr>1)&(basin_data.ros_area>0.1)&(basin_data.SCA_trend<0)].q.dropna(),
#             bins='auto',alpha=0.6, density=False)


#%%
import scipy.stats as st
floods = basin_data.quickflow.dropna()
peaks = find_peaks(floods, height=200, prominence=100, distance=12)[0]

floods = floods.iloc[peaks].sort_values().index



ros_floods = []
nonros_floods=[]
data=[]
for f in floods:
    t1,t2 = pd.to_datetime(f)-pd.Timedelta(hours=24),pd.to_datetime(f)+pd.Timedelta(hours=48)
    t1,t2 = t1.strftime("%Y-%m-%dT%H:%M:%S"),t2.strftime("%Y-%m-%dT%H:%M:%S")
    f2 = basin_data[t1:t2]
    f2.index = np.arange(len(f2))
    f2.quickflow = f2.quickflow.fillna(0)
    
    
    
    f2['hidro'] = butter_lowpass_filter(f2.quickflow,1/6,1,1)
    f2['hidro'] = f2['hidro']/f2['pr'].sum()
    f2['hieto'] = f2.pr/f2.pr.max()
    
    f2.quickflow =  butter_lowpass_filter(f2.quickflow,1/6,1,1)
    # sos = butter(1, 1e-10, 'hp', fs=1, output='sos')

    # filtered = sosfilt(sos, f2['a_hidrograph'].dropna())
    # if len(find_peaks(f2.quickflow,threshold=20)[0])==1:
    if ((f2[f2.pr>1].ros_area.dropna().max()>0.25)&(f2[f2.pr>1].SCA_trend.mean()<0)):
        ros_floods.append(f)
    # if (f2.ros_area.dropna().mean()<0.175):
    #     nonros_floods.append(f)
    else:
        nonros_floods.append(f)
    
    data.append(f2)
data = pd.concat(data,axis=1,keys=floods)
# 

ros_floods=np.array(ros_floods)
nonros_floods = np.array(nonros_floods)[-10:]
# plt.plot(f2['a_hidrograph'].dropna())
# plt.plot(butter_lowpass_filter(f2['a_hidrograph'].dropna(),1/6,1,1))
typical_flood = data[nonros_floods].T.unstack().mean().unstack()
typical_flood_ci = data[nonros_floods].T.unstack().std().unstack()/np.sqrt(len(nonros_floods))*st.t.ppf(1-0.05/2,len(nonros_floods)-1)
typical_rosflood = data[ros_floods].T.unstack().mean().unstack()

typical_rosflood_ci = data[ros_floods].T.unstack().std().unstack()/np.sqrt(len(ros_floods))*st.t.ppf(1-0.05/2,len(ros_floods)-1)
# typical_flood.a_hidrograph.plot();typical_rosflood.a_hidrograph.plot()
# %%
plt.rc('font',size=18)
fig,ax = plt.subplots(2,1, sharex=True, figsize=(9,5), sharey='row')
ax = ax.ravel(order='F')

for axis in ax:
    axis.set_xticks(np.arange(0,72+6,6))
    axis.set_xticklabels(np.arange(0,72+6,6)-24)
    axis.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(3))
    axis.grid(True,ls=":")
    axis.grid(which='minor',axis="x",ls=":")
    
    axis.set_xlim(0,60)
    # axis.set_ylim(0,1)
ax0 = ax[0].twinx()
ax1 = ax[1].twinx()
ax0.set_ylabel('Snow cover\n(-)')
ax1.set_ylabel('Pluvial area\n(-)')


ax0.set_ylim(0.3,0.8)
# ax0.set_ylim(0,1)
# ax1.set_ylim(0,1)

c1 = 'tab:blue'
c2 = 'tab:red'

ax[0].bar(typical_flood.index,typical_flood.pr,
          width=1,alpha=0.5, color=c1)
ax[0].bar(typical_rosflood.index,typical_rosflood.pr,
          width=1,alpha=0.5,
          color=c2)

# ax[0].fill_between(typical_flood.index,
#                    typical_flood.hieto+typical_flood_ci.hieto,
#                    typical_flood.hieto-typical_flood_ci.hieto,
#                    color='grey',alpha=0.2)
# ax[0].fill_between(typical_rosflood.index,
#                    typical_rosflood.hieto+typical_rosflood_ci.hieto,
#                    typical_rosflood.hieto-typical_rosflood_ci.hieto,
#                    color='grey',alpha=0.2)
# ax[0].set_ylim(0,1)

ax0.plot(typical_flood.sca, ls="--", color=c1)
ax0.plot(typical_rosflood.sca, ls="--", color=c2)


ax[1].plot(typical_flood.hidro, color=c1)
ax[1].plot(typical_rosflood.hidro, color=c2)

ax[1].fill_between(typical_flood.index,
                   typical_flood.hidro+typical_flood_ci.hidro,
                   typical_flood.hidro-typical_flood_ci.hidro,
                   color='grey',alpha=0.2)
ax[1].fill_between(typical_rosflood.index,
                   typical_rosflood.hidro+typical_rosflood_ci.hidro,
                   typical_rosflood.hidro-typical_rosflood_ci.hidro,
                   color='grey',alpha=0.2)

ax[1].plot(typical_flood.pluv_area, color=c1, ls=":")

ax[1].plot(typical_rosflood.pluv_area, color=c2, ls=":")

ax[1].set_xlabel('Hours from peak flow')
ax[0].set_ylabel('Precipitation\n(mm/h)')
ax[1].set_ylabel('Adimensional\nrunoff\n(-)')
ax[1].set_ylim()


ax[0].scatter([],[],color=c1,label='Typical Flood\n'+'N={:d}'.format(len(nonros_floods)),
              marker="s")
ax[0].scatter([],[],color=c2,label='MROS Flood\n'+'N={:d}'.format(len(ros_floods)),marker="s")

ax[0].legend(frameon=False,loc=(0,1), ncol=2)

ax[0].text(0,0.8, "(a)",transform=ax[0].transAxes)


ax[1].text(0,0.8, "(b)",transform=ax[1].transAxes)

plt.savefig('plots/COMPUESTO_HIDROLOGICO.pdf',dpi=150,bbox_inches='tight')