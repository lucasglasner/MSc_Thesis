#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 12:44:52 2021

@author: lucas

# =============================================================================
# Compute 0°C Isotherm Height from Wyoming radiosonde data
# =============================================================================

"""


import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, splev, splrep, BSpline
import matplotlib.pyplot as plt
import pandas as pd
#%%
# =============================================================================
# Read radiosonde file
# =============================================================================

radiosonde_path = "datos/stodomingo/radiosonde_19991224_20210101.npy"
soundings       = np.load(radiosonde_path,allow_pickle=True).item()

#%%
# =============================================================================
# Compute 0°C Isotherm Height
# =============================================================================

dates   = list(soundings.keys())  #Dates in dictionary
H0      = np.nan*np.empty(len(dates))    #Where to store 0°C isotherm height

#%%
for i in range(len(dates)):
    date   = dates[i]
    radiosonde = soundings[date]
    temp   = radiosonde["temperature"].values
    height = radiosonde["height"].values
    try:
        # spline = BSpline(temp, height, 3)
        # H0[i]  = spline(0)
        # spline = splrep(height, temp)
        # H0[i]  = splev(0, spline)
        H0[i]  = interp1d(temp, height, kind="slinear", bounds_error=False, fill_value=np.nan)(0) 
    except ValueError as e:
        print(e)
        
# H0 = H0.dropna()
#%%
# =============================================================================
# Format data and export output
# =============================================================================
H0 = pd.Series(H0,index=dates).replace([np.inf,-np.inf],np.nan).dropna()
H0 = H0.reindex(pd.date_range("1999-12-24","2021-01-01",freq="12h"))
# H0 = pd.Series(H0.groupby([H0.index.year,
#                             H0.index.month,
#                             H0.index.day]).mean().ravel(),
#                 index=np.unique(H0.index.date))
# H0 = pd.Series(H0,index=dates).replace([np.inf,-np.inf],np.nan)
# # H0 = H0.where(H0<6000).interpolate(method="spline", order=3)
# # H0 = H0.resample("d").mean().dropna()
# H0 = H0.groupby
# H0 = H0.dropna()
H0.name = "H0_StoDomingo"
#%%
H0.to_csv("datos/stodomingo/isoterma0.csv")
#%%

# #Verificar
# import datetime as dt
# import seaborn as sns
# def corregir_ECR2(df):
#     """
#     Función para transformar las series de tiempo del CR2 a series
#     de tiempo trabajables facilmente con los algoritmos de pandas
#     para series de tiempo.

#     Parameters
#     ----------
#     df : Data Frame de pandas con una serie de tiempo descargada
#          del explorador climático (directo de pandas.read_csv()).

#     Returns
#     -------
#     df : Serie de Pandas con los indices de la serie objetos de 
#          tiempo.

#     """
#     index=[]
#     for i in range(len(df)):
#         index.append(dt.datetime(df["agno"][i],df[" mes"][i],df[" dia"][i]))
#     df.index = index
#     df.drop(["agno"," mes"," dia"],axis=1,inplace=True)
#     df = df.resample("d").asfreq()[" valor"].rename("serie")
#     return df

# pp = pd.read_csv("datos/estaciones/pr_quintanormal.csv")
# pp = corregir_ECR2(pp)

# data = pd.concat((H0.reindex(pp.index),pp),axis=1)
# h0_conpp = data[data["serie"]>5].iloc[:,0]
# sns.kdeplot(H0,color="brown",lw=2,label="pr>0mm")
# sns.kdeplot(h0_conpp,color="royalblue",lw=2,label="pr>5mm")
# plt.hist(H0,bins="auto",alpha=0.3,density=True,color="brown")
# plt.hist(h0_conpp,density=True,bins="auto",alpha=0.3,color="royalblue")
# plt.xlabel("H0")
# plt.legend()
# plt.xlim(0,6000)
# plt.savefig("plots/DistH0_Stgo.pdf",dpi=150,bbox_inches="tight")