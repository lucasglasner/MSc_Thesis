#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 11:14:35 2021

@author: lucas

# =============================================================================
# Comparison of datasets in relationship with the best representation of 
# real snow limit and freezing level on RioMaipoEnElManzano basin.
# =============================================================================


"""
# %%


from taylorDiagram import TaylorDiagram, test1, test2
import xarray as xr
import numpy as np
from scipy.fftpack import rfft, irfft, fftfreq, fft
from scipy.stats import linregress
import matplotlib.pyplot as plt
import scipy.stats as st
import pandas as pd

import sys

sys.path.append('scripts/')
# %%

# =============================================================================
# Load data
# =============================================================================


freezinglevel = pd.read_csv("datos/isotermas0_maipomanzano.csv", index_col=0)
freezinglevel.index = pd.to_datetime(freezinglevel.index)

snowlimits = pd.read_csv("datos/snowlimits_maipomanzano.csv", index_col=0)
snowlimits.index = pd.to_datetime(snowlimits.index)

snowcovers = pd.read_csv('datos/snowcovers_maipomanzano.csv', index_col=0)
snowcovers.index = pd.to_datetime(snowcovers.index)

# %%

# =============================================================================
# Compute mean annual cycles
# =============================================================================

FL_sc = freezinglevel.resample("d").mean()
FL_sc = FL_sc.groupby(FL_sc.index.month).mean()

SL_sc = snowlimits.groupby(snowlimits.index.month).mean()

SC_sc = snowcovers.groupby(snowcovers.index.month).mean()


# %%

# =============================================================================
# Group data by seasons
# =============================================================================

summer = [12, 1, 2]
autumn = [3, 4, 5]
winter = [6, 7, 8]
spring = [9, 10, 11]

FL_summer = freezinglevel[[x in summer for x in freezinglevel.index.month]]
FL_autumn = freezinglevel[[x in autumn for x in freezinglevel.index.month]]
FL_winter = freezinglevel[[x in winter for x in freezinglevel.index.month]]
FL_spring = freezinglevel[[x in spring for x in freezinglevel.index.month]]

SL_summer = snowlimits[[x in summer for x in snowlimits.index.month]]
SL_autumn = snowlimits[[x in autumn for x in snowlimits.index.month]]
SL_winter = snowlimits[[x in winter for x in snowlimits.index.month]]
SL_spring = snowlimits[[x in spring for x in snowlimits.index.month]]

SC_summer = snowcovers[[x in summer for x in snowcovers.index.month]]
SC_autumn = snowcovers[[x in autumn for x in snowcovers.index.month]]
SC_winter = snowcovers[[x in winter for x in snowcovers.index.month]]
SC_spring = snowcovers[[x in spring for x in snowcovers.index.month]]


# %%

fig = plt.figure(figsize=(12, 8))
fig.tight_layout()
ax0 = fig.add_subplot(231)
ax1 = fig.add_subplot(232)
ax2 = fig.add_subplot(233)

ref0 = freezinglevel["STODOMINGO"].copy()
ref1 = snowlimits["IANIGLA"].copy()
ref2 = snowcovers["IANIGLA"].copy()

freezinglevel.drop("STODOMINGO", axis=1, inplace=True)
snowlimits.drop("IANIGLA", axis=1, inplace=True)
snowcovers.drop("IANIGLA", axis=1, inplace=True)


td0 = TaylorDiagram(ref0.std(), fig=fig, rect=231,
                    srange=(0, 2.0), label="STODOMINGO")
td1 = TaylorDiagram(ref1.std(), fig=fig, rect=232,
                    srange=(0, 2.5), label="IANIGLA")
td2 = TaylorDiagram(ref2.std(), fig=fig, rect=233,
                    srange=(0, 1.5), label="IANIGLA")

std0 = [freezinglevel[m].std() for m in freezinglevel]
std1 = [snowlimits[m].std() for m in snowlimits]
std2 = [snowcovers[m].std() for m in snowcovers]

corr0 = []
corr1 = []
corr2 = []

for n, var in enumerate([freezinglevel, snowlimits, snowcovers]):
    for i in range(var.shape[1]):
        index = eval("ref"+str(n)).dropna().index
        index = index.intersection(var.iloc[:, i].dropna().index)
        x = eval("ref"+str(n)).reindex(index).dropna()
        y = var.iloc[:, i].reindex(index).dropna()
        eval("corr"+str(n)).append(st.pearsonr(x.values, y.values)[0])

cmaps = [plt.cm.plasma(np.linspace(0, 1, freezinglevel.shape[1])),
         plt.cm.nipy_spectral(np.linspace(0, 1, snowlimits.shape[1])),
         plt.cm.viridis(np.linspace(0, 1, snowcovers.shape[1]))]
for n, var in enumerate([freezinglevel, snowlimits, snowcovers]):
    colors = cmaps[n]
    for i, (std, corr) in enumerate(zip(eval("std"+str(n)), eval("corr"+str(n)))):
        td = eval("td"+str(n))
        td.add_sample(std, corr,
                      marker='o', ms=5, ls='',
                      mfc=colors[i], mec="k",
                      label=var.columns[i])

for td in [td0, td1, td2]:
    td.add_grid(ls=":")
    contours = td.add_contours(colors='0.5', levels=5)
    plt.clabel(contours, inline=1, fontsize=10, fmt='%.0f')

for ax in [ax0, ax1, ax2]:
    ax.axis("off")


fig.legend(td0.samplePoints,
           [p.get_label() for p in td0.samplePoints],
           prop=dict(size='small'), loc=(0.06, 0.2),
           ncol=2)

fig.legend(td1.samplePoints,
           [p.get_label() for p in td1.samplePoints],
           prop=dict(size='small'), loc=(0.4, 0.2),
           ncol=2)

fig.legend(td2.samplePoints,
           [p.get_label() for p in td2.samplePoints],
           prop=dict(size='small'), loc=(0.68, 0.2),
           ncol=2)

plt.savefig("ej.pdf", dpi=150, bbox_inches="tight")

# %%


# fig, ax = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(10, 5))
# fig.tight_layout(pad=3)
# fig.text(0, 0.5, "Snow Limit height (m)",
#          ha="center", va="center", rotation=90)

# # box1 = []
# box2 = []

# for axis in ax.ravel()[:-1]:
#     bbox = axis.get_position()
#     box2.append(fig.add_axes([bbox.xmax, bbox.ymin, 0.05, bbox.ymax-bbox.ymin],
#                              sharey=axis))


# var = [snowlimits.iloc[:, i] for i in range(6)]
# names = ["MODIS_H20", "MODIS_H50", "MODIS_H80", "DGF", "HYPSOMETRY\nIANIGLA"]
# colors = plt.cm.get_cmap("tab10", len(names))(np.linspace(0, 1, len(names)))
# ax = ax.ravel()
# # ax[4].plot([],[],ls=":",color="k",label="Linear\nRegression")
# ax[5].plot([], [], ls=":", color="r", label="$y\sim x$")
# linmodels = []
# for i in range(len(ax)-1):
#     x, y = var[5].dropna(), var[i].dropna()
#     x = x.reindex(y.index)
#     y = y.reindex(x.index)
#     ax[i].scatter(x,
#                   y, edgecolor="k", alpha=0.6, color=colors[i], zorder=3,
#                   lw=0.4, label=names[i])
#     m = linregress(x, y)
#     ax[i].plot(np.arange(800, 7.5e3),
#                np.arange(800, 7.5e3)*m.slope+m.intercept, ls=":", color="k")
#     t = ax[i].text(x=0.5, y=0.02,
#                    s="y~" +
#                    "{:.2f}".format(m.slope)+"x"+"{:.2f}".format(m.intercept) +
#                    "\n$R^2$: "+"{:.1%}".format(m.rvalue**2),
#                    transform=ax[i].transAxes)
#     ax[i].plot([1e2, 7.5e3], [1e2, 7.5e3], color="red", ls=":")
#     ax[i].set_xlim(0, 7.5e3)
#     ax[i].set_ylim(0, 7.5e3)

#     ax[i].grid(True, ls=":")

#     box2[i].axis("off")
#     box2[i].boxplot(y, vert=True)
#     linmodels.append(m)
#     ax[i].legend(frameon=False, loc=(0, 1))


# ax[5].axis("off")
# ax[5].legend(frameon=False, loc=(0.1, 0.8), ncol=1)


# ax[4].set_xlabel("Snow Limit IANIGLA (m)")

# plt.savefig("plots/maipomanzano/datasetcomparison/scatterplots_snowlimitS.pdf",
#             dpi=150, bbox_inches="tight")
