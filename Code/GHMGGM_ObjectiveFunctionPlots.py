#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 10:00:50 2020

@author: Pau Wiersma

This script:
    Loads the objective functions and plots them 
        ND
        Ratio
        RRD
        Benchmark efficiency
        Flow duration curve efficiency
        Mass balance efficiency
        Standalone GHM evaluation
            NSE
            Calendar benchmark efficiency (Schaefli&Gupta2007)
    

To do:

Done:
    

Files needed:

    

Directories needed in run_dir:
    Files
        glacier_dailybasinsum
        glaciers_nc
    Figures
    Output (where hydrographs are )




"""
import os
from os.path import join 
import subprocess
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime
import pandas as pd
import xarray as xr
import hydroeval as he
import glob
import scipy.stats as stats

from matplotlib import rcParams, cycler
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

#%% INPUTS
RUN_DIR = r'd:\Documents\Master vakken\Thesis\Code'
os.chdir(RUN_DIR)
NC_DIR =r'd:\Documents\Master vakken\Thesis\Code\Files\glaciers_nc'
FIG_DIR             =join(RUN_DIR,'Figures')


basin_info = pd.read_csv(join(
    RUN_DIR,'Files','GHMGGM_basin_info.csv'),index_col = 0)

### 25 basins used in the paper
# BASIN_NAMES = ['RHONE']
# BASIN_NAMES = basin_info[basin_info['suitable']=='y'].index
# BASIN_NAMES = ['AMAZON','IRRAWADDY','MACKENZIE',
#                 'OB','YUKON','ALSEK', 'CLUTHA', 'COLUMBIA', 'COPPER', 'DANUBE', 'DRAMSELV',
#         'FRASER', 'GLOMA',  
#         'KUSKOKWIM', 'NASS', 'NEGRO', 'OELFUSA',
#         'RHINE', 'RHONE',  'SKAGIT', 'SKEENA', 'STIKINE','SUSITNA',
#         'TAKU', 'THJORSA'] #minus Indus, Kalixaelven, Nelson, Joekulsa, Santa Cruz, Lule
# BASIN_NAMES = ['ALSEK','COLUMBIA','RHONE','OELFUSA','MACKENZIE','RHINE']
### Basins with routing problems
# BASIN_NAMES = ['JOEKULSA','NELSON','SANTA_CRUZ','LULE','KALIXAELVEN]
###Large Basins
# BASIN_NAMES  = ['AMAZON','IRRAWADDY','MACKENZIE','OB','YUKON'] 
###Basins with a second station
# BASIN_NAMES = ['CLUTHA','COLUMBIA','RHINE','SUSITNA','DANUBE']

MODEL_SETUPS        = ['2','0','1']
MODEL_NAMES         ={'0':'Modelled (Benchmark)',
                      '1':'Modelled (Bare)',
                      '2':'Modelled (Coupled)'}


MAIN_PLOT           =['s0','s1','s2'] #Which model settings to plot
ONLY_OBSYEARS       =True  #Consider only years for which GRDC is available
PLOT_BOOL           = True
SAVE_FIGS           = False
CALENDAR_DAY        =False #Calculate Caldendar day bencmark (Schaefli&Gupta2007) 

NO_OF_BASINS      = 25

GHM_NAME          ='PCRG'
GG_GCM              ='HadGEM2-ES'
GG_rcp              ='rcp26'
Fromyear            =2001
Untilyear           =2012

#%% Load OFs
OF_sorted = pd.read_csv(join(RUN_DIR,'Output','OF_sorted'+str(NO_OF_BASINS)+'.csv'),index_col=0)

ratio_monthly= pd.read_csv(join(RUN_DIR,'Output','ratio_'+str(NO_OF_BASINS)+'.csv'),index_col=0)

normdif_means = pd.read_csv(join(RUN_DIR,'Output','ND_'+str(NO_OF_BASINS)+'.csv'),index_col=0)


RRD_means = pd.read_csv(join(RUN_DIR,'Output','RRD_means.csv'),index_col=0).transpose()



#%% Benchmark efficiency plots
OF_sorted_title_index = [OF_sorted.index[i].title() for i in range(len(OF_sorted.index))]

f1,ax=plt.subplots(1,3,figsize=(10,10),sharey=True)
f1.subplots_adjust(wspace=0.02)    
OF_sorted['TFBE'] = 1-(OF_sorted['MB2']/OF_sorted['MB0'])

OFS = ['BE','FDBE','TFBE']
for i in range(3):
    ax[i].plot(OF_sorted[OFS[i]],OF_sorted_title_index,marker='o',
             linestyle='None',color='orange',markeredgecolor='black')
    plt.xticks(rotation='horizontal')
    ax[i].axvline(0,color='black',linestyle='--')
    ax[i].set_xlim(right=1,left=-1.5)

    xmax = 1
    zz  = np.array([OF_sorted['GF99'].values]).transpose()
    cmap = plt.cm.Blues
    im = ax[i].imshow(zz,cmap=cmap,aspect='auto',
                    extent = (*ax[i].get_xlim(),len(OF_sorted)-0.5,-0.5))
    ax[i].set_xlabel(OFS[i]+' [-]')
    ax[i].grid()
    ax[i].set_xticks([-1,-0.5,0,0.5,1])
    
    for j in range(len(OF_sorted[OFS[i]])):
        OF = OF_sorted[OFS[i]]
        if OF[j]<-1.5:
            # ax[i].plot(-1.5,j,color='red',marker='x')
            ax[i].annotate(str(round(OF[j],2)),(-1.45,j)
                           ,bbox=dict( fc="0.9",alpha=0.6))

ax2 = f1.add_axes([0,0,0.95,1])
ax2.set_visible(False)
im = ax2.imshow(np.array([[0,1]]),cmap=cmap)
bar =plt.colorbar(im,fraction=0.02,label=r'$P_{99}$'+' glacier contribution [-]')
# f1.savefig(join(FIG_DIR,'Overallmetrics_99.svg'),format = 'svg',bbox_inches = 'tight')
# im = ax1.imshow(zz,cmap=palette,aspect='auto',vmin=0,vmax=1,
#                 extent = (*ax1.get_xlim(),len(RRD)-0.5,-0.5))

#%%Independent GHM evaluation
if CALENDAR_DAY==True:
    OF_sorted.index = [OF_sorted.index[i].title() for i in range(len(OF_sorted.index))]
    
    f1,ax=plt.subplots(1,2,figsize=(8,10),sharey=True)
    f1.subplots_adjust(wspace=0.02)    
    
    OFS = ['NSE0','BE0']
    for i in range(2):
        ax[i].plot(OF_sorted[OFS[i]],OF_sorted.index,marker='o',
                 linestyle='None',color='orange',markeredgecolor='black')
        plt.xticks(rotation='horizontal')
        ax[i].axvline(0,color='black',linestyle='--')
        ax[i].set_xlim(right=1,left=-1.5)
        
        xmax = 1
        zz  = np.array([OF_sorted['GF99'].values]).transpose()
        cmap = plt.cm.Blues
        im = ax[i].imshow(zz,cmap=cmap,aspect='auto',
                        extent = (*ax[i].get_xlim(),len(OF_sorted)-0.5,-0.5))
        if i==0:
            ax[i].set_xlabel(r'$NSE_{Benchmark}\/ [-]$')
        elif i==1:
            ax[i].set_xlabel(r'$CBE_{Benchmark}\/ [-]$')
        ax[i].grid()
        ax[i].set_xticks([-1,-0.5,0,0.5,1])
        
        for j in range(len(OF_sorted[OFS[i]])):
            OF = OF_sorted[OFS[i]]
            if OF[j]<-1.5:
                # ax[i].plot(-1.5,j,color='red',marker='x')
                ax[i].annotate(str(round(OF[j],2)),(-1.45,j)
                               ,bbox=dict( fc="0.9",alpha=0.6))    
    
    ax2 = f1.add_axes([0.03,0.0,0.95,1])
    ax2.set_visible(False)
    im = ax2.imshow(np.array([[0,1]]),cmap=cmap)
    bar =plt.colorbar(im,fraction=0.03,label=r'$Q_{99}$'+' glacier contribution [-]')
    # f1.savefig(join(FIG_DIR,'NSE0_BBE0.svg'),format='svg',bbox_inches = 'tight')
#%% ND plot for all basins


import matplotlib.colors
pubu = plt.get_cmap('PuBu',200)
cmap = ListedColormap(pubu(range(200))[40:])


vmin =0
norm = matplotlib.colors.Normalize(vmin=vmin, vmax=1)

OF_sorted_ascending= OF_sorted.sort_values(by=['GF99'],axis=0,ascending=True)

rcParams['axes.prop_cycle'] = cycler(color=cmap(norm(OF_sorted_ascending.GF99.values)))

# xticks = normdif_means.plot(linewidth=1.2,alpha =0,legend=False).get_xticks()

f1,ax1 = plt.subplots(figsize=(12,8))
im = ax1.plot(normdif_means.index,normdif_means.values)

quantile_alpha = 1
normdif_means.quantile(0.75,axis=1).plot(ax=ax1,color='red',
                                        linewidth=1.5,label='75',
                                        linestyle='--',alpha=quantile_alpha)
normdif_means.mean(axis=1).plot(ax=ax1,color='black',linewidth=2,
                                label='50',alpha=quantile_alpha)
normdif_means.quantile(0.25,axis=1).plot(ax=ax1,color='red',
                                          linewidth=1.5,label='25',
                                          linestyle='--',alpha=quantile_alpha)

ax1.set_ylabel('Coupled model - Benchmark\n Normalized difference [-] ')

ax1.axhline(0,color='k',linestyle='--')
monthcombis = ['January','February','March','April','May','June','July','August','September','October','November','December','']
# plt.xticks(xticks[:-1],monthcombis,rotation = 45)
plt.xticks(np.linspace(0,365,13),monthcombis,rotation=45)
ax1.grid()
ax1.set_xlim(normdif_means.index[0],normdif_means.index[-1])
ax2 = f1.add_axes([0,0,0.93,1])
ax2.set_visible(False)
im = ax2.imshow(np.array([[vmin,1]]),cmap=cmap)
bar =plt.colorbar(im,fraction=0.015,label=r'$P_{99}$'+' glacier contribution [-]')
# bar =plt.colorbar(im,fraction=0.015,label=r'FQ99 [-]')
ax1.legend(loc='upper left',title='Percentiles')

ax1.set_ylim(-0.8125,0.9957)
if SAVE_FIGS==True:
    f1.savefig(join(FIG_DIR,'ND_longnames_percentiles_P99.svg'),format='svg',bbox_inches = 'tight')


#%%



#%%RRD-plots
# With cax and R2

RRD=RRD_means
RRD['GF99']=OF_sorted.sort_values(by='GF99',ascending=True)['GF99']
RRD = RRD.sort_values(by='GF99',ascending=False)
RRD.index = [RRD.index[i].title() for i in range(len(RRD.index))]

rscores = []
for i in range(5,10):
    a,b,r,p,std = stats.linregress(x=RRD[i].values,
                                    y=RRD['GF99'].values)
    rscores.append(r**2)
    
    
markers=['o','v','^','s','D']
colors = ['tab:cyan','tab:orange','tab:green','tab:red','tab:purple']
f1,ax1 = plt.subplots(figsize=(6,10))
alpha = [1,1,1,1,1]
for month in range(5,10):
    ax1.plot(RRD[month],RRD.index,color=colors[month-5],
              linestyle='-.',marker=markers[month-5],
              label=monthcombis[month-1]+'\n'+r'$R^2$ = '+str(round(rscores[month-5],2))+'\n',
              markeredgecolor='black',markeredgewidth=0.5,markersize=7,
              alpha=alpha[month-5])
plt.xticks(rotation = 0)

ax1.grid()
ax1.set_xlim(right=1,left=-1)
ax1.set_xlabel('RRD [-]')
plt.legend(loc = (1.01,0))
xmax = 1
zz  = np.array([RRD['GF99'].values]).transpose()
palette = plt.cm.Blues
import matplotlib.colors
norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
im = ax1.imshow(zz,cmap=palette,aspect='auto',vmin=0,vmax=1,
                extent = (*ax1.get_xlim(),len(RRD)-0.5,-0.5))
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
axins = inset_axes(ax1,
                   width="5%",  # width = 5% of parent_bbox width
                   height="50%",  # height : 50%
                    loc='right',
                    bbox_to_anchor=(0.17, 0.17, 1, 1),
                    bbox_transform=ax1.transAxes,
                    # borderpad=0,
                   )
cbar = f1.colorbar(im,cax=axins,ticks=np.linspace(0,1,6),fraction=0.04
                    ,label=r'$P_{99}$'+' glacier contribution [-]')
# cbar = f1.colorbar(im,cax=axins,ticks=np.linspace(0,1,6),fraction=0.04
                    # ,label=r'FQ99 [-]')

ax1.axvline(0,linestyle='--',color='black',alpha=0.6)

if SAVE_FIGS==True:
    f1.savefig(join(FIG_DIR,'RRD_PQ99.svg'),format='svg',bbox_inches = 'tight')




#%% Calculate percentages normdif_stack = pd.concat(normdiflist,axis=1)


import matplotlib.colors
pubu = plt.get_cmap('PuBu',200)
cmap = ListedColormap(pubu(range(200))[40:])


vmin =0
norm = matplotlib.colors.Normalize(vmin=vmin, vmax=1)

OF_sorted_ascending= OF_sorted.sort_values(by=['GF99'],axis=0,ascending=True)

rcParams['axes.prop_cycle'] = cycler(color=cmap(norm(OF_sorted_ascending.GF99.values)))

# xticks = normdif_means.plot(linewidth=1.2,alpha =0,legend=False).get_xticks()




f1,(ax1,ax2) = plt.subplots(2,1,figsize=(12,10),sharex=True,
                            gridspec_kw={
                           'width_ratios': [1],
                           'height_ratios': [2.5, 1]})
plt.subplots_adjust(hspace=0.05)
im = ax1.plot(normdif_means.index,normdif_means.values,alpha=1)

normdif_means.quantile(0.75,axis=1).plot(ax=ax1,color='red',
                                        linewidth=1.5,label='75',
                                        linestyle='--',alpha=1)
normdif_means.mean(axis=1).plot(ax=ax1,color='black',linewidth=2,
                                label='50',alpha=1)
normdif_means.quantile(0.25,axis=1).plot(ax=ax1,color='red',
                                          linewidth=1.5,label='25',
                                          linestyle='--',alpha=1)

ax1.set_ylabel('Coupled model - Benchmark \n Normalized difference  [-]')
ax1.legend(loc='upper left',title='Percentiles')
ax1.axhline(0,color='k',linestyle='--')
# monthcombis = ['January','February','March','April','May','June','July','August','September','October','November','December']

# plt.xticks(xticks[:-1],monthcombis,rotation = 45)
plt.xticks(np.linspace(0,365,13),monthcombis,rotation=30)


ax1.grid()
ax1.set_xlim(normdif_means.index[0],normdif_means.index[-1])
ax1.set_ylim(-0.8125,0.9957)
ax1.text(340,0.85,'a)',size=15)
# xtickss = ax1.get_xticks()
# xtickss = np.append(xtickss,xtickss[-1]+30)
summer ={}
for Basin_name in OF_sorted_ascending.index:
    ratio= np.append(ratio_monthly[Basin_name].values,ratio_monthly[Basin_name].values[0])
    ax2.plot(np.linspace(0,365,13)[:13],ratio,alpha=1)
    summer[Basin_name]=(np.mean(ratio_monthly[6:8]))
ax2.axhline(1,linestyle='--',color='black')
# plt.xticks(np.linspace(0,365,13),monthcombis,rotation=30)
ax2.grid()
ax2.set_ylabel('Coupled model to Benchmark \n Ratio [%]')
# ax2.set_ylabel('Ratio [-]')

# ax2.text(xticks[-1]-10,3.2,'b)',size=15)
ax2.text(340,220,'b)',size=15)
# ax2.set_xticks(ax2.get_xticks(),monthcombis)
# newticks=ax2.get_xticks()-10
# ax2.set_xticks(newticks)
# ax2.set_xticklabels(monthcombis,rotation=20)



ax3 = f1.add_axes([0,0,0.93,1])
ax3.set_visible(False)
im = ax3.imshow(np.array([[vmin,1]]),cmap=cmap)
bar =plt.colorbar(im,fraction=0.015,label=r'$P_{99}$'+' glacier contribution [-]')
# bar =plt.colorbar(im,fraction=0.015,label=r'FQ99 [-]')

if SAVE_FIGS==True:
    f1.savefig(join(FIG_DIR,'ND+ratio_P99_longlabels_percentage.svg'),format='svg',bbox_inches = 'tight')


# xlabels=ax1.get_xticklabels()
# ax1.set_xticks(xticks-5,)

# ax1.set_title('Thjorsa')
# 
# bar.ax.set_yticklabels(['0','1','2','3','4','5'])
# plt.xticks(n)
# ax1.set_xlabel('Months')
# ax1.set_xticks()


