# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 13:31:16 2021

@author: Pau Wiersma

Scripts accompanying PCR-GLOBWB 2 and GloGEM coupling
Auxiliary: Plotting of the worldmap in the paper (figure 1)

This script loads the basin shapefiles and plots them on a world map 
with the glacierization degree in blue hue

Files needed: 
    GHMGGM_basin_info.csv
    basin_geojsons 
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os 
from os.path import join 
import geopandas as gp
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings
warnings.filterwarnings('ignore')
#%%

RUN_DIR = r'd:\Documents\Master vakken\Thesis\Code'
os.chdir(RUN_DIR)
NC_DIR =r'd:\Documents\Master vakken\Thesis\Code\Files\glaciers_nc'
FIG_DIR             =join(RUN_DIR,'Figures')


basin_info = pd.read_csv(join(
    RUN_DIR,'Files','GHMGGM_basin_info.csv'),index_col = 0)

### 25 basins used in the paper
# BASIN_NAMES = ['RHONE']
BASIN_NAMES = basin_info[basin_info['suitable']=='y'].index

Basins={}
for B in BASIN_NAMES: 
    Basins[B]={}
    for k in basin_info.keys():
        Basins[B][k]=basin_info.loc[B,k]
    Basins[B]['shp'] = gp.read_file(join(RUN_DIR,'Files','basin_geojsons',B+'.geojson'))
    
SAVE_FIG = False
    

#%% Create worldmap

FIGSIZE=10
SHAPEX=0.6
f1 = plt.figure(figsize=(FIGSIZE,SHAPEX*FIGSIZE))

OCEANALPHA=0.02
ax1 = f1.add_axes(projection=ccrs.PlateCarree(),
                     rect = [0,0,1,1])
ax1.set_extent((-120,118,-60,85),crs=ccrs.PlateCarree())
ax1.coastlines()
ax1.add_feature(cfeature.LAND,color='green',alpha=0.00)
# ax1.add_feature(cfeature.OCEAN,color='blue',alpha=OCEANALPHA)

#NEW-ZEALAND subaxes
ax2 = f1.add_axes(projection=ccrs.PlateCarree(),
                     rect=[0.762,0.00,0.3,0.4])
ax2.set_extent((160,180,-60,-30),crs=ccrs.PlateCarree())
ax2.coastlines()
ax2.add_feature(cfeature.LAND,color='green',alpha=00)
# ax2.add_feature(cfeature.OCEAN,color='blue',alpha=OCEANALPHA)


#North-America subaxes
ax3 = f1.add_axes(projection=ccrs.PlateCarree(),
             rect=[0.007,0.519,0.39,0.51])
ax3.set_extent((-179,-100,30,85),crs=ccrs.PlateCarree())
ax3.coastlines()
ax3.add_feature(cfeature.LAND,color='green',alpha=00)
# ax3.add_feature(cfeature.OCEAN,color='blue',alpha=OCEANALPHA)

#Iceland subaxes
ax4 = f1.add_axes(projection=ccrs.PlateCarree(),
      rect=[0.470,0.0,0.4,0.3])
ax4.set_extent((-28,-10,60,70),crs=ccrs.PlateCarree())
ax4.coastlines()
ax4.add_feature(cfeature.LAND,color='green',alpha=00)
# ax4.add_feature(cfeature.OCEAN,color='blue',alpha=OCEANALPHA)

#Large
# arrowxy = {'AMAZON':[-3,-3],
#            'IRRAWADDY':[5,5],
#            'MACKENZIE':[2,0],
#            'YUKON':[-4,10],
#            'ALSEK':[-8,-3],
#            'CLUTHA':[3,-3],
#            'COLUMBIA':[0,-8],
#            'COPPER':[-9,-5],
#            'DRAMSELV':[-12,4],
#            'FRASER':[-12,-7],
#            'GLOMA':[-5,7],
#            'KUSKOKWIM':[-20,-5],
#            'NASS':[-12,-7],
#            'NEGRO':[0,6],
#            'OELFUSA':[-3,-1],
#            'RHINE':[6,2],
#            'RHONE':[-3,-6],
#            'SKAGIT':[-10,-7],
#            'SKEENA':[-8,-5.5],
#            'STIKINE':[-12,-5],
#            'SUSITNA':[-12,-10],
#            'TAKU':[-10,-3.5],
#            'THJORSA':[-1,-2],
#            'OB'     :[-3,-3],
#            'DANUBE': [5,5]}

#Smaller figure
arrowxy = {'AMAZON':[-5,0],
       'IRRAWADDY':[2,7],
       'MACKENZIE':[-3,0],
       'YUKON':[-4,8],
       'ALSEK':[-12,-7],
       'CLUTHA':[3,-3],
       'COLUMBIA':[0,-8],
       'COPPER':[-9,-5],
       'DRAMSELV':[-14,4],
       'FRASER':[-12,-9],
       'GLOMA':[-5,7],
       'KUSKOKWIM':[-8,10], #-16,-3
       'NASS':[-17,-12],
       'NEGRO':[0,6],
       'OELFUSA':[-3,-1],
       'RHINE':[6,2],
       'RHONE':[-5,-10],
       'SKAGIT':[-10,-9],
       'SKEENA':[-18,-13.5], #-8,-7.5
       'STIKINE':[-17,-10],
       'SUSITNA':[-16,-8],
       'TAKU':[-15,-8.5],
       'THJORSA':[-1,-2],
       'OB'     :[0,0],
       'DANUBE': [8,5]}


cmap   = plt.cm.Blues
VMIN =0
import matplotlib.colors as mcolors


norm = mcolors.Normalize(vmin=VMIN, vmax=20)



for B in BASIN_NAMES:
    if B=='CLUTHA':
        ax=ax2
    elif (Basins[B]['center_lon']<-50)&(Basins[B]['center_lat']>0):
        ax=ax3
    elif B in ['THJORSA','OELFUSA']:
        ax=ax4
        ax1.add_geometries(Basins[B]['shp'].geometry,
                            crs=ccrs.PlateCarree(),
                            alpha=0.5,
                            facecolor=cmap(norm(Basins[B]['glacerization_degree'])),edgecolor='black')
    else: ax=ax1
    
    xy = (Basins[B]['center_lon'],Basins[B]['center_lat'])
    # stack_shapefiles[i].plot(ax=ax1,edgecolor='black')
    # ax1.plot(stack_shapefiles[i].geometry)
    ax.add_geometries(Basins[B]['shp'].geometry.buffer(0.01),
                        crs=ccrs.PlateCarree(),
                        facecolor=cmap(norm(Basins[B]['glacerization_degree'])))
    ax.add_geometries(Basins[B]['shp'].geometry.buffer(0.01),
                        crs=ccrs.PlateCarree(),
                        alpha=0.5,facecolor=cmap(norm(Basins[B]['glacerization_degree'])),
                        edgecolor='black')
    if B in ['AMAZON','OB','MACKENZIE']:
        ax.annotate(B.title(),xy,
                 xy+np.array(arrowxy[B])*1.3)
    else:
        ax.annotate(B.title(),xy,
                     xy+np.array(arrowxy[B])*1.3,
                arrowprops=dict(facecolor='black',
                                width=0.001,headlength=0.02,
                                headwidth=0.02,alpha=0.7))

ax5 = f1.add_axes([-0.25,0.05,0.3,0.4])
ax5.set_visible(False)
im = ax5.imshow(np.array([[VMIN,20]]),cmap=cmap)
bar =plt.colorbar(im,fraction=0.035,label='Glacierization degree [%]')

if SAVE_FIG:
    save_at = join(RUN_DIR,r'Figures\Basin_maps','worldmap25_small.svg')
    plt.savefig(save_at,bbox_inches='tight',format = 'svg')
    plt.show()