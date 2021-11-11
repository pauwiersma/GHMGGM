# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 13:31:16 2021

@author: Internet
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# import scipy as sc
import os 
from os.path import join 
# from mpl_toolkits.basemap import Basemap
import netCDF4 as nc4
import glob
# import descartes
import geopandas as gp
import rasterio as rio
import datetime as dt
import xarray as xr
from geocube.api.core import make_geocube
from geocube.rasterize import rasterize_points_griddata, rasterize_points_radial,rasterize_image
import json
from shapely.geometry import box,mapping
from functools import partial
from rasterio.enums import MergeAlg
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings
warnings.filterwarnings('ignore')
import time
from matplotlib.lines import Line2D
#%%

RUN_DIR = r'd:\Documents\Master vakken\Thesis\Code'
os.chdir(RUN_DIR)
NC_DIR =r'd:\Documents\Master vakken\Thesis\Code\Files\glaciers_nc'
FIG_DIR             =join(RUN_DIR,'Figures')


basin_info = pd.read_csv(join(
    RUN_DIR,'Files','basin_info_45min.csv'),index_col = 0)

### 25 basins used in the paper
# BASIN_NAMES = ['RHONE']
BASIN_NAMES = basin_info[basin_info['suitable']=='y'].index

Basins={}
for B in BASIN_NAMES: 
    Basins[B]={}
    for k in basin_info.keys():
        Basins[B][k]=basin_info.loc[B,k]
    Basins[B]['shp'] = gp.read_file(join(RUN_DIR,'Files','basin_geojsons',B+'.geojson'))
    
    

#%% Create worldmap
# OF30 = pd.read_csv(join(r'd:\Documents\Master vakken\Thesis\Code\Output',
#                         'OF29.csv'),index_col=0)



figsize=10
shapex=0.6
f1 = plt.figure(figsize=(figsize,shapex*figsize))

oceanalpha=0.00
ax1 = f1.add_axes(projection=ccrs.PlateCarree(),
                     rect = [0,0,1,1])
ax1.set_extent((-120,118,-60,85),crs=ccrs.PlateCarree())
ax1.coastlines()
ax1.add_feature(cfeature.LAND,color='green',alpha=0.00)
ax1.add_feature(cfeature.OCEAN,color='blue',alpha=oceanalpha)

#NEW-ZEALAND subaxes
ax2 = f1.add_axes(projection=ccrs.PlateCarree(),
                     rect=[0.762,0.00,0.3,0.4])
ax2.set_extent((160,180,-60,-30),crs=ccrs.PlateCarree())
ax2.coastlines()
ax2.add_feature(cfeature.LAND,color='green',alpha=00)
ax2.add_feature(cfeature.OCEAN,color='blue',alpha=oceanalpha)


#North-America subaxes
ax3 = f1.add_axes(projection=ccrs.PlateCarree(),
             rect=[0.007,0.519,0.39,0.51])
ax3.set_extent((-179,-100,30,85),crs=ccrs.PlateCarree())
ax3.coastlines()
ax3.add_feature(cfeature.LAND,color='green',alpha=00)
ax3.add_feature(cfeature.OCEAN,color='blue',alpha=oceanalpha)

#Iceland subaxes
ax4 = f1.add_axes(projection=ccrs.PlateCarree(),
      rect=[0.470,0.0,0.4,0.3])
ax4.set_extent((-28,-10,60,70),crs=ccrs.PlateCarree())
ax4.coastlines()
ax4.add_feature(cfeature.LAND,color='green',alpha=00)
ax4.add_feature(cfeature.OCEAN,color='blue',alpha=oceanalpha)

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

arrowxy = {'AMAZON':[-5,0],
       'IRRAWADDY':[2,7],
       'MACKENZIE':[-4,0],
       'YUKON':[1,6],
       'ALSEK':[-12,-7],
       'CLUTHA':[3,-3],
       'COLUMBIA':[0,-8],
       'COPPER':[-11,-5],
       'DRAMSELV':[-14,4],
       'FRASER':[-12,-9],
       'GLOMA':[-5,7],
       'KUSKOKWIM':[-16,-8],
       'NASS':[-17,-12],
       'NEGRO':[0,6],
       'OELFUSA':[-3,-1],
       'RHINE':[6,2],
       'RHONE':[-5,-10],
       'SKAGIT':[-10,-9],
       'SKEENA':[-19,-14],
       'STIKINE':[-17,-10],
       'SUSITNA':[-6,8],
       'TAKU':[-15,-8.5],
       'THJORSA':[-1,-2],
       'OB'     :[-3,0],
       'DANUBE': [8,5]}


cmap   = plt.cm.get_cmap('Blues')
vmin =0
import matplotlib.colors as mcolors


norm = mcolors.Normalize(vmin=vmin, vmax=20)



for B in BASIN_NAMES:
    if B=='CLUTHA':
        ax=ax2
    elif (Basins[B]['center_lon']<-50)&(Basins[B]['center_lat']>0):
        ax=ax3
    elif B in ['THJORSA','OELFUSA']:
        ax=ax4
        ax1.add_geometries(Basins[B]['shp'].geometry,
                crs=ccrs.PlateCarree(),
                facecolor=cmap(norm(Basins[B]['glac_degree'])))
# ax1.add_geometries(stack_shapefiles[i].geometry,
#                 crs=ccrs.PlateCarree(),
#                 alpha=0.5,)
        ax1.add_geometries(Basins[B]['shp'].geometry,
                            crs=ccrs.PlateCarree(),
                            alpha=0.5,
                            facecolor=cmap(norm(Basins[B]['glac_degree'])),edgecolor='black')
    else: ax=ax1
    
    xy = (Basins[B]['center_lon'],Basins[B]['center_lat'])
    # stack_shapefiles[i].plot(ax=ax1,edgecolor='black')
    # ax1.plot(stack_shapefiles[i].geometry)
    ax.add_geometries(Basins[B]['shp'].geometry,
                        crs=ccrs.PlateCarree(),
                        facecolor=cmap(norm(Basins[B]['glac_degree'])))
        # ax1.add_geometries(stack_shapefiles[i].geometry,
        #                 crs=ccrs.PlateCarree(),
        #                 alpha=0.5,)
    ax.add_geometries(Basins[B]['shp'].geometry,
                        crs=ccrs.PlateCarree(),
                        alpha=0.5,
                        facecolor=cmap(norm(Basins[B]['glac_degree'])),edgecolor='black')
    if B in ['AMAZON','OB','MACKENZIE']:
        arrows = None
    else:
        arrows = dict(facecolor='black',
                            width=0.001,headlength=0.02,
                            headwidth=0.02,alpha=0.7)
    ax.annotate(B.title(),xy,
                 xy+np.array(arrowxy[B])*1.3,
            arrowprops=arrows)

   

# ax1.set_title('25 basins with observations')

ax5 = f1.add_axes([-0.25,0.05,0.3,0.4])
ax5.set_visible(False)
im = ax5.imshow(np.array([[vmin,20]]),cmap=cmap)
bar =plt.colorbar(im,fraction=0.035,label='Glaciation degree [%]')

# red_dot = [Line2D([0],[0],marker='o',linestyle='None',
#            markerfacecolor='red',markeredgecolor='red',
#            markersize=10,label='Glaciers')]
# plt.legend(handles = red_dot)


save_at = join(RUN_DIR,r'Figures\Basin_maps','worldmap25_small.svg')
plt.savefig(save_at,bbox_inches='tight',format = 'svg')
plt.show()