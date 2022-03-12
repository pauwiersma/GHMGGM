# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 11:52:06 2021

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
import subprocess
from matplotlib.lines import Line2D
from descartes import PolygonPatch as pp
from fiona.crs import from_epsg
import fiona
from shapely.geometry import shape,mapping, Point, Polygon, MultiPolygon
import scipy.stats as stats
import rioxarray as xrio

run_dir   = r'D:\Documents\Master Vakken\Thesis\Code'
os.chdir(run_dir)
files               = join(run_dir,'Code','Files')
#%%

"""
- HH2018 annual mass balance
- daily SWE raster
    -take annual mass balance from glacfrac
    - Compare timing of snowmelt bewteen N0 and N2
- HG
    - Take difference between hg0 and hg2 and see what percentage the mass balances take in it

"""
#%% INPUTS
files               = join(run_dir,'Files')
figures             =join(run_dir,'Figures')
daily_outputs       =join(run_dir,'Output','daily_outputs')
HH_MB               = join(files,'HH2018_validation')

basin_info = pd.read_csv(join(files,'basin_info_45min.csv'),index_col = 0)
# Basin_names = basin_info[basin_info['suitable']=='y'].index


Basin_names = ['ALSEK','COLUMBIA','OELFUSA','RHONE']

model_setups        = ['0','1','2']
model_names         ={'0':'Modelled (Benchmark)',
                      '1':'Modelled (Bare)',
                      '2':'Modelled (Coupled)'}


seasonal            = False  #False for whole year, true for only summer
main_plot           =['s0','s1','s2']#,'s4']
only_obsyears       =True
save_figs           = False
new_seasonal        =False
calendar_day        =True


Model_name          ='PCRG'
GG_GCM              ='HadGEM2-ES'
GG_rcp              ='rcp26'
Fromyear            =2000
Untilyear           =2012

#%%
# Basin_names = ['ALSEK','COLUMBIA','OELFUSA','RHONE']
Basin_names = ['OELFUSA']
run_name = 'N'
SWE_diff_list = []
MB_sum_list = []
hg_diff_list = []
spilling_list = []
for Basin_name in Basin_names:
    print (Basin_name)
    if Basin_name in ['CLUTHA','NEGRO','AMAZON','SANTA_CRUZ']:
        hemisphere = 'South'
        daterange = pd.date_range(str(Fromyear)+'-04-01',
                                 str(Untilyear)+'-03-31')
    else: 
        hemisphere = 'North'
        daterange = pd.date_range(str(Fromyear)+'-10-01',
                                  str(Untilyear)+'-09-30')
    
    shp =   gp.read_file(join(files,'subbasin_geojsons',Basin_name+'.geojson'))
    Overview_path = basin_info.loc[Basin_name,'RGI_dir']+'_HadGEM2-ES_RCP26_Overview.dat'
    Overview    = pd.read_csv(join(HH_MB,Overview_path),
                            delim_whitespace=True,index_col=0)
    points              = gp.GeoSeries(gp.points_from_xy(Overview['lon'].values,Overview['lat'])
                                        ,index = Overview.index)
    
    GG_idx     = [shp.geometry.contains(point)[0] for point in points] 
    GG_points   = points[GG_idx]
    GG_area0    = Overview['Area0'][GG_idx] # in km2
        #Testplot:
    #maxarea=3223
    #Already in hydrological years
    MB_path  = basin_info.loc[Basin_name,'RGI_dir']+'_HadGEM2-ES_RCP26_Annual_balance.dat'
    MB = pd.read_csv(join(HH_MB,MB_path),delim_whitespace=True,index_col=0)
    MB =MB[[str(year) for year in np.arange(2001,2013)]]
    MB = MB[GG_idx]
    MB = MB.multiply(GG_area0,axis=0) * 1e6 #km2 to m2
    
    MB[MB==-99] = np.nan
    MB_sum =MB.sum(axis=0)
    MB_sum.index = np.arange(2001,2013)
    MB_sum_list.append(MB_sum)
    # MB_sum.plot()
    #meter, SWE is also in meter?
    
    #SWE
    SWE = xr.open_dataset(join(daily_outputs,Basin_name+'_s_d_g_N0.nc')).snow_water_equivalent
    SWE = SWE.sel(time=slice('2000','2012'))
    SWE_aletsch = SWE.isel(lat =slice(41,43),lon=60).sum(axis=1)
    
    nc_path = join(files,'glaciers_nc',Basin_name+'_2000_2012_N_R.nc')
    isglac = xr.open_dataset(nc_path).isglac
    
    
    
    cellsize =xrio.open_rasterio(join(files,'clones','cellsize05min.correct.map')).squeeze()
    cellsize = cellsize.rename({'y':'lat','x':'lon'})
    cellsize = cellsize.interp_like(SWE)
    
    SWE = SWE*cellsize * isglac #M to M3
    # with rio.open(join(files,'clones','cellsize05min.correct.map')) as cs:
    #     cs_=cs
    #     cellsize = cs.read(1)
        
    SWE_sum = SWE.sum(axis=(1,2))
    # SWE_sum.plot()
    # plt.figure()
    # SWE_sum.sel(time='2001').plot()
    # SWE_yearsum = []
    # for year in range(2001,2013):
    #     daterange = pd.date_range(str(year-1)+'-10-01',str(year)+'-09-30')
    #     SWE_yearsum.append(SWE_sum.sel(time=daterange).sum())
    #     # SWE_yearsum.append(pd.DataFrame({'Yearsum':SWE_sum.sel(time=daterange).sum().data},index=year))
    # SWE_yearsum = pd.DataFrame(SWE_yearsum,index=range(2001,2013))
    SWE_yearsum = []
    for year in range(2000,2013):
        # daterange = pd.date_range(str(year-1)+'-10-01',str(year)+'-09-30')
        if hemisphere=='North':
            SWE_yearsum.append(SWE_sum.sel(time=str(year)+'-09-30'))
        elif hemisphere=='South':
            SWE_yearsum.append(SWE_sum.sel(time=str(year)+'-03-31'))
        # SWE_yearsum.append(pd.DataFrame({'Yearsum':SWE_sum.sel(time=daterange).sum().data},index=year))
    SWE_yearsum = pd.DataFrame(SWE_yearsum,index=range(2000,2013))
    SWE_diff = SWE_yearsum.diff()[1:][0]
    SWE_diff_list.append(SWE_diff)
    
    
    
    GW = xr.open_dataset(join(daily_outputs,Basin_name+'_s_d_g_N0.nc')).groundwater_recharge
    GW2 = xr.open_dataset(join(daily_outputs,Basin_name+'_s_d_g_N2.nc')).groundwater_recharge
    GW_diff = GW - GW2 
    GW = GW.sel(time=slice('2000','2012'))
    
    
    GW = GW*cellsize * isglac #M to M3/day
    GW_sum = GW.sum(axis=(1,2))
    
    GW_wholebasin = GW*cellsize
    GW_wholebasin_sum = GW_wholebasin.sum(axis=(1,2))
    
    GW_diff = GW_diff*cellsize
    GW_diff_sum = GW_diff.sum(axis=(1,2))

    GW_yearsum = []
    GW_wholebasin_yearsum = []
    GW_diff_yearsum = []
    # for year in range(2000,2013):
    #     # daterange = pd.date_range(str(year-1)+'-10-01',str(year)+'-09-30')
    #     if hemisphere=='North':
    #         GW_yearsum.append(GW_sum.sel(time=str(year)+'-09-30'))
    #     elif hemisphere=='South':
    #         GW_yearsum.append(GW_sum.sel(time=str(year)+'-03-31'))
        # GW_yearsum.append(pd.DataFrame({'Yearsum':GW_sum.sel(time=daterange).sum().data},index=year))
    # GW_yearsum = pd.DataFrame(GW_yearsum,index=range(2000,2013))
    # GW_diff = GW_yearsum.diff()[1:][0]
    
    
    ### Discharge
    Q = xr.open_dataset(join(daily_outputs,Basin_name+'_s_d_g_N0.nc')).discharge
    Q = Q.sel(time=slice('2000','2012'))
    
    
    Q = Q #M to M3/day
    Q_sum = Q.sum(axis=(1,2))

    Q_yearsum = []
    for year in range(2000,2013):
        # daterange = pd.date_range(str(year-1)+'-10-01',str(year)+'-09-30')
        if hemisphere=='North':
            Q_yearsum.append(Q_sum.sel(time=str(year)+'-09-30'))
        elif hemisphere=='South':
            Q_yearsum.append(Q_sum.sel(time=str(year)+'-03-31'))
        # Q_yearsum.append(pd.DataFrame({'Yearsum':Q_sum.sel(time=daterange).sum().data},index=year))
    Q_yearsum = pd.DataFrame(Q_yearsum,index=range(2000,2013))
    # Q_diff = Q_yearsum.diff()[1:][0]
    
    
    
    
    
    # f1,ax1 =plt.subplots(figsize=(10,5))
    # GW_diff.plot(ax=ax1)
    # GW_diff_list.append(GW_diff)
    # SWE_diff = SWE_diff * GG_area0.sum()1
    # SWE_yearsum.plot()
    
    # SWE_yearsum = SWE_sum.groupby('time.year').sum()1
    # SWE_MB = SWE_yearsum[1:].data-SWE_yearsum[:-1]
    
    #
    hg = {}
    for setup in model_setups:
        hgdic = {'Model_name'          :Model_name,
           'Basin_name'          : Basin_name,
           'Setting_no'          :setup,
           'Fromyear'           :str(Fromyear),
           'Untilyear'          :str(Untilyear),
           'run_name'           :run_name,
           'station'            :'hg'}
        # if Basin_name == 'COLUMBIA':
        #     hgdic['station']='hg_1' 
        hg_name = '_'.join(hgdic.values()) +'.txt'
        hg['s'+setup]=pd.read_csv(
                join(run_dir,'Output','HG_2000_2012',
                hg_name)
                ,index_col=0,
                parse_dates=True).loc[daterange]
    hg_years = daterange.year.unique()[1:]
    
    hg_diff = (hg['s0']-hg['s2']) * 3600*24 
    hg_yearsum = []
    
    # NNSP = no spilling prevention 
    NNSP = pd.read_csv(
                join(run_dir,'Output','HG_2000_2012',
                'PCRG_'+Basin_name+'_2_2000_2012_NNSP_hg.txt')
                ,index_col=0,
                parse_dates=True).loc[daterange]
    spdiff = (hg['s2']-NNSP)*3600*24
    # f1,ax1 =plt.subplots(figsize=(10,5))
    # NNSP.plot(ax=ax1)
    # hg['s2'].plot(ax=ax1)
    
    # f1,ax1 =plt.subplots(figsize=(10,5))
    # spdiff.plot()
    spilling_dif=[]
    
    for year in range(2001,2013):
        # drange = pd.date_range(str(year-1)+'-10-01',str(year)+'-09-30')
        # hg_yearsum.append(hg_diff[drange].sum())
        if hemisphere=='North':
            rng = [str(year-1)+'-10-01',str(year)+'-09-30']
        elif hemisphere=='South':
            rng = [str(year-1)+'-04-01',str(year)+'-03-31']
        hg_yearsum.append(hg_diff[rng[0]:rng[1]].hg.sum())
        spilling_dif.append(spdiff[rng[0]:rng[1]].hg.sum())
        GW_yearsum.append(GW_sum.sel(time=slice(rng[0],rng[1])).sum())
        GW_wholebasin_yearsum.append(GW_wholebasin_sum.sel(time=slice(rng[0],rng[1])).sum())
        GW_diff_yearsum.append(GW_diff_sum.sel(time=slice(rng[0],rng[1])).sum())
        # Quantify percentage of spillingg prevention
        hg2_yearsum = hg['s2'][rng[0]:rng[1]].hg.sum()*3600*24
        spilling_yearsum = spdiff[rng[0]:rng[1]].hg.sum()
        print (Basin_name,spilling_yearsum/hg2_yearsum)
        
    f10  =plt.figure()
    plt.plot(spdiff/(hg['s2']*3600*24))
    plt.title('Nomalized spilling difference')
    f10.show()
    
    hg_yearsum = np.array(hg_yearsum)
    hg_diff_list.append(hg_yearsum)
    
    spilling = np.array(spilling_dif)
    spilling_list.append(spilling)
    # Spilling prevention
    GW_yearsum = pd.DataFrame(GW_yearsum,index=range(2001,2013))[0]
    GW_wholebasin_yearsum = pd.DataFrame(GW_wholebasin_yearsum,index=range(2001,2013))[0]
    GW_diff_yearsum = pd.DataFrame(GW_diff_yearsum,index=range(2001,2013))[0]
    plt.figure()
    # GW_wholebasin_yearsum.plot(title=Basin_name+' Total GW')
    GW_diff_yearsum.plot(title=Basin_name+ ' GW diff')
    GW_yearsum.plot()
    plt.show()

    ### Evap
    outmap = join(run_dir,r'Output\output_folders\output_PCRG_'+Basin_name+'_HadGEM2-ES_rcp26_0_2001_2016','netcdf')
    evap = xr.open_dataset(join(outmap,'totalEvaporation_annuaTot_output.nc'))
    evap= evap.interp_like(cellsize)*cellsize * isglac
    evap_sum = evap.sum(axis=(1,2)).total_evaporation
    f1,ax1 = plt.subplots()
    evap_sum.plot(ax=ax1)
    # evap_sum.coords['time']=evap.time
    #
    # f1,ax1 = plt.subplots()
    # SWE_diff.plot(ax=ax1,label='Benchmark snow towers')
    # MB_sum.plot(ax=ax1,label='GloGEM mass balance')
    # ax1.plot(SWE_diff.index,hg_yearsum,label='Coupled-benchmark')
    # ax1.plot(SWE_diff.index,MB_sum.values-SWE_diff.values,label='GGMB-PCRGMB')
    # plt.legend()
    
    
    
    from matplotlib.patches import Polygon
    
    f1,ax1 =plt.subplots(figsize=(10,5))
    posdif = -hg_yearsum
    ax1.plot(SWE_diff.index,posdif,label='Coupled model-Benchmark')
    
    # ix = SWE_diff.index
    # iy = SWE_diff
    # verts = [(SWE_diff.index[0], 0), *zip(ix, iy), (SWE_diff.index[-1], 0)]
    # poly = Polygon(verts, facecolor='tab:green', alpha=0.5)
    # ax1.add_patch(poly)
    
    # ix = MB_sum.index
    # iy = -MB_sum+SWE_diff
    # verts = [(MB_sum.index[0], SWE_diff[2001]), *zip(ix, iy), (MB_sum.index[-1], SWE_diff[2012])]
    # poly = Polygon(verts, facecolor='tab:red', alpha=0.5)
    # ax1.add_patch(poly)
    # ax1.stackplot(SWE_diff.index,SWE_diff,-MB_sum,spilling,
    #               colors=['tab:brown','tab:red','tab:green'],alpha=0.5,
    #               labels=['Snow towers','GloGEM mass loss','Spilling'])
  
    ax1.stackplot(SWE_diff.index,SWE_diff,-MB_sum,spilling,GW_yearsum,
                  colors=['tab:brown','tab:red','tab:green','tab:olive'],alpha=0.5,
                  labels=['Snow towers','GloGEM mass loss','Spilling','Groundwater recharge'])
    # ax1.plot(GW_yearsum,color='red')
    ax1.legend()
    ax1.grid()
    ax1.set_ylabel(r'$Mass imbalance\ (m^{3})$')
    ax1.set_title(Basin_name.title())
    ax1.set_xlim(SWE_diff.index[0],SWE_diff.index[-1])
    f1.show()
    
    print ()
    
    
    
f1,ax1 = plt.subplots(figsize=(40,8))
hg['s2'].plot(ax=ax1)
NNSP.plot(ax=ax1,label='NNSP')
ax1.set_title('NNSP vs normal')
ax1.legend()

NNSPnc = xr.open_dataset('d:\Documents\Master vakken\Thesis\Code\Files\glaciers_nc\COLUMBIA_2000_2012_NNSP_R.nc')
ncR = NNSPnc.R.mean(dim=('lat','lon'))
ncR.plot(figsize=(30,10))
#%% 
import matplotlib.ticker
from matplotlib import cm

class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = r'$\mathdefault{%s}$' % self.format


f1,axes = plt.subplots(4,1,figsize=(10,8),sharex=True)
plt.subplots_adjust(hspace=0.24)
# colors = ['grey','tab:blue','tab:red']
turbo = cm.get_cmap('turbo',10)
colors = turbo((9,1,5))
        
for i in range(len(Basin_names)):
    ax = axes[i]
    posdif = -hg_diff_list[i]
    ax.plot(SWE_diff_list[i].index,posdif,label='Coupled model - Benchmark',color='black')
    ax.stackplot(SWE_diff_list[i].index,SWE_diff_list[i],-MB_sum_list[i],spilling_list[i],
                 colors=colors,alpha=0.7,
                  labels=['PCR-GLOBWB 2 snow towers','GloGEM net mass loss','Spilling prevention'])
    ax.grid(alpha=0.6)
    if i==0:
        ax.legend(loc = (0.68,1.02))
        # handles, labels = ax.get_legend_handles_labels()
        # ax.legend(handles[::-1], labels[::-1], loc=(0.68,1.02))
    ax.set_title(Basin_names[i].title())
    # ax.set_ylabel(r'$Mass\ imbalance\ (m^{3})$')
    ax.set_xlim(SWE_diff.index[0],SWE_diff.index[-1])
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_formatter(OOMFormatter(9,"%1.i"))
f1.text(0.08, 0.5, r'$Mass\ imbalance\ [m^{3}]$', va='center', rotation='vertical',
        size='large')
# f1.savefig(join(figures,'MB_diffs_all_colorblind.svg'),format='svg',bbox_inches='tight')

