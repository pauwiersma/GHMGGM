# -*- coding: utf-8 -*-
"""
@author: Pau Wiersma

This script
    Calculates the annual runoff difference between the coupled model and the benchmark
    and the contribution of
        Snow towers
        Glacier mass loss 
        Spilling prevention
    to this annual runoff difference
    for the basins Alsek, Columbia, Oelfusa and Rhone

Files needed
    -HH2018 annual mass balance
    -Subbasin shapefile 
    -HH2018 Overview.dat
    -daily SWE raster
        -take annual mass balance from glacfrac
        - Compare timing of snowmelt bewteen N0 and N2
    - Hydrographs
        - Take difference between hg0 and hg2 and see what percentage the mass balances take in it
    -Cellsize raster
    -Glacier fraction (isglac) from nc file
    -GHMGGM_basin_info.csv

"""

import warnings
warnings.filterwarnings('ignore')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from os.path import join
import geopandas as gp
import xarray as xr

import rioxarray as xrio
import matplotlib.ticker
from matplotlib import cm

run_dir   = r'D:\Documents\Master Vakken\Thesis\Code'
os.chdir(run_dir)
files               = join(run_dir,'Code','Files')

#%% INPUTS
files               = join(run_dir,'Files')
figures             =join(run_dir,'Figures')
daily_outputs       =join(run_dir,'Output','daily_outputs')
HH_MB               = join(files,'HH2018_validation')

basin_info = pd.read_csv(join(files,'GHMGGM_basin_info.csv'),index_col = 0)

Basin_names = ['ALSEK','COLUMBIA','OELFUSA','RHONE']

model_setups        = ['0','1','2']

save_figs           = False


Model_name          ='PCRG'
GG_GCM              ='HadGEM2-ES'
GG_rcp              ='rcp26'
Fromyear            =2000
Untilyear           =2012

#%% Function for ytick formatting 
#(from https://stackoverflow.com/questions/42656139/set-scientific-notation-with-fixed-exponent-and-significant-digits-for-multiple)
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
#%%

run_name = 'N'
#Prepare lists for mass balance components per basin
SWE_diff_list = []
MB_sum_list = []
hg_diff_list = []
spilling_list = []

#Loop through each basin and assemble the mass balance components
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
    
    ###GloGEM mass loss
    # Load the GloGEM point data
    shp =   gp.read_file(join(files,'subbasin_geojsons',Basin_name+'.geojson'))
    Overview_path = basin_info.loc[Basin_name,'RGI_dir']+'_HadGEM2-ES_RCP26_Overview.dat'
    Overview    = pd.read_csv(join(HH_MB,Overview_path),
                            delim_whitespace=True,index_col=0)
    points              = gp.GeoSeries(gp.points_from_xy(Overview['lon'].values,Overview['lat'])
                                        ,index = Overview.index)

    GG_idx     = [shp.geometry.contains(point)[0] for point in points]
    GG_points   = points[GG_idx]
    GG_area0    = Overview['Area0'][GG_idx] # in km2

    #Load annual mass balance data for each individual glacier
    MB_path  = basin_info.loc[Basin_name,'RGI_dir']+'_HadGEM2-ES_RCP26_Annual_balance.dat'
    MB = pd.read_csv(join(HH_MB,MB_path),delim_whitespace=True,index_col=0)
    MB =MB[[str(year) for year in np.arange(2001,2013)]]
    MB = MB[GG_idx]
    MB = MB.multiply(GG_area0,axis=0) * 1e6 #meter to m*km2 to m*m2= m3
    MB[MB==-99] = np.nan
    #Calculate the total basin glacier mass balance
    MB_sum =MB.sum(axis=0)
    MB_sum.index = np.arange(2001,2013)
    MB_sum_list.append(MB_sum)

    ###Snow towers
    #Load SWE output of the benchmark run and calculate how much survives each summer
    SWE = xr.open_dataset(join(daily_outputs,Basin_name+'_s_d_g_N0.nc')).snow_water_equivalent
    SWE = SWE.sel(time=slice('2000','2012'))
    # SWE_aletsch = SWE.isel(lat =slice(41,43),lon=60).sum(axis=1)
    nc_path = join(files,'glaciers_nc',Basin_name+'_2000_2012_N_R.nc')
    isglac = xr.open_dataset(nc_path).isglac
    cellsize =xrio.open_rasterio(join(files,'clones','cellsize05min.correct.map')).squeeze()
    cellsize = cellsize.rename({'y':'lat','x':'lon'})
    cellsize = cellsize.interp_like(SWE)
    #M to M3 and correct for glacizeration fraction
    SWE = SWE*cellsize * isglac 
    SWE_sum = SWE.sum(axis=(1,2))
    SWE_yearsum = []
    #Store the total SWE on the glacierized area at the end of the hydrological year 
    for year in range(2000,2013):
        if hemisphere=='North':
            SWE_yearsum.append(SWE_sum.sel(time=str(year)+'-09-30'))
        elif hemisphere=='South':
            SWE_yearsum.append(SWE_sum.sel(time=str(year)+'-03-31'))
        # SWE_yearsum.append(pd.DataFrame({'Yearsum':SWE_sum.sel(time=daterange).sum().data},index=year))
    SWE_yearsum = pd.DataFrame(SWE_yearsum,index=range(2000,2013))
    SWE_diff = SWE_yearsum.diff()[1:][0]
    SWE_diff_list.append(SWE_diff)

    ### Annual runoff difference
    #Load coupled and benchmark hydrographs and calculate annual difference
    hg = {}
    for setup in model_setups:
        hgdic = {'Model_name'          :Model_name,
           'Basin_name'          : Basin_name,
           'Setting_no'          :setup,
           'Fromyear'           :str(Fromyear),
           'Untilyear'          :str(Untilyear),
           'run_name'           :run_name,
           'station'            :'hg'}
        hg_name = '_'.join(hgdic.values()) +'.txt'
        hg['s'+setup]=pd.read_csv(
                join(run_dir,'Output','HG_2000_2012',
                hg_name)
                ,index_col=0,
                parse_dates=True).loc[daterange]
    hg_years = daterange.year.unique()[1:]

    hg_diff = (hg['s0']-hg['s2']) * 3600*24 #m3/s to m3/day
    hg_yearsum = []

    ### Spilling prevention
    # NNSP = no spilling prevention
    # Load coupled model hydrograph without spilling prevention and calculate difference
    NNSP = pd.read_csv(
                join(run_dir,'Output','HG_2000_2012',
                'PCRG_'+Basin_name+'_2_2000_2012_NNSP_hg.txt')
                ,index_col=0,
                parse_dates=True).loc[daterange]
    spdiff = (hg['s2']-NNSP)*3600*24 #m3
    spilling_dif=[]

    # Store spilling difference and annual runoff difference per year 
    for year in range(2001,2013):
        if hemisphere=='North':
            rng = [str(year-1)+'-10-01',str(year)+'-09-30']
        elif hemisphere=='South':
            rng = [str(year-1)+'-04-01',str(year)+'-03-31']
        hg_yearsum.append(hg_diff[rng[0]:rng[1]].hg.sum())
        spilling_dif.append(spdiff[rng[0]:rng[1]].hg.sum())
        # Quantify percentage of spillingg prevention
        hg2_yearsum = hg['s2'][rng[0]:rng[1]].hg.sum()*3600*24
        spilling_yearsum = spdiff[rng[0]:rng[1]].hg.sum()
        print (Basin_name,year,'(spilling prevention)/(total runoff difference) =',spilling_yearsum/hg2_yearsum)

    hg_yearsum = np.array(hg_yearsum)
    hg_diff_list.append(hg_yearsum)

    spilling = np.array(spilling_dif)
    spilling_list.append(spilling)
    # Spilling prevention
    plt.figure()

    f1,ax1 =plt.subplots(figsize=(10,5))
    posdif = -hg_yearsum
    ax1.plot(SWE_diff.index,posdif,label='Coupled model-Benchmark')

    ax1.stackplot(SWE_diff.index,SWE_diff,-MB_sum,spilling,
                  colors=['tab:brown','tab:red','tab:green','tab:olive'],alpha=0.5,
                  labels=['Snow towers','GloGEM mass loss','Spilling','Groundwater recharge'])
    # ax1.plot(GW_yearsum,color='red')
    ax1.legend()
    ax1.grid()
    ax1.set_ylabel(r'$Mass imbalance\ (m^{3})$')
    ax1.set_title(Basin_name.title())
    ax1.set_xlim(SWE_diff.index[0],SWE_diff.index[-1])
    f1.show()


#%%

f1,axes = plt.subplots(4,1,figsize=(10,8),sharex=True)
plt.subplots_adjust(hspace=0.24)
# colors = ['grey','tab:blue','tab:red']
turbo = cm.get_cmap('turbo',10)
colors = turbo((9,1,5))

for i,B in enumerate(Basin_names):
    ax = axes[i]
    posdif = -hg_diff_list[i]
    ax.plot(SWE_diff_list[i].index,posdif,label='Coupled model - Benchmark',color='black')
    ax.stackplot(SWE_diff_list[i].index,SWE_diff_list[i],-MB_sum_list[i],spilling_list[i],
                 colors=colors,alpha=0.7,
                  labels=['PCR-GLOBWB 2 snow towers','GloGEM net mass loss','Spilling prevention'])
    ax.grid(alpha=0.6)
    if i==0:
        ax.legend(loc = (0.68,1.02))
    ax.set_title(B.title())
    ax.set_xlim(SWE_diff.index[0],SWE_diff.index[-1])
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_formatter(OOMFormatter(9,"%1.i"))
f1.text(0.08, 0.5, r'$Runoff\ difference\ [m^{3}/year]$', va='center', rotation='vertical',
        size='large')
f1.savefig(join(figures,'MB_diffs_all_colorblind_runoffdifference_m3year.svg'),format='svg',bbox_inches='tight')

