# -*- coding: utf-8 -*-

#%%
"""
Created on Mon Jan 25 16:27:45 2021

@author: Pau Wiersma

Scripts accompanying PCR-GLOBWB 2 and GloGEM coupling
(1/3): Preprocessing



This script:
    - helps to select GRDC stations and to create subbasins and clonemaps
            based on their location
    - rasterizes and resamples GloGEM data
    - prevents their spilling into other basins
    - saves the result in a NetCDF file

Files needed:
    basin_info_45min.csv
    Huss & Hock (2018) data per basin:
        Overview.dat
    Hydrosheds data per global region
        hybas (lvl 00)
        bas
    GRDC metadata
    GRDC runoff observations

    
To do:
    isbasins
    List all necessary files
    Sync rasterize functions
"""
#%% Packages
import datetime as dt
import glob
import json
import os
import warnings
from functools import partial
from os.path import join

import shutil
import subprocess
import sys
import time

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from geocube.api.core import make_geocube
from geocube.rasterize import rasterize_image
from rasterio.enums import MergeAlg
from shapely.geometry import box, mapping

warnings.filterwarnings('ignore')


from matplotlib.lines import Line2D





#%% Options and settings

#Set running directory
run_dir = r'D:\Documents\Master Vakken\Thesis\Code'
os.chdir(run_dir)
files               = join(run_dir,'Files')


GG_GCM              = 'HadGEM2-ES'  #GloGEM GCM (irrelevant for historical sims)
GG_RCP              = 'rcp26'       #GloGEM rcp scenario (irrelevant for historical sims)
GG_FROMYEAR         = 2000  #First hydrological year (so start from October/April of GG_FROMYEAR-1)
GG_UNTILYEAR        = 2012 #Last hydrological year
GHM_RESOLUTION      = 0.0833333



#Other inputs
FIND_BASINS         =False #True=create basin shapefiles, false=load
FIND_SINK           =False #True=find, false=load from basin_info.csv

# If FIND_GRDC_STATIONS==True, then the subbasins will not be searched and the run will 
# stop after the plots. The plot will contain the GRDC stations that have observations 
# within OBS_MINYEAR and OBS_MAXYEAR with a minimum length of OBS_MINLEN years.
# Select the station, download the data from the website, add it to basin_info.csv and
#set FIND_GRDC_STATIONS to False again
FIND_GRDC_STATIONS   =False
OBS_MINYEAR         =2000
OBS_MAXYEAR         =2012
OBS_MINLEN          =5

SAVE_FIG = False

#Find_subbasin==True will find the subbasin upstream of the selected GRDC station
FIND_SUBBASIN       =False #True=find based on GRDC_station, False=Load
FIND_CLONE          =False #True= find clone with margins that are compatible with PCRGLOB, False = Load
PLOT_BASIN          =False



FIND_ISOUT_ISGLAC   =False #True=Create, False = Load
FIND_QM             =True #True=Create, False = Load
FIND_QD             =True #True=Create, False = Load
ERA_RESAMPLE        =True #True= T-based weight resampling, False= Equal weight resampling
SPILLING_PREVENTION =True #False=skip this step
INCLUDE_SP_AUX      =True #True=include spilling prevention aux data to nc-file: isout,isglac,glacfrac
INCLUDE_HG          =True #True= include GRDC observations in nc-file
SAVE_NC             =False #Set to False for testruns
CREATE_WORLDMAP     =False
RUN_NAME            ='Na0' #run ID for savefiles

#Variables for temperature weighted resampling of GloGEM data
RESAMPLING_METHOD   = 'Backfill' #or Linear (doesn't work)
WEIGHT_FACTOR       = 20 #the higher the flashier, 30 is probably best
T_THRESHOLD         =268.15 #Daily mean T threshold above which weight is given to a day in the ERA resampling


#Load basin_info.csv with all basin info
basin_info = pd.read_csv(join(files,'basin_info_45min.csv'),index_col = 0)
basin_names = basin_info[basin_info['suitable']=='y'].index
basin_names = ['RHONE','ALSEK']


#%% Create dictionary for each basin in loop
Basins={}
for B in basin_names:
    Basins[B]={}
    for k in basin_info.keys():
        Basins[B][k]=basin_info.loc[B,k]

#%% Load paths, both for loading and for saving
for B in basin_names:
    
    RGI_dir          = Basins[B]['RGI_dir']
        
    p     ={'basin_shp'     :join(files,'basin_geojsons',B+'.geojson'),
            'subbasin_shp'  :join(files,'subbasin_geojsons',B+'.geojson'),
            'isout'         :join(files,'isout','isout_'+B+'.nc'),
            'isglac'        :join(files,'isglac','isglac_'+B+'.nc'),
            'isbasin'       :join(files,'isbasin','isbasin_'+B+'.nc'),
            'clonemaps'     :join(files,'clones'),
            'globglacfrac'  :join(files,'globglacfrac','glob_glac_frac.shp'),
            'glob_mask'     :join(files,'globglacfrac','global_masked_vectorgrid.shp'),
            'Qcatch'        :join(files,'HH2018_validation',RGI_dir+'_'+GG_GCM+'_'+GG_RCP.upper()+'_Discharge_catchment.dat'),
            'Qm'            :join(files,'GG_rasterized_monthly',B+'_Qm_'+str(GG_FROMYEAR)+'_'+str(GG_UNTILYEAR)+'.nc'),
            'Qd'            :join(files,'GG_resampled_daily',B+'_Qd_'+str(GG_FROMYEAR)+'_'+str(GG_UNTILYEAR)+'.nc'),
            'Qsp'           :join(files,'GG_spilling_prevented',B+'_Qsp_'+str(GG_FROMYEAR)+'_'+str(GG_UNTILYEAR)+'.nc'),
            'Qfinal'        :join(files,'GG_final',B+'_Qfinal_'+str(GG_FROMYEAR)+'_'+str(GG_UNTILYEAR)+'.nc'),
            'Overview'      :join(files,'HH2018_validation',RGI_dir+'_'+GG_GCM+'_'+GG_RCP.upper()+'_Overview.dat'),
            'ERA_tas'       :join(files,'ERA_Interim','pcrglobwb_OBS6_ERA-Interim_reanaly_1_day_tas_1999-2016_'+B+'.nc'),
            'ERA_pr'        :join(files,'ERA_Interim','pcrglobwb_OBS6_ERA-Interim_reanaly_1_day_pr_1999-2016_'+B+'.nc'),
            'nc_out'        :join(files,'glaciers_nc'),
            'Hydrosheds'    :join(files,'HydroSHEDS'),
            'yan_rivers'    :join(files,r'Yan2019rivers')}
    
    #Load GRDC observation paths to include in NC-file
    if B=='RHONE':
        p['Beaucaire']  = join(files,'Q_obs','Beaucaire_1980_2014.csv')
        p['Tarrascon']  =join(files,'Q_obs','Tarrascon_2010_2020.csv')
    elif B =='INDUS':
        p['Tarbela']  = join(files,'Q_obs','Tarbela2015_2016.csv')
    else:
        GRDC_no         = str(int(Basins[B]['grdc_no1']))
        p['GRDC']=glob.glob(join(files,'Q_obs',GRDC_no+'*.txt'))
    Basins[B]['p']=p

#%% Find main basins with Hydrosheds basins
if FIND_BASINS:
    for B in basin_names:
        bas_cont        = Basins[B]['continent']
        bas_cont_dir    = bas_cont+'_bas_dissolved'
        bas_dir         = join(p['Hydrosheds'],'bas_files',bas_cont_dir,bas_cont_dir+'.shp')
        Basins[B]['shp'] = gp.read_file(bas_dir,
                bbox=(Basins[B]['center_lon'],Basins[B]['center_lat'],
                      Basins[B]['center_lon'],Basins[B]['center_lat']))
        Basins[B]['shp'].to_file(Basins[B]['p']['basin_shp'],driver='GeoJSON')
else:
    for B in basin_names:
        Basins[B]['shp'] = gp.read_file(Basins[B]['p']['basin_shp'])

#%% Load overview.dat with glacier points and areas and store in dictionary
for B in basin_names:
    print (B + ' loading glacier points')
    Overview    = pd.read_csv(Basins[B]['p']['Overview'],
                            delim_whitespace=True,index_col=0)
    points              = gp.GeoSeries(gp.points_from_xy(Overview['lon'].values,Overview['lat'])
                                        ,index = Overview.index)
    
    Basins[B]['GG_idxs']     =[Basins[B]['shp'].geometry.contains(point)[0] for point in points]
    Basins[B]['GG_points']   = points[Basins[B]['GG_idxs']]
    Basins[B]['GG_area0']    = Overview['Area0'][Basins[B]['GG_idxs']] # in km2
    #Testplot:
        #Basins['RHONE']['GG_points'].plot(markersize=Basins['RHONE']['GG_area0'])

#%% Load Hybas
def load_hybas(Basins_,B_):
    bas_cont        = Basins_[B_]['continent']

    hybas_cont_dir      = 'hybas_'+bas_cont+'_lev00_v1c'
    hybas_dir           =os.path.join(Basins_[B_]['p']['Hydrosheds'],
                                      'hybas_files',hybas_cont_dir,hybas_cont_dir+'.shp')
    hybas_full               = gp.read_file(hybas_dir) #Full hybas for continent
    
    main_bas   = hybas_full['MAIN_BAS'][hybas_full.geometry.contains(
        gp.points_from_xy([Basins_[B_]['center_lon']],[Basins_[B_]['center_lat']])[0])].values[0]
    return hybas_full[hybas_full['MAIN_BAS']==main_bas]


#%% Find/Load glacier sinks (most downstream point to include all glacier runoff)
t = time.time()

# Glacier sink and GRDC location

if FIND_SINK:
    for B in basin_names:
        print ('Searching for glacier sink in ', B)
        hybas = load_hybas(Basins,B) #Hydrosheds subcatchments
        hybas_glac      = [] #Create list with subcatchments containing RGI glaciers
        for  point in Basins[B]['GG_points'].values:
            hybas_glac.append(hybas[hybas.geometry.contains(point)])
        hybas_glac = pd.concat(hybas_glac)
        
        #Determine routing using the Hydrosheds Pfafstetter coding system 
        #https://en.wikipedia.org/wiki/Pfafstetter_Coding_System & https://www.hydrosheds.org/page/hydrobasins
        for i in range(1,13):
            key = 'PFAF_'+str(i)
            #look for first PFAF that's not the same for all
            unq = np.unique(hybas_glac[key])
            if len(unq)>1:
                un = np.min(unq)
                if un%10!=1: # If the lowest PFAF hybas the sink of the subbasin
                    levelx = hybas[hybas[key]==un]
                    next_down = levelx.sort_values(by='DIST_MAIN').head(1).NEXT_DOWN
                    glac_sink=hybas[hybas.HYBAS_ID.values==next_down.values[0]].centroid
                else: 
                    level1 = hybas[hybas[key]==un]
                    glac_sink = level1.sort_values(by='DIST_MAIN').head(1).centroid
                break
        Basins[B]['glac_sink'] = glac_sink
        basin_info.glac_sink_lon[B] = glac_sink.x.values[0]
        basin_info.glac_sink_lat[B] = glac_sink.y.values[0]
    basin_info.to_csv(join(files,'basin_info_45min.csv'))
else:
    for B in basin_names:
        x = [Basins[B]['glac_sink_lon']]
        y = [Basins[B]['glac_sink_lat']]
        Basins[B]['glac_sink']=gp.GeoSeries(gp.points_from_xy(x,y))

#%% Find GRDC stations
if FIND_GRDC_STATIONS:
    #Load GRDC station metadata
    GRDC_meta       = pd.read_csv(os.path.join(files,'GRDC_Stations.csv')
                                  ,na_values = 'n.a.', index_col = 0)
    
    for B in basin_names:
        GRDC_name = Basins[B]['GRDC_name']
        
        GR_stats = GRDC_meta[(GRDC_meta.river==GRDC_name)]
        if B=='NEGRO':
            GR_stats = GR_stats[GR_stats.country=='AR'] #There are multiple
        
        if GR_stats.size==0: 
            print(GRDC_name+' doens"t appear in GRDC metadata')
            continue
        
        #Find all stations that fulfill the OBS requirements
        GR_day = GR_stats[(GR_stats.d_start<=OBS_MAXYEAR)&(GR_stats.d_end>=OBS_MINYEAR)]
        GR_day = GR_day[GR_day.d_end-GR_day.d_start>=OBS_MINLEN]
        if GR_day.size==0:
            print ('No GRDC station with full daterange for '+GRDC_name)
            continue
        
        print (B+' GRDC: '+str(GR_day.shape[0])+' stations loaded with day data')
        Basins[B]['GRDC_stationlist']=GR_day
    ##### Add GRDC station info to basin_info45min.csv

#%% Find subbasin
if not FIND_GRDC_STATIONS:
    if FIND_SUBBASIN:
        for B in basin_names:
            print (B+' Creating subbasin')
            
            if B in ['OELFUSA','THJORSA']: #Station is just downstream of basin borders
                print ('Station falls outside of basin, subbasin equals basin')
                Basins[B]['subshp'] = Basins[B]['shp']
                Basins[B]['subshp'].to_file(Basins[B]['p']['subbasin_shp'],driver='GeoJSON')
                continue
            
            hybas = load_hybas(Basins,B)
            
            station_point = gp.GeoSeries(gp.points_from_xy([Basins[B]['station_lon1']],
                                                          [Basins[B]['station_lat1']]))[0]
            hybas_station      = hybas[hybas.geometry.contains(station_point)]
            
            if B in ['COPPER','OB']:
                hybas_station=hybas[hybas.HYBAS_ID.values==hybas_station.NEXT_DOWN.values]
        
            #Select all stations with larger distance to the river sink
            larger_dist     = hybas[hybas.DIST_MAIN>hybas_station.DIST_MAIN.values[0]]
            t = time.time()
            for i in range(1,13):
                key = 'PFAF_'+str(i) 
                statpfaf = hybas_station[key].values[0]
                if not np.all(larger_dist[key]==statpfaf):
                    key = 'PFAF_'+str(i-1) #The last PFAF that all have in common
                    statpfaf = hybas_station[key].values[0] 
                    # print (key)
                    selection = larger_dist[larger_dist[key]==statpfaf]
                    # Now we have a selection of hybas that share the highest common PFAF, 
                    # but there may be hybas that don't drain in hybas_station
                    # A first check is to see if there are hybas with higher UP_AREA, that would mean
                    # hybas_station is not in the main branch of the river
                    if np.any(selection.UP_AREA>hybas_station.UP_AREA.values[0]):
                        #If station_hybas is not in the main branch, then the whole subbasin will be inside
                        # this PFAF
                        keyplusone = 'PFAF_'+str(i+1)
                        statpfafplusone=hybas_station[keyplusone].values[0]
                        selection = selection[selection[keyplusone]==statpfafplusone]
                        
                    # Create a blacklist, loop through hybas and add any basins that don't eventually drain into hybas_station
                    blacklist = np.zeros(len(selection),dtype=bool)
                    for j in range(len(selection)):
                        #Start from hybas and follow the river to see if it passes hybas_station
                        next_down = selection.iloc[j,:]
                        found = False
                        while not found:
                            if next_down.NEXT_DOWN==hybas_station.HYBAS_ID.values[0]:
                                found=True
                            else:
                                next_down = larger_dist[larger_dist.HYBAS_ID==next_down.NEXT_DOWN]
                                if len(next_down)==0:
                                    # print (selection.iloc[j,:].HYBAS_ID,' does not drain in hybas_station')
                                    blacklist[j]=1
                                    break
                                next_down = next_down.iloc[0,:]
                    break
            sub_basin = selection.iloc[~blacklist,:].unary_union
            Basins[B]['subshp']=gp.GeoSeries(sub_basin)
            Basins[B]['subshp'].to_file(Basins[B]['p']['subbasin_shp'],driver='GeoJSON')
    else:
        for B in basin_names:
            Basins[B]['subshp'] = gp.read_file(Basins[B]['p']['subbasin_shp'])
                              

#%% Find extents
if not FIND_GRDC_STATIONS:
    for B in basin_names:
        # Add a margin of 0.5 and then find the extents to match the 45 min resolution of ERA-Interim
        if B=='COPPER': #exception
            LLC=[-149,60]
            URC=[-137,66]
        else:
            LLC       = np.int32(np.floor(Basins[B]['subshp'].total_bounds[:2]-0.5))
            URC       = np.int32(np.ceil(Basins[B]['subshp'].total_bounds[2:]+0.5))
            #If the extents are not divisable by 3 they're not suited for 45 arcmin forcing
            while abs(URC[0]-LLC[0])%3!=0:
                URC[0]+=1
            while abs(URC[1]-LLC[1])%3!=0:
                URC[1]+=1
        
        Basins[B]['llc_lon']    =LLC[0]
        Basins[B]['llc_lat']    =LLC[1]
        Basins[B]['urc_lon']    =URC[0]
        Basins[B]['urc_lat']    =URC[1]
        
        basin_info.loc[B,'llc_lon']  =LLC[0]
        basin_info.loc[B,'llc_lat']    =LLC[1]
        basin_info.loc[B,'urc_lon']    =URC[0]
        basin_info.loc[B,'urc_lat']    =URC[1]
        
    basin_info.to_csv(join(files,'basin_info_45min.csv'))

#%% Make clonemaps

# This section contains a workaround due to spaces in the directory, 
# Rewrite if your directory doesn't contain spaces
if FIND_CLONE:
    for B in basin_names:
        print ('Clone '+B)
        source = join(r'd:\Documents\clonemaps','clone_global_05min.map')
        temp_dest = join(r'd:\Documents\clonemaps',B+'_05min.clone.map')
        dest   = join(Basins[B]['p']['clonemaps'],B+'_05min.clone.map')
        
        def coord2string(coord):
            return Basins[B][coord].astype(float).astype(str)
        coord2string('llc_lon')
        
        command     = 'gdal_translate -of PCRASTER -projwin '+\
            coord2string('llc_lon')+' '+coord2string('urc_lat')+\
                ' '+coord2string('urc_lon')+' '+coord2string('llc_lat')+' '\
            +source+' '\
            + temp_dest
        subprocess.Popen(command,shell=True)
        time.sleep(5)
        shutil.move(temp_dest,dest)
#%% Plot
if PLOT_BASIN:
    for B in basin_names:
        print ('Plot '+B)
        FIGSIZE=10
        SHAPEX=0.8
        f1 = plt.figure(figsize=(FIGSIZE,SHAPEX*FIGSIZE))
        ax1 = f1.add_subplot(projection=ccrs.PlateCarree())
        ax1.set_extent((Basins[B]['llc_lon'],Basins[B]['urc_lon'],
                        Basins[B]['llc_lat'],Basins[B]['urc_lat']))
        
        #Plot basemap
        ax1.coastlines()
        ax1.gridlines(draw_labels=True)
        ax1.add_feature(cfeature.LAND,color='green',alpha=0.3)
        ax1.add_feature(cfeature.OCEAN)
        ax1.add_feature(cfeature.BORDERS,linestyle=':',linewidth=1.5)
        # ax1.add_feature(cfeature.RIVERS,color='b',linewidth=1)
        ax1.add_feature(cfeature.LAKES,color='b')
        
        cont = Basins[B]['continent']
        if cont=='ar':
            cont='na'
        if cont =='si':
            cont='as'
        
        #Plot rivers
        #Rivers of Yan et al. 2019 are higher resolution than the cartopy set (https://doi.org/10.1038/s41597-019-0243-y)
        yan_river_path = join(p['yan_rivers'],cont+'_river.shp')
        yan_river = gp.read_file(yan_river_path,bbox=
                                 (Basins[B]['llc_lon'],Basins[B]['llc_lat'],
                                  Basins[B]['urc_lon'],Basins[B]['urc_lat']))
        yan_river.plot(ax=ax1,color='b',alpha=0.7)
        Basins[B]['shp'].plot(ax=ax1,alpha=0.3,label=B,color='tab:blue')
        
        #Test with hybas:
        # hybas =load_hybas(Basins,B)
        # hybas.plot(ax=ax1,facecolor=None,edgecolor='black',alpha=0.5)
        
        #Station plot
        if not FIND_GRDC_STATIONS: #Plot the selected GRDC station
            Basins[B]['subshp'].plot(ax=ax1,color='tab:blue',alpha=0.5,edgecolor='black')
            GRDC = ax1.plot(Basins[B]['station_lon1'],Basins[B]['station_lat1'],
                            'v',color='yellow',label='Gauging station',
                            markeredgecolor='black')[0]
        else:         #Plot all stations found to match the date requirements
            stations=Basins[B]['GRDC_stationlist']
            stations.plot(x='long',y='lat',kind='scatter',
                          ax=ax1,color='y',s=30)
            xtext = 0.1*(abs(ax1.get_xlim()[1]-ax1.get_xlim()[0]))

            for idx in stations.index:
                xy = np.array([stations.loc[idx,'long'],
                                  stations.loc[idx,'lat']])
                ax1.annotate(idx,xy,xy-np.array([xtext,0]),
                             bbox=dict( fc="0.9",alpha=0.5))
            GRDC =Line2D([0],[0],marker='o',linestyle='None',
                   markerfacecolor='y',markeredgecolor='y',
                   markersize=5,label='GRDC stations') 
        
        #Plot all RGI glacier points in the basin, with the point size scaled based on the glacier area
        Basins[B]['GG_points'].plot(ax=ax1,color='red',
                        markersize=Basins[B]['GG_area0']*2,label='Glaciers')
        
        #Plot glacier sink
        if len(Basins[B]['glac_sink'])==0:
            print ('No glacier sink found')
        else:
            Basins[B]['glac_sink'].plot(ax=ax1,marker='^',
                            color='lime',label='Glacier sink',
                            zorder = 3,edgecolor='black')
            green_dot =  Line2D([0],[0],marker='^',linestyle='None',
                   markerfacecolor='lime',markeredgecolor='black',
                   markersize=6,label='Glacier sink')
        ax1.set_title(B.title())
        red_dot = Line2D([0],[0],marker='o',linestyle='None',
                   markerfacecolor='red',markeredgecolor='red',
                   markersize=5,label='Glaciers')
        ax1.legend(handles = [red_dot,GRDC,green_dot])
        # if SAVE_FIG == True:
        #     save_at = join(run_dir,r'Code\Figures\Basin_maps\subbasins',B+'_map')
        #     plt.savefig(save_at,bbox_inches='tight')
        #     plt.show()
        
        #save as vectorplot 
        if SAVE_FIG:
            save_at = join(run_dir,r'Code\Figures\Basin_maps\subbasins',B+'_map.svg')
            plt.savefig(save_at,format = 'svg',bbox_inches='tight')
            plt.show()
            
if FIND_GRDC_STATIONS:
    print ('Exit to find GRDC stations')
    sys.exit()

#%% Rasterize function using geocube
def rasterize(vector,var_key,merge_algorithm='replace'):
    vector['bool'] = 1
    llc_lon,llc_lat,urc_lon,urc_lat = [
        Basins[B][i] for i in ['llc_lon','llc_lat','urc_lon','urc_lat']]
    lon360 = 360 * (llc_lon<0)
    latflip = 1-2*(llc_lat<0)
    vector.geometry = vector.geometry.affine_transform([1,0,0,latflip,lon360,0])
    res = 0.0833333
    
    #rasterio.enums.MergeAlg algorithms, default is replace
    if merge_algorithm == 'add':
        MERGE_ALG=MergeAlg.add
    elif merge_algorithm =='replace':
        MERGE_ALG=MergeAlg.replace
    geo_grid=make_geocube(vector_data=vector,
                           measurements = [var_key],
                          resolution = (res,res),
                          geom = json.dumps(mapping(box(llc_lon+lon360,
                                                        llc_lat*latflip,
                                                        urc_lon+lon360,
                                                        urc_lat*latflip))),
        rasterize_function=partial(rasterize_image,merge_alg = MERGE_ALG),
        fill = 0)
    geo_grid['x'] = geo_grid['x']-lon360
    geo_grid['y'] = geo_grid['y']*latflip
    geo_grid= geo_grid.rename({'y':'lat','x':'lon'})
    geo_grid = geo_grid.where(~xr.ufuncs.isnan(geo_grid),0,drop=False)
    return geo_grid


#%% Creat isglac and isout maps
#Isglac is a raster with a value of 1 for all glaciated pixels within the basin
#Isout is a raster with a value of 1 for all glaciated pixels that lie for >50% outside the basin boundary

if FIND_ISOUT_ISGLAC:
    for B in basin_names:
        bbox = (Basins[B]['llc_lon'],Basins[B]['llc_lat'],
                  Basins[B]['urc_lon'],Basins[B]['urc_lat'])
        #Load global glacier fraction vectorized grid on the GHM resolution
        glac_frac = gp.read_file(Basins[B]['p']['globglacfrac'],bbox=bbox) # basin_geo = shape(basin['geometry']).buffer(0.0001) #buffer against invalid shape error
        glob_mask = gp.read_file(Basins[B]['p']['glob_mask'],
                                 bbox=bbox)  #Global vectorized landmask on the GHM resolution
    
        print (B+' isout & isglac')
        basin_geo = Basins[B]['subshp'].geometry[0]
        # basin_geo = basin_geo.buffer(-0.001)
        
        #append all grid cells with glaciers in it
        isglac_idx = [] 
        #append all grid cells with glaciers & grid cell area >50% outside basin boundaries
        isout_idx=[]   
        for i in range(len(glac_frac)):
            cell = glac_frac.geometry[i]
            if cell.intersects(basin_geo):
                isglac_idx.append(i)
            # if cell.overlaps(basin_geo):
            #     isout_idx.append(i) ##wrong
            # #
            if cell.overlaps(basin_geo):
                diff = cell.difference(basin_geo)
                if (diff.area/cell.area)>0.5:
                    isout_idx.append(i)

        isglac_vec = glac_frac.loc[isglac_idx]
        isout_vec  = glac_frac.loc[isout_idx]
        
        isbasin_idx = []
        for i in range(len(glob_mask)):
            cell=glob_mask.geometry[i]
            # if basin_geo.contains(cell):
            if cell.intersects(basin_geo):
                diff = cell.difference(basin_geo)
                if (diff.area/cell.area)<0.5:
                    isbasin_idx.append(i)
        isbasin_vec = glob_mask.loc[isbasin_idx]
        
        f1,ax1=plt.subplots(figsize=(10,7))
        isbasin_vec.plot(ax=ax1,color='tab:red',alpha=0.7,label='Basin grid')
        isglac_vec.plot(ax=ax1,alpha=0.7,color='tab:blue',label='Glacier grid')
        isout_vec.plot(ax=ax1,color='black',label='Spill grid')
        Basins[B]['subshp'].plot(ax=ax1,facecolor='None',edgecolor='tab:blue',label='Basin shapefile')
        # ax1.legend(['a','b','c','d'])
        ax1.set_title(B.title())
        
        
        #Rasterize isglac & isout

        isglac = rasterize(isglac_vec,'glac_frac').glac_frac
        isout  = rasterize(isout_vec,'bool').bool
        isbasin=rasterize(isbasin_vec,'bool').bool
        
        Basins[B]['isglac'] = isglac
        Basins[B]['isout']  = isout
        Basins[B]['isbasin']=isbasin
        
        isout.to_netcdf(Basins[B]['p']['isout'])
        isglac.to_netcdf(Basins[B]['p']['isglac'])
        isbasin.to_netcdf(Basins[B]['p']['isbasin'])
        
        #Save as raster tif
        # isglac.rename({'lat':'y','lon':'x'}).rio.to_raster(Basins[B]['p']['isglac'])
        # isout.rename({'lat':'y','lon':'x'}).rio.to_raster(Basins[B]['p']['isout'])
else:
    for B in basin_names:
        with xr.open_dataset(Basins[B]['p']['isout']) as isout_,\
             xr.open_dataset(Basins[B]['p']['isglac']) as isglac_,\
                xr.open_dataset(Basins[B]['p']['isbasin']) as isbasin_:
            Basins[B]['isglac'] = isglac_.glac_frac
            Basins[B]['isout']  = isout_.bool
            Basins[B]['isbasin']=isbasin_.bool


#%% Establish dates
#Q_months is made with an additional month at the start,
# otherwise xr.resample.backfill() doesn't work
for B in basin_names:
    if Basins[B]['llc_lat']>0:
        startdate   =str(GG_FROMYEAR-1)+'-10-01'
        enddate     =str(GG_UNTILYEAR)+'-09-30'
        Basins[B]['Q_months']    =pd.date_range(str(GG_FROMYEAR-1)+'-09',
                                           str(GG_UNTILYEAR)+'-10',
                                           freq='M') #.strftime('%Y-%m')
        Basins[B]['Q_days']      = pd.date_range(str(GG_FROMYEAR-1)+'-09-01',
                                    str(GG_UNTILYEAR)+'-09-30')
    else: 
        startdate   =str(GG_FROMYEAR-1)+'-04-01'
        enddate     =str(GG_UNTILYEAR)+'-03-31'
        Basins[B]['Q_months']    =pd.date_range(str(GG_FROMYEAR-1)+'-03',
                                   str(GG_UNTILYEAR)+'-04',
                                   freq='M')#.strftime('%Y-%m')
        Basins[B]['Q_days']      = pd.date_range(str(GG_FROMYEAR-1)+'-03-01',
                                    str(GG_UNTILYEAR)+'-03-31')
        
    Basins[B]['daterange'] =pd.date_range(startdate,enddate) 
    

#%% Load GloGEM timeseries
for B in basin_names:
    print (B,' loading GloGEM data')
    #Southern Andes and LowLatitudes are corrupted in some way
    #So i manually converted the original space delimited to comma seperated using excel
    if (Basins[B]['RGI_dir'] == 'RGI17_SouthernAndes') or \
        (Basins[B]['RGI_dir'] == 'RGI16_LowLatitudes'):    delimiter=','
    else:    delimiter = '   '
    
    Q_catch      = pd.read_csv(Basins[B]['p']['Qcatch'],
                        skiprows=[0],
                        usecols=[0,*range((GG_FROMYEAR-1980)*12,(GG_UNTILYEAR+1-1980)*12+1)],
                        header=None,
                        delimiter=delimiter,
                        skipinitialspace=True,
                        index_col=0,
                        engine='python', #not sure why this option is needed
                        dtype=float)
    Q_catch     =Q_catch.transpose()
    Q_catch     =Q_catch.set_index(Basins[B]['Q_months'])
    
    Q_catch  = Q_catch.loc[:,Basins[B]['GG_points'].index] # in m
    Basins[B]['Q_volume'] = (Basins[B]['GG_area0']*1e6)*Q_catch #in m3

#%% Rasterize GloGEM from point data
if FIND_QM:
    # def rasterize(month,Basins,B,res = GHM_RESOLUTION):
    #     """ This function converts the monthly runoff data per individual glacier into 
    #     monthly runoff grid data with Geocube 
    #     """
    #     gdf        = gp.GeoDataFrame({'geometry':Basins[B]['GG_points'],
    #             'R':Basins[B]['Q_volume'][month.strftime('%Y-%m')].values[0]})
    #     #Perform rasterize in positive coordinates and then flip back
    #     lon360   = 360 *(Basins[B]['llc_lon']<0)
    #     latflip   = 1-2* (Basins[B]['llc_lat']<0)
    #     gdf.geometry = gdf.geometry.affine_transform([1,0,0,latflip,lon360,0])    
    #     geo_grid    = make_geocube(vector_data = gdf,
    #             measurements = ['R'],
    #             resolution = (res,res),
    #             geom = json.dumps(mapping(box(Basins[B]['llc_lon']+lon360,
    #                                           Basins[B]['llc_lat']*latflip,
    #                                           Basins[B]['urc_lon']+lon360,
    #                                           Basins[B]['urc_lat']*latflip))),
    #             rasterize_function=partial(rasterize_image,merge_alg = MergeAlg.add),
    #             fill = 0)
    #     geo_grid['x'] = geo_grid['x']-lon360
    #     geo_grid['y'] = geo_grid['y']*latflip
    #     return geo_grid
    
    

    for B in basin_names:
        print ('Rasterize '+B)
        def find_gdf(month):
            return gp.GeoDataFrame({'geometry':Basins[B]['GG_points'],
                'R':Basins[B]['Q_volume'][month.strftime('%Y-%m')].values[0]})
        Qlist       = [rasterize(find_gdf(m),'R','add') for m in Basins[B]['Q_months']]
        # Qlist       = [rasterize(m,Basins,B) for m in Basins[B]['Q_months']]
        
        
        Qm          = xr.concat(Qlist,dim=pd.Index(Basins[B]['Q_months'],name='time'))
        
        
        
        #Mass conservation check
        pointsum  = Basins[B]['Q_volume'].sum().sum()
        rastersum = Qm.R.sum().data
        print ('Monthly rastersum is '+str(round((rastersum/pointsum)*100,3))+' % of pointsum')
        
        Qm.to_netcdf(Basins[B]['p']['Qm'])
        Qm.close()
else:
    for B in basin_names:
        if not os.path.isfile(Basins[B]['p']['Qm']):
            print ('No Qm for '+B+' this daterange')
    

#%% Resampling of GloGEM data from monthly to daily resolution

if FIND_QD:
    for B in basin_names:
        print ('Resampling '+B)
        #Q_catch is given at the last day of each month, backfill gives the same value to the whole month
        with xr.open_dataset(Basins[B]['p']['Qm']) as Qm:
            Qm_=Qm
        # Qm_          = xr.open_dataset(Basins[B]['p']['Qm_'])
        Qd          = Qm_.resample(time='D').bfill() 
        Qd          = Qd.sel(time=Basins[B]['daterange']) # Get rid of extra months at start and finish
        Qd          = Qd.rename({'x':'lon','y':'lat'})
        
        ### Setup dates and temperature series
        #Hydrolical year is from October 1st to September 30th (april 1st to March 31st in Southern Hemisphere (6 months earlier))
        Q_mlen      = [m.day for m in Basins[B]['Q_months']] #Lengths of months over the daterange
        Q_halfidx   = np.cumsum(Q_mlen)-16      #Index of Dates halfway the months
        Q_mhalf     = Basins[B]['Q_days'][Q_halfidx]         #Dates halfway the months
        
        
        if ERA_RESAMPLE == True: #Resample with a temperature weight function
        
            ### LOAD ERA-INTERIM DATA
            ERA_full        =xr.open_dataset(Basins[B]['p']['ERA_tas']).load()
            ERA_full['lon'] = xr.where(ERA_full.lon>180,ERA_full.lon-360,ERA_full.lon)
            ERA_full['lon_bnds'] = xr.where(ERA_full.lon_bnds>180,ERA_full.lon_bnds-360,ERA_full.lon_bnds)
            ERA                 = ERA_full.sel(time=Basins[B]['daterange'])
            
            if not np.all(Qd.time==ERA.time):
                print ('Dates do not match!')        
                
            w_15    = xr.zeros_like(ERA.tas)
            day     = 0
            #Loop through months, D=number of days in month
            
            t=time.time()
            for D in Q_mlen:   #nvm #[1:-1]: #Get rid of extra months at start and finish
                for lats in ERA.lat:
                    for lons in ERA.lon:
                        T = ERA.isel(time=slice(day,day+D)).sel(lat=lats,lon=lons).tas
                        #Weights excluding subzero days
                        positive_days = T>T_THRESHOLD
                        pd_count = np.count_nonzero(positive_days)
                        if pd_count==0:
                            w_15.loc[w_15.time[day:day+D],lats,lons] = 1/D#to make weights sum to 1 and not lose any melt on subzero days
                        else:
                            Tp      = T.where(positive_days,np.nan)
                            Tm      = Tp.mean()
                            Tsum    = Tp.sum()
                            # w = (1+WEIGHT_FACTOR*(x-xmean)/xmean)n
                            w_15.loc[w_15.time[day:day+D],lats,lons] = \
                                xr.where(positive_days,
                                          (1+ WEIGHT_FACTOR*(Tp-Tm)/Tm)/pd_count,
                                          0)
                day+=D
            print (time.time()-t,'seconds have passed in resampling loop')
            
            #Upscale to 5 arcmin 
            w_05    = w_15.interp_like(Qd,method='nearest',
                                        kwargs={'fill_value':'extrapolate'})
            
            if not np.all(np.round(w_05.sum(axis=0))==len(Basins[B]['Q_months'])-1):
                print('Weights do not add up to 1!') 
                
                
        else:
            #Simply divide by number of days in month
            w_05 =xr.zeros_like(Qd.R)
            day     = 0
            for D in Q_mlen: #nvm the [1:]
                w_05.loc[w_05.time[day:day+D]]=1/D
                day+=D
            if not np.all(np.round(w_05.sum(axis=0))==len(Basins[B]['Q_months'])-1):
                print('Weights do not add up to 1!')         
        
        Qd    = Qd*w_05
        
        pointsum  = Basins[B]['Q_volume'][1:].sum().sum().data
        rastersum = Qd.R.sum().data
        print ('Rastersum is '+str(round((rastersum/pointsum)*100,3))+' % of pointsum')
        
        Qm_.close()
        Qd.to_netcdf(Basins[B]['p']['Qd'])
else:
    for B in basin_names:
        if not os.path.isfile(Basins[B]['p']['Qd']):
            print ('No Qd for '+B+' this daterange')

    
    
#%% Spilling prevention
# Take all runoff from isout grid cells and transfer to nearby safe grid cells
if SPILLING_PREVENTION:
    for B in basin_names:
        print (B+' spilling prevention')
        with xr.open_dataset(Basins[B]['p']['Qd']) as Qd_:
            Qd=Qd_
        Qd = Qd.where(~xr.ufuncs.isnan(Qd),0) #Crucial to remove nans
        
        R_out=0 #Total runoff that would've spilled
        spilled = 0
        
        #Way too complicated loop to make sure all runoff falls inside basin
        for i in range(Qd.dims['lat']):
            for j in range(Qd.dims['lon']):
                lats = Qd.lat[i]
                lons = Qd.lon[j]
                b = Basins[B]['isout'].sel(lon=lons,lat =lats,method='nearest',tolerance=0.05)
                SUCCESS=False
                if b==1: #runoff falls outside basin
                    for ii in [-1,0,1]:
                        if (i+ii<-1)or(i+ii>=Qd.dims['lat']):
                            continue
                        for jj in [-1,0,1]:
                            if (j+jj<-1)or(j+jj>=Qd.dims['lon']):
                                continue
                            latss = Qd.lat[i+ii]
                            lonss = Qd.lon[j+jj]
                            c = Basins[B]['isout'].sel(lon=lonss,lat =latss,method='nearest',tolerance=0.01)
                            isgl = Basins[B]['isglac'].sel(lon=lonss,lat =latss,method='nearest',tolerance=0.01)
                            isbasin=Basins[B]['isbasin'].sel(lon=lonss,lat =latss,method='nearest',tolerance=0.01)
                            if (isgl>0)&(c==0):
                                Qd.R.loc[:,latss,lonss]=Qd.R.loc[:,latss,lonss] +Qd.R.loc[:,lats,lons]
                                R_out    =R_out+ Qd.R.loc[:,lats,lons].sum()
                                Qd.R.loc[:,lats,lons]=0
                                SUCCESS=True
                                break            
                            if ii==jj==1:
                                spilled = spilled + Qd.R.loc[:,lats,lons].sum()
                                print ('Runoff falls out of basin at',latss.data,lonss.data)
                        if SUCCESS:break
        
        R_total = Qd.R.sum() 
        outfrac = R_out/R_total
        spillfrac = spilled/R_total
        print (str(np.round(float(outfrac)*100,3))+'% of the total runoff would fall outside basin and is thrown back in')
        print (str(np.round(float(spillfrac)*100,3))+'% of the total runoff is spilled')
        
        pointsum  = Basins[B]['Q_volume'][1:].sum().sum()
        rastersum = Qd.R.sum().data
        print ('Rastersum is '+str(round((rastersum/pointsum)*100,3))+' % of pointsum')
        
        basin_info.loc[B,'pct_spilled']    =outfrac.data
        
        #Create output for Aletsch (&Rhone) glacier observation comparison
        if B=='RHONE': 
            Aletsch_Q = Qd.R[:,41:43,60].sum(axis=1)
            Aletsch_Q.to_netcdf(join(files,'Aletsch_Q_'+str(WEIGHT_FACTOR)+'.nc'))
            Rhoneglac_Q = Qd.R[:,42,64]
            Rhoneglac_Q2=Qd.R[:,41,63]
            Rhoneglac_Q.to_netcdf(join(files,'Rhoneglac_Q_'+str(WEIGHT_FACTOR)+'.nc'))
            Rhoneglac_Q2.to_netcdf(join(files,'Rhoneglac_Q2_'+str(WEIGHT_FACTOR)+'.nc'))

        Qd.to_netcdf(Basins[B]['p']['Qsp'])
    basin_info.to_csv(join(files,'basin_info_45min.csv'))
    # 
#%% Add Hydrographs to the netcdf files
for B in basin_names:
    if SPILLING_PREVENTION:
        with xr.open_dataset(Basins[B]['p']['Qsp']) as Qfinal_:
            Qfinal=Qfinal_
    else: 
        with xr.open_dataset(Basins[B]['p']['Qd']) as Qfinal_:
            Qfinal=Qfinal_    
    if INCLUDE_HG:
        print (B+' load basin runoff observations')
        #Load Basin discharge observations

            
            
        if B=='RHONE':
            beau_total  = pd.read_csv(Basins[B]['p']['Beaucaire']
                                      ,index_col=0,
                                      parse_dates = True)
            tarr_total = pd.read_csv(Basins[B]['p']['Tarrascon']
                                      ,index_col=0,
                                      parse_dates = True)
            beautarr = pd.concat([beau_total[:'2009'],tarr_total['2010':]])
            beautarr = beautarr.rename(columns={'0':'Q'})
            hg = beautarr[str(GG_FROMYEAR-1)+'-10':str(GG_UNTILYEAR)+'-09']
            Qfinal['hg']=(('time'),hg['Q'])
            Qfinal.hg.attrs['long_var_name'] = 'Hydrograph of observations downstream'
            Qfinal.hg.attrs['unit'] = 'm3/s'
            Qfinal.hg.attrs['source']='Hydrobanque France'
            Qfinal.hg.attrs['station']='Beaucaire 1980-2009, Tarrascon 2010-2020'
        elif B =='INDUS':
            tarbela = pd.read_csv(Basins[B]['p']['Tarbela'],
                                  index_col=0,
                                  parse_dates=True,
                                  skiprows=1,dayfirst=True)
            tarbela=tarbela.to_xarray()
            Qfinal = xr.merge([Qfinal,tarbela])
            Qfinal.hg.attrs['long_var_name'] = 'Hydrograph of observations downstream'
            Qfinal.hg.attrs['unit'] = 'm3/s'
            Qfinal.hg.attrs['source']='WAPDA'
            Qfinal.hg.attrs['station']='Tarbela dam inflow'
        else:
            stations=1
            for path in Basins[B]['p']['GRDC']:
                hg = pd.read_csv(path,
                                 delimiter=';',
                                 skiprows=36,
                                 index_col=0,
                                 parse_dates=True,
                                 skipinitialspace=True,
                                 usecols=[0,2])
                hg = hg.rename(columns  = {'YYYY-MM-DD':'Date','Value':'Q'})
                hg_full = pd.DataFrame(data=hg,index=Basins[B]['daterange'])
                if stations==1:
                    Qfinal['hg']=(('time'),hg_full['Q'])
                elif stations==2:
                    Qfinal['hg2']=(('time'),hg_full['Q'])
                Qfinal.hg.attrs['observation_date_range'] = (dt.datetime.strftime(hg.index[0].date(),'%Y-%m-%d'),\
                    dt.datetime.strftime(hg.index[-1].date(),'%Y-%m-%d'))
                Qfinal.hg.attrs['long_var_name'] = 'Hydrograph of observations downstream'
                Qfinal.hg.attrs['unit'] = 'm3/s'
                Qfinal.hg.attrs['source']='GRDC'
                Qfinal.hg.attrs['station']=os.path.basename(path).split('_')[1]
                stations+=1
        Qfinal.attrs['station_coords1'] = (Basins[B]['station_lon1'],
                                              Basins[B]['station_lat1'])
        # if ~np.isnan(Basins[B]['station_lon2']):
        #     Qfinal.attrs['station_coords2'] = (Basins[B]['station_lon2'],
        #                                       Basins[B]['station_lat2'])
    if INCLUDE_SP_AUX:
        Qfinal['isout']     =(('lat','lon'),Basins[B]['isout']) #Needs to be DataArray
        Qfinal['isglac']    =(('lat','lon'),Basins[B]['isglac'])
        Qfinal['isbasin']   =(('lat','lon'),Basins[B]['isbasin'])
    Qfinal.to_netcdf(Basins[B]['p']['Qfinal'])

#%% Add attributes and save 
if SAVE_NC:
    for B in basin_names:
        with xr.open_dataset(Basins[B]['p']['Qfinal']) as Qfinal_:
            Qfinal=Qfinal_
        Qfinal.attrs['llc']=(Basins[B]['llc_lon'],Basins[B]['llc_lat'])
        Qfinal.attrs['urc']=(Basins[B]['urc_lon'],Basins[B]['urc_lat'])
        Qfinal.attrs['Basin name'] = B
        Qfinal.attrs['resolution'] = GHM_RESOLUTION
        Qfinal.R.attrs['unit'] = 'm3/day'
        Qfinal.R.attrs['long_var_name']='Glacier Runoff from GloGEM'
        Qfinal.R.attrs['glaciaction_degree'] = Basins[B]['glac_degree']
        # Qfinal.to_netcdf(Basins[B]['p']['Qfinal'])
        Qfinal.to_netcdf(join(
            Basins[B]['p']['nc_out'],'_'.join([B,
                                        str(GG_FROMYEAR),
                                        str(GG_UNTILYEAR),
                                        RUN_NAME,
                                        'R.nc'])))


#%% Create worldmap with all basins 
# OF30 = pd.read_csv(join(r'd:\Documents\Master vakken\Thesis\Code\Output',
#                         'OF29.csv'),index_col=0)

if CREATE_WORLDMAP:
    
    FIGSIZE=10
    SHAPEX=0.6
    f1 = plt.figure(figsize=(FIGSIZE,SHAPEX*FIGSIZE))
    
    OCEANALPHA=0.02
    ax1 = f1.add_axes(projection=ccrs.PlateCarree(),
                         rect = [0,0,1,1])
    ax1.set_extent((-120,118,-60,85),crs=ccrs.PlateCarree())
    ax1.coastlines()
    ax1.add_feature(cfeature.LAND,color='green',alpha=0.00)
    ax1.add_feature(cfeature.OCEAN,color='blue',alpha=OCEANALPHA)
    
    #NEW-ZEALAND subaxes
    ax2 = f1.add_axes(projection=ccrs.PlateCarree(),
                         rect=[0.762,0.00,0.3,0.4])
    ax2.set_extent((160,180,-60,-30),crs=ccrs.PlateCarree())
    ax2.coastlines()
    ax2.add_feature(cfeature.LAND,color='green',alpha=00)
    ax2.add_feature(cfeature.OCEAN,color='blue',alpha=OCEANALPHA)

    
    #North-America subaxes
    ax3 = f1.add_axes(projection=ccrs.PlateCarree(),
                 rect=[0.007,0.519,0.39,0.51])
    ax3.set_extent((-179,-100,30,85),crs=ccrs.PlateCarree())
    ax3.coastlines()
    ax3.add_feature(cfeature.LAND,color='green',alpha=00)
    ax3.add_feature(cfeature.OCEAN,color='blue',alpha=OCEANALPHA)
    
    #Iceland subaxes
    ax4 = f1.add_axes(projection=ccrs.PlateCarree(),
          rect=[0.470,0.0,0.4,0.3])
    ax4.set_extent((-28,-10,60,70),crs=ccrs.PlateCarree())
    ax4.coastlines()
    ax4.add_feature(cfeature.LAND,color='green',alpha=00)
    ax4.add_feature(cfeature.OCEAN,color='blue',alpha=OCEANALPHA)

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
    arrowxy = {'AMAZON':[-3,-3],
           'IRRAWADDY':[2,7],
           'MACKENZIE':[2,0],
           'YUKON':[-4,10],
           'ALSEK':[-12,-7],
           'CLUTHA':[3,-3],
           'COLUMBIA':[0,-8],
           'COPPER':[-11,-5],
           'DRAMSELV':[-14,4],
           'FRASER':[-12,-9],
           'GLOMA':[-5,7],
           'KUSKOKWIM':[-16,-3],
           'NASS':[-17,-12],
           'NEGRO':[0,6],
           'OELFUSA':[-3,-1],
           'RHINE':[6,2],
           'RHONE':[-5,-10],
           'SKAGIT':[-10,-9],
           'SKEENA':[-8,-7.5],
           'STIKINE':[-17,-10],
           'SUSITNA':[-18,-8],
           'TAKU':[-15,-8.5],
           'THJORSA':[-1,-2],
           'OB'     :[-3,-3],
           'DANUBE': [8,5]}
    

    cmap   = plt.cm.Blues
    VMIN =0
    import matplotlib.colors as mcolors
    
    
    norm = mcolors.Normalize(vmin=VMIN, vmax=20)
    

    
    for B in basin_names:
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
        ax.annotate(B.title(),xy,
                     xy+np.array(arrowxy[B])*1.3,
                arrowprops=dict(facecolor='black',
                                width=0.001,headlength=0.02,
                                headwidth=0.02,alpha=0.7))
    
   
    
    # ax1.set_title('25 basins with observations')
    
    ax5 = f1.add_axes([-0.25,0.05,0.3,0.4])
    ax5.set_visible(False)
    im = ax5.imshow(np.array([[VMIN,20]]),cmap=cmap)
    bar =plt.colorbar(im,fraction=0.035,label='Glaciation degree [%]')
    
    # red_dot = [Line2D([0],[0],marker='o',linestyle='None',
    #            markerfacecolor='red',markeredgecolor='red',
    #            markersize=10,label='Glaciers')]
    # plt.legend(handles = red_dot)
    
    
    save_at = join(run_dir,r'Code\Figures\Basin_maps','worldmap25_small.svg')
    plt.savefig(save_at,bbox_inches='tight',format = 'svg')
    plt.show()
