# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 16:27:45 2021

@author: Pau Wiersma

Scripts accompanying PCR-GLOBWB 2 and GloGEM coupling
(2/3): ewatercycle_runs



This script:
    - Loads the NetCDF-files of (1/3): Preprocessing and couples them to PCR-GLOBWB
    - Loads and adjusts the PCRGLOB inifile
    - Calls the PCRGLOB docker and initializes the model
    - Has a setting for doing only spinup
    - Runs PRCGLOB coupled or uncoupled for one basin
    - Saves modeled discarge at gauging station 
    - Saves daily parameter fields into netcdf


    - Can be made into a function and run together with eval_inputs.py 
        for job arrays
    - 

Files needed:
    NetCDF file from preprocessing
    setup_05min_non-natural.ini
    clone_global_05min.map
    ERA-Interim T&P nc files
    Adjusted landcover maps!

    
To do:
    Make logging optional for use with eval_inputs.py
    Make script function optional
    Remove Proportional settings?
    Station 2? Is it used in preprocessing?
    daily_outputs too large?
"""
from ewatercycle.parametersetdb import build_from_urls
from os.path import join 
import subprocess
import numpy as np
import time
import datetime
import pandas as pd
import xarray as xr
import sys
import logging
import os
# import rasterio
# import glob
import os
# from glaciers_ewcfunc import *
import sys

#%% Inputs 
run_dir = r'/media/pau/FREECOM/Thesis/Data/PCRGLOB/Global/'
# run_dir       = r'/projects/0/wtrcycle/users/pwiersma/pcrg_glaciers/'
#run_dir = r'/scratch-shared/pwiersma/pcrg_glaciers/'
os.chdir(run_dir)


Basin_name  = 'RHONE'

run_name = 'test'
nc_run_name = 'test'

#model_setup = 3
# model_setup  =0 #int(sys.argv[2])
setup = 0
# 0 = default model
# 1 = no coupling, Grasslands landcover
# 2 = coupling, grasslands landcover
# 3 = no coupling, proportional landcover
# 4 = coupling, proportional landcover

#Test_run = True or False
# test_run_string = sys.argv[3]
# test_run = test_run_string=='True'
test_run = True

#HH2018 inputs
GG_GCM              = 'HadGEM2-ES'  
GG_rcp              = 'rcp26'
GG_fromyear         =2001
GG_untilyear        = 2002

spinup=True




#%%Functions


def day2date(day):
    date = datetime.date(1901,1,1)+datetime.timedelta(days=day)
    return date

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


#%% Make function of script for eval_inputs.py
#Validaiton function
# def ewc_val(Basin_name,
#             setup,
#             GG_fromyear,
#             GG_untilyear,
#             run_dir,
#             spinup,
#             test_run,
#             GG_GCM = 'HadGEM2-ES',
#             GG_rcp = 'rcp26'):
    
#%%

timer1      = time.time()
run_name = '_'.join(['PCRG',
                        Basin_name,
                        str(setup),
                        str(GG_fromyear),
                        str(GG_untilyear),
                        run_name])


#%% Setup logging in case eval_inputs.py is used

# logging.basicConfig(filename=join(r'/nfs/home3/pwiersma/pcrg_glaciers','monitoring'
#                                   ,Basin_name+str(setup)+'.log'),
#                     level=logging.INFO)


#%% Prepare runs
print(time.strftime('%Y-%m-%d %H:%M'))
print('Start of '+Basin_name+' with setting '+str(setup))
print('Runtime: '+str(GG_fromyear)+' - '+str(GG_untilyear))
  
# Adjust settings depending on model setup
if setup==0:
    couple_GG       =False
    adjust_landcov   =False
elif setup==1:
    couple_GG       =False
    adjust_landcov  = 'GRASSLANDS'
elif setup==2: 
    couple_GG       = True
    adjust_landcov  ='GRASSLANDS'
elif setup==3:
    couple_GG       =False
    adjust_landcov  ='PROPORTIONAL'
elif setup==4:
    couple_GG       =True
    adjust_landcov  ='PROPORTIONAL'
else: 
    print('Not a valid model setup')
    sys.exit()


# Prepare glacier input

nc_path = join(run_dir,'glaciers_ncfiles','_'.join([Basin_name,
                                    '2000',
                                    '2001',
                                    nc_run_name,
                                    'R.nc'])) 

nc = xr.open_dataset(nc_path)


# station2   = False
# if 'station_coords2' in nc.attrs.keys():
#   station2 = True

# isglac = np.any(nc.R>0,axis=0) #create xarray with ones for glacier pixels
  


### load pmset

glob_ini = join(run_dir,'glaciers_inifiles','setup_05min_non-natural.ini')
glob_ini_new = join(run_dir,'glaciers_inifiles','setup_05min_non_natural_'+Basin_name+str(setup)+'.ini')

parameter_set = build_from_urls(config_format = 'ini', 
                                config_url = 'file://'+glob_ini,
                                  datafiles_format='svn', 
                                  datafiles_url='file://'+run_dir)

# parameter_set.save_datafiles('input1') #werkt niet
### make clonemap
print('make clonemap and setup config')
clone_path  = join(run_dir,'global_05min/cloneMaps')
orig_clone  = 'clone_global_05min.map'
# dest_clip   = Basin_name+'05min_clone.map'
dest_clip   = Basin_name+'_05min.map'


ullon,lrlat     = nc.llc.astype(float).astype(str)
lrlon,ullat     = nc.urc.astype(float).astype(str)
nlat,nlon       = nc.dims['lat'],nc.dims['lon']

command     = 'gdal_translate -of PCRASTER -projwin '+ullon+' '+ullat+' '+lrlon+' '+lrlat+' '\
    +join(clone_path,orig_clone)+' '\
    + join(clone_path,dest_clip)
subprocess.Popen(command,shell=True)

if int(float(lrlat))>0:
  startT   = str(GG_fromyear-1)+'-01-01' #From 1st of January to catch on Spinup
  endT     = str(GG_untilyear)+'-09-30'
else:
  startT   = str(GG_fromyear-1)+'-01-01'
  endT     = str(GG_untilyear)+'-03-31'      

# setup config
precip_name = 'pcrglobwb_OBS6_ERA-Interim_reanaly_1_day_pr_1999-2016_'+Basin_name+'.nc'
temp_name  ='pcrglobwb_OBS6_ERA-Interim_reanaly_1_day_tas_1999-2016_'+Basin_name+'.nc'

#input en output_dir moeten /data/input en /data/output zijn
parameter_set.config['globalOptions']['inputDir']       = '/data/input'
parameter_set.config['globalOptions']['outputDir']      = '/data/output'
parameter_set.config['globalOptions']['cloneMap']       = join(r'global_05min/cloneMaps',dest_clip)
parameter_set.config['globalOptions']['startTime']      = startT
parameter_set.config['globalOptions']['endTime']        = endT
parameter_set.config['meteoOptions']['precipitationNC'] = join(r'global_05min/meteo',precip_name)
parameter_set.config['meteoOptions']['temperatureNC']   = join(r'global_05min/meteo',temp_name)
parameter_set.config['meteoOptions']['referenceETPotMethod'] = 'Hamon'
parameter_set.config['meteoOptions']['refETPotFileNC']  = 'None'
parameter_set.config['routingOptions']['routingMethod']             ='kinematicWave' #not kinematicWave
parameter_set.config['routingOptions']['dynamicFloodPlain']         ='False' #True doesn't work
parameter_set.config['meteoOptions']['precipitationVariableName']   = 'pr'
parameter_set.config['meteoOptions']['temperatureVariableName']     = 'tas'  #temperature at surface
parameter_set.config['meteoOptions']['referenceEPotVariableName']   ='evspsblpot'
parameter_set.config['meteoOptions']['temperatureConstant']         = '-273.15' #check if this is necessary
parameter_set.config['meteoDownscalingOptions']['downscalePrecipitation']='True'

if adjust_landcov!=False:
    if adjust_landcov=='GRASSLANDS':
        frac_prefix = 'grassf'
    elif adjust_landcov=='PROPORTIONAL':
        frac_prefix = 'propf'
    else:
        print('Wrong input adjust_landcov')
#     parameter_set.config['landSurfaceOptions']['noLandCoverFractionCorrection']='True'
#     parameter_set.config['forestOptions']['fracVegCover']   ='global_05min/landSurface/landCover/naturalTall/'+frac_prefix+'_tall.map'
#     parameter_set.config['grasslandOptions']['fracVegCover']='global_05min/landSurface/landCover/naturalShort/'+frac_prefix+'_short.map'
#     parameter_set.config['irrPaddyOptions']['fracVegCover'] ='global_05min/landSurface/landCover/irrPaddy/'+frac_prefix+'_pad.map'
#     parameter_set.config['irrNonPaddyOptions']['fracVegCover']='global_05min/landSurface/landCover/irrNonPaddy/'+frac_prefix+'_nonpad.map'
    parameter_set.config['landSurfaceOptions']['noLandCoverFractionCorrection']='True'
    parameter_set.config['forestOptions']['fracVegCover']   ='glaciers_landcov/'+frac_prefix+'_tall.map'
    parameter_set.config['grasslandOptions']['fracVegCover']='glaciers_landcov/'+frac_prefix+'_short.map'
    parameter_set.config['irrPaddyOptions']['fracVegCover'] ='glaciers_landcov/'+frac_prefix+'_pad.map'
    parameter_set.config['irrNonPaddyOptions']['fracVegCover']='glaciers_landcov/'+frac_prefix+'_nonpad.map'
elif adjust_landcov==False:
    parameter_set.config['landSurfaceOptions']['noLandCoverFractionCorrection']='False'
    parameter_set.config['forestOptions']['fracVegCover']   ='global_05min/landSurface/landCover/naturalTall/vegf_tall.map'
    parameter_set.config['grasslandOptions']['fracVegCover']='global_05min/landSurface/landCover/naturalShort/vegf_short.map'
    parameter_set.config['irrPaddyOptions']['fracVegCover'] ='global_05min/landSurface/landCover/irrPaddy/fractionPaddy.map'
    parameter_set.config['irrNonPaddyOptions']['fracVegCover']='global_05min/landSurface/landCover/irrNonPaddy/fractionNonPaddy.map'
else:
    print('Wrong input adjust_landcov')

parameter_set.save_config(glob_ini_new)

if spinup==True:
    import fileinput
    with fileinput.FileInput(glob_ini_new, inplace=True) as file:
        for line in file:
            print(line.replace('initialConditions','initialConditions/'+Basin_name+'_spinup'), end='')
    

print('call docker')
time.sleep(5)
if test_run==True:
    ti='test'+Basin_name+str(time.clock())
else:
    ti='_'+run_name
output_dir1  =join(run_dir,'output_maps','output'+ti)
if 'FREECOM' in run_dir:
    from grpc4bmi.bmi_client_docker import BmiClientDocker
    pcrg = BmiClientDocker(image='ewatercycle/pcrg-grpc4bmi:setters', image_port=55555, 
                   input_dir=run_dir, 
                   output_dir=output_dir1)
elif ('scratch' in run_dir) or ('project' in run_dir):
    from grpc4bmi.bmi_client_singularity import BmiClientSingularity
    pcrg = BmiClientSingularity(image='ewatercycle-pcrg-grpc4bmi.sif', 
                       input_dir=run_dir, 
                       output_dir=output_dir1)

print('input: '+pcrg.input_dir)
print('output: '+pcrg.output_dir)
# initialize
start = time.time()
time.sleep(18)
print ('Initialize...')
pcrg.initialize(glob_ini_new)
print (str((time.time()-start )/60)+'minutes of initialize')
print ('\a')
# Usually 3 minutes for the Rhone







tstart  = pcrg.get_start_time()
tend    = pcrg.get_end_time()
tstep   = pcrg.get_time_step()
print('tstart = '+str(day2date(tstart)))
if (day2date(tstart)!=nc.time.to_index()[0])&(couple_GG==True):
    print('Ini-startdate and GloGEM startdate do not match!')
    print('Adjust GloGEM timespan...')
    nc = nc.sel(time=slice(day2date(tstart),day2date(tend)))
    if (day2date(tstart)!=nc.time.to_index()[0])&(couple_GG==True):
      print('Ini_startdate and GloGEM still do not match!')
      sys.exit()

#Set up variables
latsize,lonsize = pcrg.get_grid_shape(1)
if (latsize!=nc.dims['lat']) or lonsize!=nc.dims['lon']:
    print ('Shapes GloGEM & PCRG do not match !')
    print('PCRG lat,lon dimensions: '+str(latsize)+' , '+str(lonsize))
    print('GloGEM lat,lon dimensions: '+str(nc.dims['lat'])+ ' , '+str(nc.dims['lon']))
    if Basin_name =='AMAZON':
      print('AMAZON is an exception, reduce GloGEM file')
      nc = nc.where(nc.lat<6,drop=True)
      isglac =np.any(nc.R>0,axis=0) 
      nlat -= 1
      if (latsize!=nc.dims['lat']) or lonsize!=nc.dims['lon']:
        print('Dimensions still not right, exit')
        sys.exit()
    else:
      sys.exit()
#Nodes are middle of gridcells 
lon_nodes       = pcrg.get_grid_y(1)
lat_nodes       = pcrg.get_grid_x(1)
#x=lat, y=lon 

# Get hydrograph
# hg_lonlat       = nc.hg.station_coords #Beaucaire
hg_lonlat       = nc.station_coords1
hg_idx          =  (find_nearest(lat_nodes,hg_lonlat[1]),
                      find_nearest(lon_nodes,hg_lonlat[0]))
if Basin_name =='THJORSA':
    hg_idx = (9,18)
hg_idx_flat     = np.ravel_multi_index(hg_idx,(latsize,lonsize))
print('hg_idx='+str(hg_idx))

Q_station1      = []
# if station2==True:
#     hg_lonlat2       = nc.station_coords2
#     hg_idx2          =  (find_nearest(lat_nodes,hg_lonlat2[1]),
#                           find_nearest(lon_nodes,hg_lonlat2[0]))
#     hg_idx_flat2     = np.ravel_multi_index(hg_idx2,(latsize,lonsize))
    
#     Q_station2     = []
Q_date          = []
# snowsum = []

# print ('Script before pcrg run takes',time.time()-timer1,'seconds')
#Timer1 = 206 seconds = 3.5 minutes (of which 3 for initialize)

# SWE_list = []
variables = ['snow_water_equivalent','discharge','groundwater_recharge']
daily_outputs   = xr.Dataset({var:xr.zeros_like(nc.R) for var in variables })

#Start Run
i=0
start = time.time()
while pcrg.get_current_time()!=pcrg.get_end_time():
    full_loop_timer = time.time()
    
    #if spilling_prevention ==True:
        # if i>0:
        #     #chanstor = m3,  discharge = m3/s
        #     chanstordif = chanstor-pcrg.get_value('discharge')*60*60*24
        #     chanstordif = np.reshape(chanstordif,(latsize,lonsize))
        
 
    
    
    if couple_GG:
        chanstor    = pcrg.get_value('channel_storage') #Gives flat array
        chanstor    =np.reshape(chanstor,(latsize,lonsize))
        # chanstor[isglac]=nc.R.isel(time=i).data[isglac]
        # chanstor = rastertotal/30 #take daily from monthly value
        
        
        if nc.lat.data[0]>nc.lat.data[-1]:
            chanstor  += nc.R.data[i,::-1,:]
        else:
            chanstor    += nc.R.data[i,:,:]
        chanstor    = chanstor.flatten()
        pcrg.set_value('channel_storage',chanstor)
        
        # Set values with flat indices
        
        #pcrg.set_value_at_indices('channel_storage',
         #                         np.flatnonzero(isglac.data),
          #                        nc.R.isel(time=i).data[isglac])
        
    pcrg_timer  = time.time()
    ##
    pcrg.update()
    ##
    # print ('PCRG_timer:', time.time()-pcrg_timer)
    date = day2date(pcrg.get_current_time())
    print(str(date)) 
    
    #Hydrograph
    if Basin_name in ['GLOMA', 'KALIXAELVEN' ,'NASS', 'SANTA_CRUZ' ,'SKAGIT', 'STIKINE', 'THJORSA','IRRAWADDY','YUKON','OB','SUSITNA']:
        discharge_reshaped =np.reshape(pcrg.get_value('discharge'),(nlat,nlon))
        localmax = np.nanmax(discharge_reshaped[hg_idx[0]-1:hg_idx[0]+2,
                                         hg_idx[1]-1:hg_idx[1]+2])
                                         
        Q_station1.append(localmax)
        # Q_station1.append(pcrg.get_value_at_indices('discharge', hg_idx_flat)[0])
        # if station2==True:
        #     localmax2 = np.nanmax(discharge_reshaped[hg_idx2[0]-1:hg_idx2[0]+2,
        #                                  hg_idx2[1]-1:hg_idx2[1]+2])
        #     Q_station2.append(pcrg.get_value_at_indices('discharge',hg_idx_flat2)[0])
    else:
        Q_station1.append(pcrg.get_value_at_indices('discharge', hg_idx_flat)[0])
        # if station2==True:
        #     Q_station2.append(pcrg.get_value_at_indices('discharge', hg_idx_flat2)[0])
    Q_date.append(date)
    # snowsum.append(np.sum(pcrg.get_value_at_indices('snow_melt',np.flatnonzero)))
    
    #Variables: 'snow_water_equivalent', 'snow_melt'
    
    for var in variables:
        vals    = pcrg.get_value(var) #Gives flat array
        vals   =np.reshape(vals,(latsize,lonsize))
        if nc.lat.data[0]>nc.lat.data[-1]:
            vals    *= nc.isglac.data[::-1,:]
        else:
            vals    *= nc.isglac.data[:,:]
        
        # SWE_list.append(np.nansum(SWE)) 
        if nc.lat.data[0]>nc.lat.data[-1]:
            daily_outputs[var].loc[np.datetime64(date),::-1,:] = vals
        else:
            daily_outputs[var].loc[np.datetime64(date),:,:] = vals
        daily_outputs[var].attrs['unit'] = pcrg.get_var_units(var)    

    
    
    #Plot
    # f1,ax1=plt.subplots()
    # variable= 'discharge'
    # vals = pcrg.get_value(variable)
    # unit = pcrg.get_var_units(variable)

    # missval = -999.
    # Z = np.reshape(ma.masked_where(vals == np.nan, vals), (latsize,lonsize))
    # ax1.set_title(ti+'\n'+variable + '[' + unit + '], t='+str(day2date(pcrg.get_current_time())))
    # im = ax1.pcolormesh(lon_nodes,lat_nodes,Z)            
    # cbar = plt.colorbar(im,ax=ax1)
    # cbar.set_label('Discharge [m3/s]')
    
    i+=1
    # print ('Full loop timer:', time.time()-full_loop_timer)
    if (test_run==True)&(i==30):
        break
print(str((time.time()-start )/60)+'minutes of run') 
print('Maximum daily discharge: '+str(round(np.max(Q_station1),2)))
#
#f2,ax2=plt.subplots()
#ax2.plot(Q_date,Q_station1,label='Station 1')
#plt.ylabel('Discharge [m3/s]')
#plt.grid()
#if not np.all(np.isnan(nc.station_coords2)):
#    ax2.plot(Q_date,Q_station2,label='Station 2')
#plt.legend()

# timesum = SWE_nc.sum(axis=(1,2)) 
# timesum = timesum[timesum!=0]
# timesum.plot()

 #SAve hydrographs to csv
if (test_run==False):
    hg_name     = run_name+'_hg1.txt' # TAKE TIS WAYYYYYYYYYYYYY
    hg_dir      = join(run_dir,'Hydrographs',hg_name)
    hg_time     = pd.date_range(day2date(pcrg.get_start_time()),day2date(pcrg.get_current_time()))
    hg_df       = pd.DataFrame(Q_station1,index=hg_time,columns = ['hg'])
    hg_df.index.name = 'time'
    hg_df.to_csv(hg_dir)
    # if station2==True: #If no second station is given don't plot the second hg
    #     hg_name     = run_name+'_hg2.txt'
    #     hg_dir      = join(run_dir,'Hydrographs',hg_name)
    #     hg_time     = pd.date_range(day2date(pcrg.get_start_time()),day2date(pcrg.get_current_time()))
    #     hg_df       = pd.DataFrame(Q_station2,index=hg_time,columns = ['hg'])
    #     hg_df.index.name = 'time'
    #     hg_df.to_csv(hg_dir)
    daily_outputs.to_netcdf(join(run_dir,'daily_outputs',Basin_name+'_'+'_'.join([var[0] for var in variables])+'.nc'))
    
    # return

