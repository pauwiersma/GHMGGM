# -*- coding: utf-8 -*-
"""
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
    - Saves daily parameter fields into netcdf (optional)


    - Can be made into a function and run together with eval_inputs.py 
        for job arrays
    - 

to do
- Follow PIP8 naming conventions
- Mention you make a function out of this script
    Make logging optional for use with eval_inputs.py
    Make script function optional
    Station 2? Is it used in preprocessing?
    daily_outputs too large?
    Remove also the 1 behind station?

Files Needed:
    NetCDF file from preprocessing
    setup_05min_non-natural.ini
    clone_global_05min.map
    ERA-Interim T&P nc files
    Adjusted landcover maps!


- 
"""
import datetime
import os
import subprocess
# from glaciers_ewcfunc import *
import sys
import time
from os.path import join

import numpy as np
import pandas as pd
import xarray as xr
from ewatercycle.parametersetdb import build_from_urls

#%% Inputs 
BASIN_NAME  = 'RHONE'

#### Model_setup
# 0 = default model
# 1 = no coupling, Grasslands landcover
# 2 = coupling, grasslands landcover
SETUP = 2
# model_setup  =0 #int(sys.argv[2])

##### TEST_RUN = True or False
# test_run_string = sys.argv[3]
# TEST_RUN = test_run_string=='True'
TEST_RUN = True

#HussHock2018 inputs
GG_GCM              = 'HadGEM2-ES'  
GG_RCP              = 'rcp26'
GG_FROMYEAR         =2001
GG_UNTILYEAR        = 2003


## 
SPINUP=False #only changes year in inifile from 2001 to 1999



CLUSTER = False #Determines whether to use docker(False) or Singularity(cluster)

### Specifiy directory 
#run_dir = r'/scratch-shared/pwiersma/pcrg_glaciers/'
# run_dir = r'/media/pau/FREECOM/Thesis/Data/PCRGLOB/Global/'
# run_dir       = r'/projects/0/wtrcycle/users/pwiersma/pcrg_glaciers/'
RUN_DIR  = r'/home/pau/Documents/Data/PCRGLOB/Global/'
os.chdir(RUN_DIR)
#%% Functions
def day2date(day):
    """Converts "day since 1901-1-1" to an actual date"""
    date = datetime.date(1901,1,1)+datetime.timedelta(days=day)
    return date

def find_nearest(array, value): 
    """Finds the value of the nearest point within an array"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# In case you need to automatize all basins this whole script can be made into a function
# def ewc_val(BASIN_NAME,
#             SETUP,
#             GG_FROMYEAR,
#             GG_UNTILYEAR,
#             RUN_DIR,
#             SPINUP,
#             TEST_RUN,
#             GG_GCM = 'HadGEM2-ES',
#             GG_RCP = 'rcp26'):

    
#%% Setup
timer1      = time.time()
### Standard naming for runs
run_name = '_'.join(['PCRG',
                        BASIN_NAME,
                        GG_GCM,
                        GG_RCP,
                        str(SETUP),
                        str(GG_FROMYEAR),
                        str(GG_UNTILYEAR)])

# logging.basicConfig(filename=join(r'/nfs/home3/pwiersma/pcrg_glaciers','monitoring'
#                                   ,BASIN_NAME+str(SETUP)+'.log'),
#                     level=logging.INFO)
print(time.strftime('%Y-%m-%d %H:%M'))
print('Start of '+BASIN_NAME+' with setting '+str(SETUP))
print('Runtime: '+str(GG_FROMYEAR)+' - '+str(GG_UNTILYEAR))
  

if SETUP==0:
    couple_GG       =False
    adjust_landcov   =False
elif SETUP==1:
    couple_GG       =False
    adjust_landcov  = 'GRASSLANDS'
elif SETUP==2: 
    couple_GG       = True
    adjust_landcov  ='GRASSLANDS'
else: 
    print('Not a valid model setup')
    sys.exit()



### Load GloGEM runoff prepared in preprocessing script

nc_path = join(RUN_DIR,'glaciers_ncfiles','_'.join([BASIN_NAME,
                                    GG_GCM,
                                    GG_RCP,
                                    '2000',
                                    '2016',
                                    'R.nc'])) 
nc = xr.open_dataset(nc_path)

#create xarray with ones for glacier pixels
isglac = np.any(nc.R>0,axis=0) 



#%% load pmset

#Load standard .ini file
glob_ini = join(RUN_DIR,'glaciers_inifiles','setup_05min_non-natural.ini')
#Rename to make Basin+setup-specific
glob_ini_new = join(RUN_DIR,'glaciers_inifiles','setup_05min_non_natural_'+BASIN_NAME+str(SETUP)+'.ini')

parameter_set = build_from_urls(config_format = 'ini', 
                                config_url = 'file://'+glob_ini,
                                  datafiles_format='svn', 
                                  datafiles_url='file://'+RUN_DIR)

# parameter_set.save_datafiles('input1') #werkt niet


#%% make clonemap from global clonemap and setup daterange
print('make clonemap and setup config')
clone_path  = join(RUN_DIR,'global_05min/cloneMaps')
orig_clone  = 'clone_global_05min.map'
dest_clip   = BASIN_NAME+'_05min.map'


ullon,lrlat     = nc.llc.astype(float).astype(str)
lrlon,ullat     = nc.urc.astype(float).astype(str)
nlat,nlon       = nc.dims['lat'],nc.dims['lon']


## Make clip of global clonemap
command     = 'gdal_translate -of PCRASTER -projwin '+ullon+' '+ullat+' '+lrlon+' '+lrlat+' '\
    +join(clone_path,orig_clone)+' '\
    + join(clone_path,dest_clip)
subprocess.Popen(command,shell=True)



##Check if basin in Northern or Southern hemisphere and set daterange accordingly
if int(float(lrlat))>0:
  startT   = str(GG_FROMYEAR-1)+'-10-01'
  endT     = str(GG_UNTILYEAR)+'-09-30'
else:
  startT   = str(GG_FROMYEAR-1)+'-04-01'
  endT     = str(GG_UNTILYEAR)+'-03-31'      

#%% setup config
precip_name = 'pcrglobwb_OBS6_ERA-Interim_reanaly_1_day_pr_1999-2016_'+BASIN_NAME+'.nc'
temp_name  ='pcrglobwb_OBS6_ERA-Interim_reanaly_1_day_tas_1999-2016_'+BASIN_NAME+'.nc'

#input and output_dir should be /data/input and /data/output 
parameter_set.config['globalOptions']['inputDir']       = '/data/input'
parameter_set.config['globalOptions']['outputDir']      = '/data/output'
parameter_set.config['globalOptions']['cloneMap']       = join(r'global_05min/cloneMaps',dest_clip)
parameter_set.config['globalOptions']['startTime']      = startT
parameter_set.config['globalOptions']['endTime']        = endT
parameter_set.config['meteoOptions']['precipitationNC'] = join(r'global_05min/meteo',precip_name)
parameter_set.config['meteoOptions']['temperatureNC']   = join(r'global_05min/meteo',temp_name)
parameter_set.config['meteoOptions']['referenceETPotMethod'] = 'Hamon'
parameter_set.config['meteoOptions']['refETPotFileNC']  = 'None'
parameter_set.config['routingOptions']['routingMethod']             ='accuTravelTime' #not kinematicWave
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

if SPINUP==True:
    import fileinput
    with fileinput.FileInput(glob_ini_new, inplace=True) as file:
        for line in file:
            print(line.replace('initialConditions','initialConditions/'+BASIN_NAME+'_spinup'), end='')
            print(line.replace('2001','1999'), end='')

    
 
#%% Call docker and initiliaze model
time.sleep(5)
if TEST_RUN==True:
    ti='test'+BASIN_NAME+str(time.clock())
else:
    ti='_'+run_name
output_dir1  =join(RUN_DIR,'output_maps','output'+ti)


print('call docker')
if CLUSTER == True:
    from grpc4bmi.bmi_client_singularity import BmiClientSingularity
    pcrg = BmiClientSingularity(image='ewatercycle-pcrg-grpc4bmi.sif', 
                            input_dir=RUN_DIR, 
                          output_dir=output_dir1)
else: 
    from grpc4bmi.bmi_client_docker import BmiClientDocker
    pcrg = BmiClientDocker(image='ewatercycle/pcrg-grpc4bmi:setters', image_port=55555, 
                    input_dir=RUN_DIR, 
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


#%%Set up variables
tstart  = pcrg.get_start_time()
tend    = pcrg.get_end_time()
tstep   = pcrg.get_time_step()
print('tstart = '+str(day2date(tstart)))

#Adust GloGEM timerange if necessary
if (day2date(tstart)!=nc.time.to_index()[0])&(couple_GG==True):
    print('Ini-startdate and GloGEM startdate do not match!')
    print('Adjust GloGEM timespan...')
    nc = nc.sel(time=slice(day2date(tstart),day2date(tend)))
    if (day2date(tstart)!=nc.time.to_index()[0])&(couple_GG==True):
      print('Ini_startdate and GloGEM still do not match!')
      sys.exit()

#
latsize,lonsize = pcrg.get_grid_shape(1)
if (latsize!=nc.dims['lat']) or lonsize!=nc.dims['lon']:
    print ('Shapes GloGEM & PCRG do not match !')
    print('PCRG lat,lon dimensions: '+str(latsize)+' , '+str(lonsize))
    print('GloGEM lat,lon dimensions: '+str(nc.dims['lat'])+ ' , '+str(nc.dims['lon']))
    if BASIN_NAME =='AMAZON':
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

# Get hydrograph
# hg_lonlat       = nc.hg.station_coords #Beaucaire
hg_lonlat       = nc.station_coords1
hg_idx          =  (find_nearest(lat_nodes,hg_lonlat[1]),
                      find_nearest(lon_nodes,hg_lonlat[0]))


if BASIN_NAME =='THJORSA': #Thjorsa needs to be fixed manually
    hg_idx = (9,18)
hg_idx_flat   = np.ravel_multi_index(hg_idx,(latsize,lonsize))
print('hg_idx='+str(hg_idx))

Q_station1      = []
# if station2==True:
#     hg_lonlat2       = nc.station_coords2
#     hg_idx2          =  (find_nearest(lat_nodes,hg_lonlat2[1]),
#                           find_nearest(lon_nodes,hg_lonlat2[0]))
#     hg_idx_flat2     = np.ravel_multi_index(hg_idx2,(latsize,lonsize))
    
#     Q_station2     = []
Q_date          = []

# print ('Script before pcrg run takes',time.time()-timer1,'seconds')
#Timer1 = 206 seconds = 3.5 minutes (of which 3 for initialize)

SWE = []

#%% Start Run
i=0
start = time.time()
while pcrg.get_current_time()!=pcrg.get_end_time(): #Testrun=True gives another end condition
    full_loop_timer = time.time()
    if couple_GG:
        #Coupling is done by adding the glogem runoff to the channel_storage
        chanstor    = pcrg.get_value('channel_storage') #Gives flat array
        chanstor    = np.reshape(chanstor,(latsize,lonsize))
        
        # chanstor[isglac]=nc.R.isel(time=i).data[isglac]
        
        if nc.lat.data[0]>nc.lat.data[-1]: #(if Southern Hemisphere)
            chanstor  += nc.R.data[i,::-1,:]
        else:
            chanstor    += nc.R.data[i,:,:]
            
        chanstor    = chanstor.flatten()
        
        
        ####
        pcrg.set_value('channel_storage',chanstor)
        ####
        
        
        # Set values with flat indices
        #pcrg.set_value_at_indices('channel_storage',
         #                         np.flatnonzero(isglac.data),
          #                        nc.R.isel(time=i).data[isglac])
        
    pcrg_timer  = time.time()
    ##
    pcrg.update()
    ##
    # print ('PCRG_timer:', time.time()-pcrg_timer)
    print(str(day2date(pcrg.get_current_time()))) 
    
    #Hydrograph
    if BASIN_NAME in ['GLOMA', 'KALIXAELVEN' ,'NASS', 'SANTA_CRUZ' ,'SKAGIT', 'STIKINE', 'THJORSA','IRRAWADDY','YUKON','OB','SUSITNA']:
        #In these basins the station-pixel does not correspond to the river pixels
        #So the maximum discharge around the station-pixel is searched and assumed to be river
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
    Q_date.append(day2date(pcrg.get_current_time()))
    # snowsum.append(np.sum(pcrg.get_value_at_indices('snow_melt',np.flatnonzero)))
    
    #Part to save additional variables
    #Variables: 'snow_water_equivalent', 'snow_melt'    
    # for var in variables:
    #     vals    = pcrg.get_value(var) #Gives flat array
    #     vals   =np.reshape(vals,(latsize,lonsize))
    #     if nc.lat.data[0]>nc.lat.data[-1]:
    #         vals    *= nc.isglac.data[::-1,:]
    #     else:
    #         vals    *= nc.isglac.data[:,:]
    #     # SWE_list.append(np.nansum(SWE)) 
    #     if nc.lat.data[0]>nc.lat.data[-1]:
    #         daily_outputs[var].loc[np.datetime64(date),::-1,:] = vals
    #     else:
    #         daily_outputs[var].loc[np.datetime64(date),:,:] = vals
    #     daily_outputs[var].attrs['unit'] = pcrg.get_var_units(var)
        
    
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
    if (TEST_RUN==True)&(i==30):
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

 #SAve hydrographs to csv
if (TEST_RUN==False):
    hg_name     = run_name+'_hg1.txt' # 
    hg_dir      = join(RUN_DIR,'Hydrographs',hg_name)
    hg_time     = pd.date_range(day2date(pcrg.get_start_time()),day2date(pcrg.get_current_time()))
    hg_df       = pd.DataFrame(Q_station1,index=hg_time,columns = ['hg'])
    hg_df.index.name = 'time'
    hg_df.to_csv(hg_dir)
    # if station2==True: #If no second station is given don't plot the second hg
    #     hg_name     = run_name+'_hg2.txt'
    #     hg_dir      = join(RUN_DIR,'Hydrographs',hg_name)
    #     hg_time     = pd.date_range(day2date(pcrg.get_start_time()),day2date(pcrg.get_current_time()))
    #     hg_df       = pd.DataFrame(Q_station2,index=hg_time,columns = ['hg'])
    #     hg_df.index.name = 'time'
    #     hg_df.to_csv(hg_dir)
    
    # return
