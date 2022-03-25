# -*- coding: utf-8 -*-
"""
@author: Pau Wiersma

Scripts accompanying PCR-GLOBWB 2 and GloGEM coupling
(2/5): ewatercycle_runs

This script:
    - Loads the NetCDF-files of (1/3): Preprocessing and couples them to PCR-GLOBWB 2
    - Loads and adjusts the PCRGLOB inifile
    - Calls the PCRGLOB docker and initializes the model
    - Has a setting for doing only spinup
    - Runs PRCGLOB coupled or uncoupled for one basin
    - Saves modeled discarge at gauging station
    - Saves daily parameter fields into netcdf (optional)


    - Can be made into a function and run together with eval_inputs.py
        for job arrays

Files Needed:
    NetCDF file from preprocessing
    setup_05min_non-natural.ini
    clone_global_05min.map
    ERA-Interim T&P nc files
    Adjusted landcover maps
    
Output:
    Hydrographs
    Daily variable maps (optional)

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
RUN_NAME = '_'.join(['PCRG',
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
    COUPLE_GLOGEM       =False
    ADJUST_LANDCOVER   =False
elif SETUP==1:
    COUPLE_GLOGEM       =False
    ADJUST_LANDCOVER  = 'GRASSLANDS'
elif SETUP==2:
    COUPLE_GLOGEM       = True
    ADJUST_LANDCOVER  ='GRASSLANDS'
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
CLONE_PATH  = join(RUN_DIR,'global_05min/cloneMaps')
ORIG_CLONE  = 'clone_global_05min.map'
DEST_CLIP   = BASIN_NAME+'_05min.map'


UL_LON,LR_LAT     = nc.llc.astype(float).astype(str)
LR_LON,UL_LAT     = nc.urc.astype(float).astype(str)
N_NLAT,N_LON       = nc.dims['lat'],nc.dims['lon']


## Make clip of global clonemap
command     = 'gdal_translate -of PCRASTER -projwin '+UL_LON+' '+UL_LAT+' '+LR_LON+' '+LR_LAT+' '\
    +join(CLONE_PATH,ORIG_CLONE)+' '\
    + join(CLONE_PATH,DEST_CLIP)
subprocess.Popen(command,shell=True)



##Check if basin in Northern or Southern hemisphere and set daterange accordingly
if int(float(LR_LAT))>0:
    START_TIME   = str(GG_FROMYEAR-1)+'-10-01'
    END_TIME     = str(GG_UNTILYEAR)+'-09-30'
else:
    START_TIME   = str(GG_FROMYEAR-1)+'-04-01'
    END_TIME     = str(GG_UNTILYEAR)+'-03-31'

#%% setup config
PRECIP_NAME = 'pcrglobwb_OBS6_ERA-Interim_reanaly_1_day_pr_1999-2016_'+BASIN_NAME+'.nc'
TEMP_NAME  ='pcrglobwb_OBS6_ERA-Interim_reanaly_1_day_tas_1999-2016_'+BASIN_NAME+'.nc'

#input and output_dir should be /data/input and /data/output
parameter_set.config['globalOptions']['inputDir']       = '/data/input'
parameter_set.config['globalOptions']['outputDir']      = '/data/output'
parameter_set.config['globalOptions']['cloneMap']       = join(r'global_05min/cloneMaps',DEST_CLIP)
parameter_set.config['globalOptions']['startTime']      = START_TIME
parameter_set.config['globalOptions']['endTime']        = END_TIME
parameter_set.config['meteoOptions']['precipitationNC'] = join(r'global_05min/meteo',PRECIP_NAME)
parameter_set.config['meteoOptions']['temperatureNC']   = join(r'global_05min/meteo',TEMP_NAME)
parameter_set.config['meteoOptions']['referenceETPotMethod'] = 'Hamon'
parameter_set.config['meteoOptions']['refETPotFileNC']  = 'None'
parameter_set.config['routingOptions']['routingMethod']             ='accuTravelTime' #not kinematicWave
parameter_set.config['routingOptions']['dynamicFloodPlain']         ='False' #True doesn't work
parameter_set.config['meteoOptions']['precipitationVariableName']   = 'pr'
parameter_set.config['meteoOptions']['temperatureVariableName']     = 'tas'  #temperature at surface
parameter_set.config['meteoOptions']['referenceEPotVariableName']   ='evspsblpot'
parameter_set.config['meteoOptions']['temperatureConstant']         = '-273.15' #check if this is necessary
parameter_set.config['meteoDownscalingOptions']['downscalePrecipitation']='True'

if ADJUST_LANDCOVER!=False:
    if ADJUST_LANDCOVER=='GRASSLANDS':
        FRAC_PREFIX = 'grassf'
    elif ADJUST_LANDCOVER=='PROPORTIONAL':
        FRAC_PREFIX = 'propf'
    else:
        print('Wrong input ADJUST_LANDCOVER')
    parameter_set.config['landSurfaceOptions']['noLandCoverFractionCorrection']='True'
    parameter_set.config['forestOptions']['fracVegCover']   ='glaciers_landcov/'+FRAC_PREFIX+'_tall.map'
    parameter_set.config['grasslandOptions']['fracVegCover']='glaciers_landcov/'+FRAC_PREFIX+'_short.map'
    parameter_set.config['irrPaddyOptions']['fracVegCover'] ='glaciers_landcov/'+FRAC_PREFIX+'_pad.map'
    parameter_set.config['irrNonPaddyOptions']['fracVegCover']='glaciers_landcov/'+FRAC_PREFIX+'_nonpad.map'
elif ADJUST_LANDCOVER==False:
    parameter_set.config['landSurfaceOptions']['noLandCoverFractionCorrection']='False'
    parameter_set.config['forestOptions']['fracVegCover']   ='global_05min/landSurface/landCover/naturalTall/vegf_tall.map'
    parameter_set.config['grasslandOptions']['fracVegCover']='global_05min/landSurface/landCover/naturalShort/vegf_short.map'
    parameter_set.config['irrPaddyOptions']['fracVegCover'] ='global_05min/landSurface/landCover/irrPaddy/fractionPaddy.map'
    parameter_set.config['irrNonPaddyOptions']['fracVegCover']='global_05min/landSurface/landCover/irrNonPaddy/fractionNonPaddy.map'
else:
    print('Wrong input ADJUST_LANDCOVER')

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
    SUFFIX='test'+BASIN_NAME+str(time.clock())
else:
    SUFFIX='_'+RUN_NAME
OUT_DIR  =join(RUN_DIR,'output_maps','output'+SUFFIX)


print('call docker')
if CLUSTER == True:
    from grpc4bmi.bmi_client_singularity import BmiClientSingularity
    pcrg = BmiClientSingularity(image='ewatercycle-pcrg-grpc4bmi.sif',
                            input_dir=RUN_DIR,
                          output_dir=OUT_DIR)
else:
    from grpc4bmi.bmi_client_docker import BmiClientDocker
    pcrg = BmiClientDocker(image='ewatercycle/pcrg-grpc4bmi:setters', image_port=55555,
                    input_dir=RUN_DIR,
                    output_dir=OUT_DIR)



print('input: '+pcrg.input_dir)
print('output: '+pcrg.output_dir)
# initialize
start_time = time.time()
time.sleep(18)
print ('Initialize...')
pcrg.initialize(glob_ini_new)
print (str((time.time()-start_time )/60)+'minutes of initialize')
print ('\a')
# Usually 3 minutes for the Rhone


#%%Set up variables
T_START  = pcrg.get_start_time()
T_END    = pcrg.get_end_time()
T_STEP   = pcrg.get_time_step()
print('T_START = '+str(day2date(T_START)))

#Adust GloGEM timerange if necessary
if (day2date(T_START)!=nc.time.to_index()[0])&(COUPLE_GLOGEM==True):
    print('Ini-startdate and GloGEM startdate do not match!')
    print('Adjust GloGEM timespan...')
    nc = nc.sel(time=slice(day2date(T_START),day2date(T_END)))
    if (day2date(T_START)!=nc.time.to_index()[0])&(COUPLE_GLOGEM==True):
        print('Ini_startdate and GloGEM still do not match!')
        sys.exit()

#
LATSIZE,LONSIZE = pcrg.get_grid_shape(1)
if (LATSIZE!=nc.dims['lat']) or LONSIZE!=nc.dims['lon']:
    print ('Shapes GloGEM & PCRG do not match !')
    print('PCRG lat,lon dimensions: '+str(LATSIZE)+' , '+str(LONSIZE))
    print('GloGEM lat,lon dimensions: '+str(nc.dims['lat'])+ ' , '+str(nc.dims['lon']))
    if BASIN_NAME =='AMAZON':
        print('AMAZON is an exception, reduce GloGEM file')
        nc = nc.where(nc.lat<6,drop=True)
        isglac =np.any(nc.R>0,axis=0)
        N_NLAT -= 1
        if (LATSIZE!=nc.dims['lat']) or LONSIZE!=nc.dims['lon']:
            print('Dimensions still not right, exit')
            sys.exit()
    else:
        sys.exit()


#Nodes are middle of gridcells
LON_NODES       = pcrg.get_grid_y(1)
LAT_NODES       = pcrg.get_grid_x(1)

# Get hydrograph
# HG_LONLAT       = nc.hg.station_coords #Beaucaire
HG_LONLAT       = nc.station_coords1
HG_IDX          =  (find_nearest(LAT_NODES,HG_LONLAT[1]),
                      find_nearest(LON_NODES,HG_LONLAT[0]))


if BASIN_NAME =='THJORSA': #Thjorsa needs to be fixed manually
    HG_IDX = (9,18)
HG_IDX_flat   = np.ravel_multi_index(HG_IDX,(LATSIZE,LONSIZE))
print('HG_IDX='+str(HG_IDX))

Q_station      = []
Q_date          = []

# print ('Script before pcrg run takes',time.time()-timer1,'seconds')
#Timer1 = 206 seconds = 3.5 minutes (of which 3 for initialize)

SWE = []

#%% Start Run
i=0
start_time = time.time()
while pcrg.get_current_time()!=pcrg.get_end_time(): #Testrun=True gives another end condition
    full_loop_timer = time.time()
    if COUPLE_GLOGEM:
        #Coupling is done by adding the glogem runoff to the channel_storage
        chan_stor    = pcrg.get_value('channel_storage') #Gives flat array
        chan_stor    = np.reshape(chan_stor,(LATSIZE,LONSIZE))

        # chan_stor[isglac]=nc.R.isel(time=i).data[isglac]

        if nc.lat.data[0]>nc.lat.data[-1]: #(if Southern Hemisphere)
            chan_stor  += nc.R.data[i,::-1,:]
        else:
            chan_stor    += nc.R.data[i,:,:]

        chan_stor    = chan_stor.flatten()


        ####
        pcrg.set_value('channel_storage',chan_stor)
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
        discharge_reshaped =np.reshape(pcrg.get_value('discharge'),(N_NLAT,N_LON))
        localmax = np.nanmax(discharge_reshaped[HG_IDX[0]-1:HG_IDX[0]+2,
                                         HG_IDX[1]-1:HG_IDX[1]+2])

        Q_station.append(localmax)

    else:
        Q_station.append(pcrg.get_value_at_indices('discharge', HG_IDX_flat)[0])
    Q_date.append(day2date(pcrg.get_current_time()))
    # snowsum.append(np.sum(pcrg.get_value_at_indices('snow_melt',np.flatnonzero)))

    #Part to save additional variables
    #Variables: 'snow_water_equivalent', 'snow_melt'
    # for var in variables:
    #     vals    = pcrg.get_value(var) #Gives flat array
    #     vals   =np.reshape(vals,(LATSIZE,LONSIZE))
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
    # Z = np.reshape(ma.masked_where(vals == np.nan, vals), (LATSIZE,LONSIZE))
    # ax1.set_title(SUFFIX+'\n'+variable + '[' + unit + '], t='+str(day2date(pcrg.get_current_time())))
    # im = ax1.pcolormesh(LON_NODES,LAT_NODES,Z)
    # cbar = plt.colorbar(im,ax=ax1)
    # cbar.set_label('Discharge [m3/s]')

    i+=1
    # print ('Full loop timer:', time.time()-full_loop_timer)
    if (TEST_RUN==True)&(i==30):
        break




print(str((time.time()-start_time )/60)+'minutes of run')
print('Maximum daily discharge: '+str(round(np.max(Q_station),2)))


 #SAve hydrographs to csv
if (TEST_RUN==False):
    hg_name     = RUN_NAME+'_hg.txt' #
    hg_dir      = join(RUN_DIR,'Hydrographs',hg_name)
    hg_time     = pd.date_range(day2date(pcrg.get_start_time()),day2date(pcrg.get_current_time()))
    hg_df       = pd.DataFrame(Q_station,index=hg_time,columns = ['hg'])
    hg_df.index.name = 'time'
    hg_df.to_csv(hg_dir)

    # return
