# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 11:23:39 2021

Scripts accompanying PCR-GLOBWB 2 and GloGEM coupling
Auxiliary script: Creation of the adjusted PCR-GLOBWB 2 landcover maps


This script takes the PCRGLOB grassland (short vegetation) landcover fraction map
and subtracts the RGI glacier cover fraction from it.

Files Needed:
    landmask_global_05min.map
    glob_glac_frac.tif
    vegf_short.map
    vegf_tall.map
    fractionNonPaddy.map
    fractionPaddy.map

Other landcovers are needed because the vegetation fractions are expressed as fractions
of the pristine landcover (forest+grassland)

Output:
    short vegetation global landcover minus global glacier landcover 


@author: Internet
"""

# Adjust fractions
import rasterio
import glob
import os
import numpy as np

#
frac_path = r'd:\Documents\Master vakken\Thesis\Data\PCRGLOB\Adjust_fracs'

#Reads all 6 files in alphabetical order (lazy coding)
frac_paths = glob.glob(os.path.join(frac_path,'inputs','*'))
[f_nonpadf,f_padf,f_glacf,landmaskf,f_shortf,f_tallf]=\
    [rasterio.open(frac_paths[i],'r') for i in range(6)]

#Add PCRASTER_VALLUESCALE to prevent errors later on
irr_profile=f_nonpadf.profile
irr_profile.update({'PCRASTER_VALUESCALE':'VS_SCALAR'})
irr_profile.update({'driver':'GTiff'})
pris_profile=f_shortf.profile
pris_profile.update({'PCRASTER_VALUESCALE':'VS_SCALAR'})
pris_profile.update({'driver':'GTiff'})


# mask = [landmask==255]
f_nonpad = f_nonpadf.read(1)
f_pad    = f_padf.read(1)
f_glac   = f_glacf.read(1)
landmask = landmaskf.read(1)
f_short  = f_shortf.read(1)
f_tall   = f_tallf.read(1)

# f_glac = np.ones_like(f_glac)*0.9

#f_pad = -3.4028e38 for nodata
mask = f_pad<0
for raster in f_nonpad,f_pad,f_glac,f_short,f_tall:
    raster[mask]=np.nan

#if Sometimes paddy+non_paddy>1, iterate to make sure they're 1
for i in 0,1:
    f_irrigated             = f_pad+f_nonpad
    f_pad[f_irrigated>1]    =f_pad[f_irrigated>1]/f_irrigated[f_irrigated>1]
    f_nonpad[f_irrigated>1]   =f_nonpad[f_irrigated>1]/f_irrigated[f_irrigated>1]


#Take only from grassland:
f_pristine = 1-(f_pad+f_nonpad)
#Convert to absolute fraction, subtract f_glac and convert back to fraction of pristine landcover
f_short_new = np.where(f_pristine>0,
                       ((f_short*f_pristine)-f_glac)/f_pristine,
                       f_short)


print (np.nansum(f_short-f_short_new),np.nansum(f_glac)) # Shouldn't be too far apartt

out_path = os.path.join(frac_path,'PCRaster_test2')
filename ='vegf_short_minglac2'

# with rasterio.open(os.path.join(out_path,filename+'.map'),'w',**pris_profile) as zf_short:
#     zf_short.write_band(1,f_short_new)

#Save array as .txt and convert to .map using vegf_short as a clone
np.savetxt(os.path.join(out_path,filename+'.txt'),f_short_new)

import subprocess
os.chdir(out_path)
command = 'asc2map --clone vegf_short.map '+filename+'.txt '+filename+'.map'
print (command)
subprocess.Popen(command,shell=True)

"""
Profile should look like this
{'driver': 'PCRaster', 'dtype': 'float32', 'nodata': -3.4028234663852886e+38, 'width': 4320, 'height': 2160, 'count': 1, 'crs': None, 'transform': Affine(0.083333, 0.0, -180.0,
       0.0, -0.083333, 90.0), 'tiled': False}
"""


#%%
# #### WRITE WITH RIOXARRAY if above doesn-t work
# import xarray as xr
# import rioxarray as riox
# xr_file = xr.open_rasterio(frac_paths[0])
# riox_file  = riox.open_rasterio(frac_paths[0])
# riox_file.rio.to_raster(os.path.join(out_path,'test.map'))

# f_shortn = f_shortf.where(f_shortf==None,f_short)

# from os.path import join
# outpath = frac_path = r'/media/sf_Thesis/Data/PCRGLOB/Adjust_fracs/evenredig'
# f_shortn.rio.to_raster(join(outpath,'f_short_test.map'))
