#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 10:00:50 2020

@author: pau
"""
# from ewatercycle.parametersetdb import build_from_urls
# import ewatercycle
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
from matplotlib.colors import ListedColormap, LinearSegmentedColormap



# run_dir     = r'/media/pau/FREECOM/Thesis/Data/PCRGLOB/Global'# 
run_dir = r'd:\Documents\Master vakken\Thesis\Code'
os.chdir(run_dir)
# nc_dir = r'd:\Documents\Master vakken\Thesis\Data\NC_files'
nc_dir =r'd:\Documents\Master vakken\Thesis\Code\Files\glaciers_nc'

def months_slice(year,hemisphere,seasonal):
    if seasonal ==True:
        if hemisphere == 'North':
            frommonth = '03'
            tomonth   = '10'
            return slice(str(year)+'-'+frommonth,str(year)+'-'+tomonth)
        elif hemisphere == 'South':
            frommonth='09'
            tomonth  ='04'
            return slice(str(year-1)+'-'+frommonth,str(year)+'-'+tomonth)
    elif seasonal ==False:
        if hemisphere =='North':
            return slice(str(year-1)+'-10',str(year)+'-09')
        if hemisphere =='South':
            return slice(str(year-1)+'-04',str(year)+'-03')

def make_space_above(axes, topmargin=1):
    """ increase figure size to make topmargin (in inches) space for 
        titles, without changing the axes sizes"""
    fig = axes.flatten()[0].figure
    s = fig.subplotpars
    w, h = fig.get_size_inches()

    figh = h - (1-s.top)*h  + topmargin
    fig.subplots_adjust(bottom=s.bottom*h/figh, top=1-topmargin/figh)
    fig.set_figheight(figh)
    
def load_hg2(GRDC_no, full_daterange):
    hg2_path = join(run_dir,'Files','Q_obs',GRDC_no+'_Q_Day.Cmd.txt')
    hg = pd.read_csv(hg2_path,
                                 delimiter=';',
                                 skiprows=36,
                                 index_col=0,
                                 parse_dates=True,
                                 skipinitialspace=True,
                                 usecols=[0,2])
    hg = hg.rename(columns  = {'Value':'hg'})
    hg_full = pd.DataFrame(data=hg,index=full_date_range)
    hg_full.index = hg_full.index.rename('time')
    return hg_full.hg.to_xarray()

def FD_curve(data):
    y = data.sort_values(ascending=False).values
    x= np.cumsum(y)/np.nansum(y) *100
    return x,y

def NSEalt(model,benchmark,observation):
    return 1-np.sum((model-observation)**2)/np.sum((observation-benchmark)**2)

def RMSE(a,b):
    return np.sqrt(((a-b)**2).mean())

def BE(mod,obs,benchmark):
    return 1-(np.nansum((mod-obs)**2)/np.nansum((obs-benchmark)**2))
#%% INPUTS
files               = join(run_dir,'Files')
figures             =join(run_dir,'Figures')

basin_info = pd.read_csv(join(files,'basin_info_45min_old.csv'),index_col = 0)
# Basin_names = basin_info[basin_info['suitable']=='y'].index


# Basin_name          = 'INDUS'
# Basin_names = ['ALSEK','CLUTHA','DRAMSELV','FRASER','RHONE']
# Basin_names = ['CLUTHA','DANUBE','RHINE'] #station2
# Basin_names  = ['AMAZON','IRRAWADDY','MACKENZIE','OB','YUKON']
# Basin_names = ['ALSEK', 'CLUTHA', 'COLUMBIA', 'COPPER', 'DANUBE', 'DRAMSELV',
#         'FRASER', 'GLOMA', 'JOEKULSA', 
#         'KUSKOKWIM', 'LULE',  'NASS', 'NEGRO', 'OELFUSA',
#         'RHINE', 'RHONE', 'SANTA_CRUZ', 'SKAGIT', 'SKEENA', 'STIKINE',
#         'SUSITNA', 'TAKU', 'THJORSA', 'NELSON'] #minus large basins, Indus , KKalixaelven
Basin_names = ['AMAZON','IRRAWADDY','MACKENZIE',
                'OB','YUKON','ALSEK', 'CLUTHA', 'COLUMBIA', 'COPPER', 'DANUBE', 'DRAMSELV',
        'FRASER', 'GLOMA',  
        'KUSKOKWIM', 'NASS', 'NEGRO', 'OELFUSA',
        'RHINE', 'RHONE',  'SKAGIT', 'SKEENA', 'STIKINE','SUSITNA',
        'TAKU', 'THJORSA'] #minus Indus, Kalixaelven, Nelson, Joekulsa, Santa Cruz, Lule
# Basin_names = ['COPPER','ALSEK','COLUMBIA','DANUBE','RHINE','RHONE','OELFUSA']
# Basin_names = ['JOEKULSA','NELSON','SANTA_CRUZ','LULE']
# Basin_names = ['ALSEK','AMAZON','CLUTHA','COLUMBIA','NEGRO','SANTA_CRUZ','SUSITNA','THJORSA']
# Basin_names = ['SANTA_CRUZ']
model_setups        = ['0','1','2']
model_names         ={'0':'Modelled (Benchmark)',
                      '1':'Modelled (Bare)',
                      '2':'Modelled (Coupled)'}


seasonal            = False  #False for whole year, true for only summer
main_plot           =['s0','s1','s2']#,'s4']
only_obsyears       =True
save_figs           = False
new_seasonal        =False
calendar_day        =False


Model_name          ='PCRG'
GG_GCM              ='HadGEM2-ES'
GG_rcp              ='rcp26'
Fromyear            =2001
Untilyear           =2012
# station             =1

OF_list = []
normdiflist=[]
FD_list = []
NRD_list     = []
MBE_list = []
NRD2_list = []
NRD3_list = []
HG_list = []
glacier_sum_list =[]
Qobs_list = []

#%%
for Basin_name in Basin_names:
    if Basin_name in ['CLUTHA','COLUMBIA','RHINE','SUSITNA','DANUBE']:
        station_loop = [2]
    else: station_loop = [1]
    for station in station_loop:
        # Load hydrographs
        print (Basin_name)
        if Basin_name in ['CLUTHA','NEGRO','AMAZON','SANTA_CRUZ']:
            hemisphere = 'South'
            daterange = pd.date_range(str(Fromyear-1)+'-04-01',
                                     str(Untilyear)+'-03-31')
        else: 
            hemisphere = 'North'
            daterange = pd.date_range(str(Fromyear-1)+'-10-01',
                                      str(Untilyear)+'-09-30')
        # if Basin_name in ['CLUTHA','DANUBE','RHINE']:
        #     station = 2
        
        nc_path = join(nc_dir,'_'.join([Basin_name,
                                    GG_GCM,
                                    GG_rcp,
                                    '2000',
                                    '2016',
                                    'R.nc']))
        nc_obs          =xr.open_dataset(nc_path).sel(time=daterange)
        glacier_sum = nc_obs.R.sum(axis=(1,2))/(24*60*60)
        glacier_sum_list.append(glacier_sum)
        
        # if new_seasonal==True:
        #     nc_df = glacier_sum.to_dataframe()
        #     nc_df.pop('spatial_ref')
        #     nc_df.pop('height')
        #     R_months = nc_df.groupby(nc_df.index.month).sum()
        #     R_fracs  = R_months/R_months.sum()>0.001
        #     min_month = R_fracs.index[R_fracs.R].min()
        #     max_month = R_fracs.index[R_fracs.R].max()
        #     daterange = daterange[(daterange.month>min_month)&
        #                           (daterange.month<max_month)]
        #     nc_obs          =nc_obs.sel(time=daterange)
        #     glacier_sum = nc_obs.R.sum(axis=(1,2))/(24*60*60)
        
        if station ==1:
            Qobs       = nc_obs.hg
        elif station ==2:
            full_date_range= pd.to_datetime(nc_obs.time.data)
            Qobs = load_hg2(str(int(basin_info.loc[Basin_name,'grdc_no2'])),
                            full_date_range)
        Qobs = Qobs.where(Qobs!=-999,np.nan,drop=False)
        Qobs_list.append(Qobs)
        Q_time     = Qobs.where(~xr.ufuncs.isnan(Qobs),
                                drop=True).time
        Q_dates     = pd.to_datetime(Q_time.data)
        Q_years     = Q_dates.year.unique()
        
        
        hg = {}
        for setup in model_setups:
            hgdic = {'Model_name'          :Model_name,
               'Basin_name'          : Basin_name,
               'GG_GCM'              : GG_GCM ,
               'GG_rcp'              : GG_rcp,
               'Setting_no'          :setup,
               'Fromyear'           :'2001',
               'Untilyear'          :'2016',
               'station'            :'hg'+str(station)}
            hg_name = '_'.join(hgdic.values()) +'.txt'
            hg['s'+setup]=pd.read_csv(
                    join(run_dir,'Output','HG_2001_2016',
                    hg_name)
                    ,index_col=0,
                    parse_dates=True).loc[daterange]
        hg_years = daterange.year.unique()[1:]
        
        #make calendar day benchmark CDB
        if calendar_day == True:
            CDB_grouped = Qobs.groupby("time.dayofyear").mean().to_dataframe()['hg']
            dayofyear = Qobs['time.dayofyear'].to_dataframe()
            CDB_list = []
            t=time.time()
            for i in range(len(dayofyear)):
                for j in range(len(CDB_grouped)):
                    if dayofyear.dayofyear[i]==CDB_grouped.index[j]:
                        # print(dayofyear.dayofyear[i],CDB.index[j])
                        CDB_list.append(CDB_grouped.values[j])
            print (time.time()-t)
            CDB = pd.DataFrame({'CDB':CDB_list},index=dayofyear.index)
        
        #Calculate OF's
        OF = {}
        nanmask = ~xr.ufuncs.isnan(Qobs).data
        Obs = Qobs.data[nanmask]
        for setup in ['0','1','2']:
            Mod = hg['s'+setup].hg.values[nanmask]
            if calendar_day ==True:
                OF['BE'+setup]  =BE(Mod,Obs,CDB.CDB.values[nanmask])
            OF['RMSE'+setup]=he.rmse(Mod,Obs)
            OF['NSE'+setup] =he.nse(Mod,Obs)
            # OF['KGE'+setup] =he.kge(Mod,Obs)[0][0]
            OF['MB'+setup]  =Mod.sum() - Obs.sum()
        # OF['NSEalt']   = NSEalt(hg['s2'].hg.values[nanmask],
        #                         hg['s0'].hg.values[nanmask],Obs)
        OF['BE']   = BE(hg['s2'].hg.values[nanmask],
                                Obs,hg['s0'].hg.values[nanmask])
        
        
        OF_plot = ''
        for i in main_plot:
            F=''
            for key,item in OF.items():
                if i[1] in key:
                    F+= key+' = '+'{:8.3f}'.format(item)+'  '
            F+= '\n'
            OF_plot +=F
            # OF_plot+= key+'{8.3f}'.format(item)
        if 's1' in hg.keys():
            routed_R    = hg['s2']-hg['s1'] #routed glacier runoff
            glac_frac = routed_R.hg / hg['s2'].hg
            max_frac  = glac_frac.max() #.max() or .quantile()
            OF['GFmax'] = max_frac
            OF['GF99']  = glac_frac.quantile(0.99)
            OF['GFmean'] = glac_frac.mean()
        
        # Rsum = (nc_obs.R.sum(axis=(1,2))/(24*60*60))
        # frac = Rsum[0:-2]/Q2.hg[1:-1]
        
        # PCRG_glacR = Q0-Q1
        # routed_R   = Q2-Q1
        # Q02dif      =Q2-Q0
        
        # Plotting parameters
        
        if only_obsyears==True:
            years = Q_years
        else: 
            years=hg_years
        if years[0]==2000:
            years = years[1:]
        N = len(years)
        
    
        #% Regular plot
        #%
        multiplier = 1.7
        figheigth = multiplier*N #was 2.5
        topmargin = 1
        
        
        f1,(axes) = plt.subplots(N,1,figsize=(10,figheigth),sharey=True)
        f1.subplots_adjust(top=0.8)
        plt.subplots_adjust(hspace=0.22) #was 0.21
        # axes[0].set_title(Basin_name)
        f1.suptitle(Basin_name[0]+Basin_name[1:].lower(),
                    size='xx-large',y=(figheigth-0.15*topmargin)/figheigth)
        colors = ['tab:orange','tab:red','tab:green']
        # f1.text(0.1,(figheigth-0.8*topmargin)/figheigth,OF_plot)
        for i in range(N):
            ax = axes[i]
            months = months_slice(years[i],hemisphere,seasonal)
            Qo = Qobs.sel(time=months)

            c=0
            zorders  =[2,1,2]
            for key,item in hg.items():
                if key in main_plot:
                    Qm = item[months]
                    ax.plot(Qm.index,Qm,
                            linewidth=1.5,
                            label=model_names[key[1]],
                            color=colors[c],
                            zorder = zorders[c])
                    c+=1
            ax.plot(Qo.time,Qo,linewidth=1.5,alpha=1,label='Observed',
                    color='tab:blue')
                        # CDB1 = CDB[months]
            # ax.plot(CDB1.index,CDB1,label='CDB')
            # ax.legend() 
            # ax.set_ylim(0,np.nanpercentile(Qobs.data,99.9))
            ax.set_ylim(0,ax.get_yticks()[-1])
            xticks = ax.get_xticks()
            ax.tick_params(axis='both', which='major', labelsize=8)
            # ax.tick_params(axis='both', which='minor', labelsize=8)
            ax.set_xlim(xticks[0]-30,xticks[-1]+29)
            ax.set_ylabel(r'$Q_{basin}\/[m^3/s]$')
            ax.grid()
    
            twin = ax.twinx()
            glacR = glacier_sum.sel(time=months)
            twin.bar(glacR.time,glacR,alpha=0.4,
                      width = 1.5, label='GloGEM',color='tab:blue')
            twin.set_ylabel(r'$Q_{glacier}\/ [m^3/s]$')
            # twin.set_ylim(0,twin.get_ylim()[1]*2.5)
            twin.set_ylim(0,glacier_sum.data.max()*1.5)
            # twin.invert_yaxis()
            
            axticks = ax.get_yticks()
            twinticks = twin.get_yticks()
            while len(twinticks)<len(axticks):
                twinticks = np.append(twinticks,twinticks[-1]+twinticks[1])
            twin.set_ylim(0,twinticks[len(axticks)-1])
            twin.invert_yaxis()
            
            if i==0:
                lines, labels = ax.get_legend_handles_labels()
                lines2, labels2 = twin.get_legend_handles_labels()
                twin.legend(
                    lines + lines2, labels + labels2, loc=(0.8, 1.1), prop={"size": 8},
                )
                # ax.text(0,1.1,OF_plot,transform = ax.transAxes)
        # f1.legend(lines + lines2, labels + labels2, loc=(0.8, 
        #                         (figheigth-0.8*topmargin)/figheigth), prop={"size": 8})
        # f1.add_subplot(111,frameon=False)
        # plt.tick_params(labelcolor='none',top=False,bottom=False,left=False,right=False)
        # plt.ylabel('Common Y')
        make_space_above(axes,topmargin=topmargin)
        # if save_figs==True:
        f1.savefig(join(figures,'Main_plots_xticks',
                        '_'.join([Basin_name,'main'])+'small_zorder.png'),format='png',
                    bbox_inches = 'tight')
        
        f1.savefig(join(figures,'Main_plots_plusbare_svg',
                        '_'.join([Basin_name,'main'])+'small_zorder.svg'),format='svg',
                    bbox_inches = 'tight')
        # # # continue
        # # #
        
        # # ## Error plots
        # # f1,(axes) = plt.subplots(N,1,figsize=(10,figheigth),sharey=True)
        # # f1.subplots_adjust(top=0.8)
        # # f1.suptitle(Basin_name.title()+' - Error',
        # #             size='xx-large',y=(figheigth-0.15*topmargin)/figheigth)
        # # for i in range(N):
        # #     ax = axes[i]
        # #     months = months_slice(years[i],hemisphere,seasonal)
        # #     Qo = Qobs.sel(time=months)
        # #     for key,item in hg.items():
        # #         if key in main_plot:
        # #             Qm = item[months]
        # #             ax.plot(Qm.index,abs(Qm.hg-Qo.data),linewidth=1.5,label='Setup '+key[1]+' - obs')
    
        # #     ax.set_ylabel('Basin runoff [m3/s]')
        # #     ax.grid()
            
        # #     if i==0:
        # #         lines, labels = ax.get_legend_handles_labels()
        # #         ax.legend(
        # #             lines , labels, loc=(0.8, 1.1), prop={"size": 8},
        # #         )
        # #         ax.text(0,1.1,OF_plot,transform = ax.transAxes)
        # # make_space_above(axes,topmargin=topmargin)
        # # if save_figs==True:
        # #     f1.savefig(join(run_dir,'Figures','Error_plots',
        # #                 '_'.join([Basin_name,'error',str(station)])+'.png'),
        # #             bbox_inches = 'tight')
        
           
        
        # # #%% Glacier fraction of runoff
        # # f1,ax1 = plt.subplots(figsize=(10,6))
        # # # ax1.plot(Q2.index[1:-1],frac,linewidth=1,alpha=1,label='Observed')
        # # ax1.plot(Q2.index,Q2,label='Coupled model')
        # # ax1.plot(Q2.index,Rsum,label='Glacier output')
        # # plt.grid()
        # # plt.legend()
        # # plt.title('Rhone glacier output vs coupled PCRGlob')
        # # plt.ylabel('Discharge [m3/s]')
        
        # #%Difference between set1 and coupled
        # # if ('s1' in hg.keys())&('s2' in hg.keys()):
        # #     f1,ax1 = plt.subplots(figsize=(10,6))
        # #     ax1.plot(routed_R.index,routed_R,label='setup2 - setup1')
        # #     ax1.axhline(0,color='k',linestyle='--')
        # #     plt.grid()
        # #     plt.legend()
        # #     # plt.ylim(-750,750)
        # #     plt.title(Basin_name.title()+ ' - Routed glacier discharge approximation')
        # #     plt.ylabel('Discharge [m3/s]')
        
        # ##% Difference between Default and 1 
        # if ('s0' in hg.keys())&('s1' in hg.keys())&('s2' in hg.keys()):
        #     PCRG_glacR = hg['s1']-hg['s0']
        #     f1,ax1 = plt.subplots(figsize=(10,6))
        #     ax1.plot(PCRG_glacR.index,PCRG_glacR,label=model_names['1']+ ' - '+model_names['0'])
        #     ax1.plot(routed_R.index,routed_R,label=model_names['2']+ ' - '+model_names['1'])
        #     ax1.axhline(0,color='k',linestyle='--')
        #     plt.grid()
        #     plt.legend()
        #     # plt.ylim(-750,750)
        #     plt.title(Basin_name.title()+' - PCRGLOB vs GloGEM glacier runoff')
        #     plt.ylabel('Discharge [m3/s]')
        #     if save_figs==True:
        #         plt.savefig(join(run_dir,'Figures','glacialR_plots',
        #                 '_'.join([Basin_name,'glacialR',str(station)])+'.png'),
        #                 bbox_inches = 'tight')
            
        
        # if ('s0' in hg.keys())&('s2' in hg.keys()):
        #     f1,(axes) = plt.subplots(N,1,figsize=(10,figheigth),sharey=True)
        #     f1.subplots_adjust(top=0.8)
        #     f1.suptitle(Basin_name.title()+' - Benchmark comparison',
        #                 size='xx-large',y=(figheigth-0.15*topmargin)/figheigth)
        #     Q02dif = hg['s2']-hg['s0']
        #     for i in range(N):
        #         ax = axes[i]
        #         months = months_slice(years[i],hemisphere,seasonal)
        #         Qdif = Q02dif[months]
        #         ax.plot(Qdif.index,Qdif,linewidth=1.5,label=model_names['2']+' - '+model_names['0'])    
        #         ax.axhline(0,color='k',linestyle='--')
    
        #         ax.set_ylabel('dQ [m3/s]')
        #         ax.grid()
                
        #         if i==0:
        #             lines, labels = ax.get_legend_handles_labels()
        #             ax.legend(
        #                 lines , labels, loc=(0.8, 1.1), prop={"size": 8},
        #             )
        #             ax.text(0,1.1,OF_plot,transform = ax.transAxes)
        #     make_space_above(axes,topmargin=topmargin)
        #     if save_figs==True:
        #         f1.savefig(join(run_dir,'Figures','Seasonality_plots',
        #                     '_'.join([Basin_name,'Seasonality',str(station)])+'.png'),
        #                 bbox_inches = 'tight')
        
        if ('s0' in hg.keys())&('s2' in hg.keys()):
            Q02dif = hg['s2']-hg['s0']
            Qnormdif = Q02dif/Q02dif.quantile(0.99)
            # Qnormdif = Q02dif/Q02dif.max()
            normdiflist.append(Qnormdif.rename(columns={'hg':Basin_name}))

            # f1,(axes) = plt.subplots(N,1,figsize=(10,figheigth),sharey=True)
            # f1.subplots_adjust(top=0.8)
            # f1.suptitle(Basin_name.title()+' - Normalized difference by maximum',
            #             size='xx-large',y=(figheigth-0.15*topmargin)/figheigth)

            # for i in range(N):
            #     ax = axes[i]
            #     months = months_slice(years[i],hemisphere,seasonal)
            #     Qdif = Qnormdif[months]
            #     ax.plot(Qdif.index,Qdif,linewidth=1.5,label='dQ/max(dQ)')    
            #     ax.axhline(0,color='k',linestyle='--')
            #     ax.set_ylabel('[-]')
            #     ax.grid()
                
            #     if i==0:
            #         lines, labels = ax.get_legend_handles_labels()
            #         ax.legend(
            #             lines , labels, loc=(0.8, 1.1), prop={"size": 8},
            #         )
            #         ax.text(0,1.1,OF_plot,transform = ax.transAxes)
            # make_space_above(axes,topmargin=topmargin)
            # if save_figs==True:
            #     f1.savefig(join(run_dir,'Figures','normdif',
            #                 '_'.join([Basin_name,'normdif',str(station)])+'.png'),
            #             bbox_inches = 'tight')
        
        if ('s0' in hg.keys())&('s2' in hg.keys()):
            FD_bm = FD_curve(hg['s0'].hg)
            FD_cpl= FD_curve(hg['s2'].hg)
            FD_obs=FD_curve(Qobs.to_dataframe().hg)
            # f1,ax1=plt.subplots(figsize=(8,5))
            # ax1.plot(*FD_obs,label='GRDC')
            # ax1.plot(*FD_bm,label=model_names['0'])
            # ax1.plot(*FD_cpl,label=model_names['2'])
            # ax1.grid()
            # ax1.legend()
            # ax1.set_title(Basin_name+' - flow-duration curves')
            # ax1.set_ylabel(r'$Q_{basin} [m^3/s]$')
            # ax1.set_xlabel('Exceedance probability [%]')
            # # ax1.semilogy()
            # if save_figs==True:
            #     f1.savefig(join(run_dir,'Figures','logCDFs',
            #                 '_'.join([Basin_name,'logCDF',str(station)])+'.png'),
            #             bbox_inches = 'tight')
            p = np.linspace(0,100,1001)
            FD_bm_int = np.interp(p,FD_bm[0],FD_bm[1])
            FD_cpl_int = np.interp(p,FD_cpl[0],FD_cpl[1])
            FD_obs_int = np.interp(p,FD_obs[0],FD_obs[1])
            FD_diff   =1- FD_bm_int/FD_cpl_int
            # FD_diff = (FD_bm_int-FD_cpl_int)/np.max(np.abs(FD_bm_int-FD_cpl_int))
            FD_list.append(pd.DataFrame({Basin_name:FD_diff},index=p))
            OF['FDBE'] = BE(FD_cpl_int,FD_obs_int,FD_bm_int)

        
        if ('s0' in hg.keys())&('s1' in hg.keys())&('s2' in hg.keys()):
            DF = pd.concat([hg['s0'],hg['s1'],hg['s2']],axis=1)
            DF.columns =  ['s0','s1','s2']
            DF['Obs'] = Qobs
            # DF['dR01']=RMSE(DF['s0'],DF['Obs'])-RMSE(DF['s1'],DF['Obs'])
            # NRD = DF.groupby([DF.index.year,DF.index.month]).apply(
            #     lambda df,s0,s1,s2,obs:
            #         (RMSE(df[s0],df[obs])-RMSE(df[s2],df[obs]))/
            #          RMSE(df[s2],df[s1]),
            #             's0','s1','s2','Obs')
            NRD = DF.groupby([DF.index.year,DF.index.month]).apply(
                lambda df,s0,s1,s2,obs:
                    (RMSE(df[s0],df[obs])-RMSE(df[s2],df[obs]))/
                    np.where(df[s2].sum()>df[s0].sum(),RMSE(df[s2],df[s1]),RMSE(df[s0],df[s1])),
                        's0','s1','s2','Obs')
            NRD2 = DF.groupby([DF.index.year,DF.index.month]).apply(
                lambda df,s0,s1,s2,obs:
                    (np.sqrt(np.sum((df[s0]-df[obs])**2))-np.sqrt(np.sum((df[s2]-df[obs])**2)))/
                    np.sqrt(np.sum((np.max([df[s0],df[s2]],axis=0)-df[s1])**2)),
                        's0','s1','s2','Obs')
            # NRD3 = DF.groupby([DF.index.year,DF.index.month]).apply(
            #     lambda df,s0,s1,s2,obs:
            #         (np.sqrt(np.sum((df[s0]-df[obs])**2))-np.sqrt(np.sum((df[s2]-df[obs])**2)))/
            #         np.where(df[s2].sum()>df[s0].sum(),
            #                   np.sqrt(np.sum((df[s1]-df[s2])**2)),
            #                   np.sqrt(np.sum((df[s1]-df[s0])**2))),
                        # 's0','s1','s2','Obs') 
            NRD3 = DF.groupby([DF.index.year,DF.index.month]).apply(
                lambda df,s0,s1,s2,obs:
                    (np.sqrt(np.sum((df[s0]-df[obs])**2))-np.sqrt(np.sum((df[s2]-df[obs])**2)))/
                    (np.sqrt(np.sum((np.max([df[s0],df[s2]],axis=0)-df[obs])**2))-np.sqrt(np.sum((np.max([df[s0],df[s2]],axis=0)-df[s1])**2))),
                        's0','s1','s2','Obs')
            # NRD = DF.groupby([DF.index.year,DF.index.month]).apply(
            #     lambda df,s0,s1,s2,obs:
            #         ((np.sum((df[s0]-df[obs])**2))-(np.sum((df[s2]-df[obs])**2)))/
            #         (np.sum((np.max([df[s0],df[s2]],axis=0)-df[s1])**2)),
            #             's0','s1','s2','Obs')
            #Calculate monthly BE 
            MBE = DF.groupby([DF.index.year,DF.index.month]).apply(
                lambda df,s0,s1,s2,obs:
                    BE(df[s2],df[obs],df[s0]),
                        's0','s1','s2','Obs')
                
            # NRD.index = [pd.datetime.datetime(*NRD.index[i],1) for i in np.arange(len(NRD))]
            NRD.index = [pd.datetime.datetime(*NRD.index[i],1) for i in np.arange(len(NRD))]
            NRD[NRD>1] = 1
            NRD[NRD<-1]=-1
            
            MBE.index = [pd.datetime.datetime(*MBE.index[i],1) for i in np.arange(len(MBE))]
            
            NRD2.index = [pd.datetime.datetime(*NRD2.index[i],1) for i in np.arange(len(NRD2))]
            NRD2[NRD2>1] = 1
            NRD2[NRD2<-1]=-1
            
            NRD3.index = [pd.datetime.datetime(*NRD3.index[i],1) for i in np.arange(len(NRD3))]
            NRD3[NRD3>1] = 1
            NRD3[NRD3<-1]=-1
            
            mask = ((glacier_sum>0).resample(time='m').sum()>10).data
            #Mask all months that have on average less than 1 percent of annual glacier runoff
            glacsumDF = glacier_sum.to_dataframe()['R']
            pctmask = (glacsumDF.groupby([glacsumDF.index.month]).sum()/glacsumDF.sum())>0.01
            mask = [pctmask.loc[NRD.index[i].month] for i in range(len(NRD))]
            NRD=NRD.where(mask,np.nan)
            NRD_list.append(pd.DataFrame({Basin_name:NRD}))
            NRD2=NRD2.where(mask,np.nan)
            NRD2_list.append(pd.DataFrame({Basin_name:NRD2}))
            NRD3=NRD3.where(mask,np.nan)
            NRD3_list.append(pd.DataFrame({Basin_name:NRD3}))            
            MBE = MBE.where(mask,np.nan)
            MBE_list.append(pd.DataFrame({Basin_name:MBE}))
            
            # f1,(axes) = plt.subplots(N,1,figsize=(10,figheigth),sharey=True)
            # f1.subplots_adjust(top=0.8)
            # f1.suptitle(Basin_name.title()+' - Normalized difference by maximum',
            #             size='xx-large',y=(figheigth-0.15*topmargin)/figheigth)
            # for i in range(N):
            #     ax = axes[i]
            #     months = months_slice(years[i],hemisphere,seasonal)
            #     DF[months].plot(ax=ax,color=['orange','red','green','blue'])
            #     ax2=ax.twinx()
            #     NRD[months].plot(ax=ax2)
            #     ax.grid()
            #     ax2.axhline(0,color='black',linestyle='--')
            #     ax.set_ylim(0,ax.get_yticks()[-1])

            # make_space_above(axes,topmargin=0.4)

        OF_list.append(pd.DataFrame(OF,index=[Basin_name]))
        HG_list.append(hg)
        continue

#%%

#%%
OF_df = pd.concat(OF_list)
OF_df = OF_df.reindex(sorted(OF_df.columns,reverse=False),axis=1)
OF_df.to_csv(join(run_dir,'Output','OF_BE'+str(len(OF_list))+'.csv'))


OF_df = OF_df[OF_df.index!='SANTA_CRUZ']
# f1,ax1 = plt.subplots()
# # ax1.plot(OF_df['GFmean']*100,OF_df['NSEalt'],'ro',label='NSEalt')
# RMSEdif = 1-OF_df['RMSE2']/OF_df['RMSE0']
# ax1.plot(OF_df['GFmax'][RMSEdif>0]*100,RMSEdif[RMSEdif>0],'bo',
#          label=str(sum(RMSEdif>0))+'/'+str(len(OF_df))+' improved')
# ax1.plot(OF_df['GFmax'][RMSEdif<=0]*100,RMSEdif[RMSEdif<=0],'ro',
#     label=str(sum(RMSEdif<=0))+'/'+str(len(OF_df))+' deteriorated')
# ax1.grid()
# ax1.legend(loc='lower left')
# ax1.set_xlabel('max glacier runoff fraction [%]')
# ax1.set_title('1 - (RMSE_Coupled / RMSE_Benchmark)')
# ax1.axhline(0,color='k',linestyle='--')
# ax1.set_ylim(-1,1)

OF_sorted= OF_df.sort_values(by=['GF99'],axis=0,ascending=False)
# OF_sorted= OF_df.sort_values(by=['GFmax'],axis=0,ascending=False)
#%% Benchmark efficiency plots
OF_sorted.index = [OF_sorted.index[i].title() for i in range(len(OF_sorted.index))]

f1,ax=plt.subplots(1,3,figsize=(10,10),sharey=True)
f1.subplots_adjust(wspace=0.02)    
OF_sorted['TFBE'] = 1-(OF_sorted['MB2']/OF_sorted['MB0'])

OFS = ['BE','FDBE','TFBE']
for i in range(3):
    ax[i].plot(OF_sorted[OFS[i]],OF_sorted.index,marker='o',
             linestyle='None',color='orange',markeredgecolor='black')
    # for i in range(len(OF_sorted)):
    #     ax[i].arrow(x=i,y=OF_sorted['BE0'][i],
    #               dx = 0 ,dy =OF_sorted['BE2'][i]-OF_sorted['BE0'][i],
    #               width=0.2,length_includes_head=True,head_length=0.1)
    plt.xticks(rotation='horizontal')
    ax[i].axvline(0,color='black',linestyle='--')
    ax[i].set_xlim(right=1,left=-1.5)
    
    # x = range(len(OF_sorted))
    # y = np.array([1])
    # xx,yy = np.meshgrid(x,y)
    xmax = 1
    zz  = np.array([OF_sorted['GF99'].values]).transpose()
    cmap = plt.cm.Blues
    im = ax[i].imshow(zz,cmap=cmap,aspect='auto',
                    extent = (*ax[i].get_xlim(),len(OF_sorted)-0.5,-0.5))
    # cbar = f1.colorbar(im,ticks=np.linspace(0,1,6),fraction=0.05
    #                    ,label='Maximum glacier contribution')
    # ax[i].set_title('Flow-duration benchmark efficiency')
    ax[i].set_xlabel(OFS[i]+' [-]')
    ax[i].grid()
    ax[i].set_xticks([-1,-0.5,0,0.5,1])
    
    for j in range(len(OF_sorted[OFS[i]])):
        OF = OF_sorted[OFS[i]]
        if OF[j]<-1.5:
            # ax[i].plot(-1.5,j,color='red',marker='x')
            ax[i].annotate(str(round(OF[j],2)),(-1.45,j)
                           ,bbox=dict( fc="0.9",alpha=0.6))
    
    # ax[i].set_xticks(rotation=45)
    # plt.setp(ax[i].xaxis.get_majorticklabels(), rotation=45)

# ax1.set_xlabel('FDBE [-]')
# ax1.set_ylim(-15,1)
# plt.setp(axa.xaxis.get_majorticklabels(), rotation=45)


ax2 = f1.add_axes([0,0,0.95,1])
ax2.set_visible(False)
im = ax2.imshow(np.array([[0,1]]),cmap=cmap)
bar =plt.colorbar(im,fraction=0.02,label=r'$Q_{99}$'+' glacier contribution [-]')
# f1.savefig(join(figures,'Overallmetrics_99.svg'),format = 'svg',bbox_inches = 'tight')
# im = ax1.imshow(zz,cmap=palette,aspect='auto',vmin=0,vmax=1,
#                 extent = (*ax1.get_xlim(),len(NRD)-0.5,-0.5))

#%%NSE plot
OF_sorted.index = [OF_sorted.index[i].title() for i in range(len(OF_sorted.index))]

f1,ax=plt.subplots(1,2,figsize=(8,10),sharey=True)
f1.subplots_adjust(wspace=0.02)    

OFS = ['NSE0','BE0']
for i in range(2):
    ax[i].plot(OF_sorted[OFS[i]],OF_sorted.index,marker='o',
             linestyle='None',color='orange',markeredgecolor='black')
    # for i in range(len(OF_sorted)):
    #     ax[i].arrow(x=i,y=OF_sorted['BE0'][i],
    #               dx = 0 ,dy =OF_sorted['BE2'][i]-OF_sorted['BE0'][i],
    #               width=0.2,length_includes_head=True,head_length=0.1)
    plt.xticks(rotation='horizontal')
    ax[i].axvline(0,color='black',linestyle='--')
    ax[i].set_xlim(right=1,left=-1.5)
    
    # x = range(len(OF_sorted))
    # y = np.array([1])
    # xx,yy = np.meshgrid(x,y)
    xmax = 1
    zz  = np.array([OF_sorted['GF99'].values]).transpose()
    cmap = plt.cm.Blues
    im = ax[i].imshow(zz,cmap=cmap,aspect='auto',
                    extent = (*ax[i].get_xlim(),len(OF_sorted)-0.5,-0.5))
    # cbar = f1.colorbar(im,ticks=np.linspace(0,1,6),fraction=0.05
    #                    ,label='Maximum glacier contribution')
    # ax[i].set_title('Flow-duration benchmark efficiency')
    # ax[i].set_xlabel(OFS[i]+' [-]')
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
    
    # ax.set_xticks(rotation=45)
    # plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)

# ax1.set_xlabel('FDBE [-]')
# ax1.set_ylim(-15,1)
# plt.setp(axa.xaxis.get_majorticklabels(), rotation=45)


ax2 = f1.add_axes([0.03,0.0,0.95,1])
ax2.set_visible(False)
im = ax2.imshow(np.array([[0,1]]),cmap=cmap)
bar =plt.colorbar(im,fraction=0.03,label=r'$Q_{99}$'+' glacier contribution [-]')
# f1.savefig(join(figures,'NSE0_BBE0.svg'),format='svg',bbox_inches = 'tight')
#%%
normdif_stack = pd.concat(normdiflist,axis=1)
for b in ['SANTA_CRUZ','NEGRO','AMAZON','CLUTHA']: #NEGRO, AMAZON, SANTA CRUZ, CLUTHA
    if b in normdif_stack.keys():    
        normdif_stack[b]=normdif_stack[b].shift(periods=182,freq='D')

normdif_means = normdif_stack.groupby([normdif_stack.index.month, normdif_stack.index.day]).mean()
normdif_means.index=normdif_stack['2004'].index


# normdif_popped=pd.concat([normdif_means.pop(b) for b in [
#     'AMAZON','OB','NEGRO','IRRAWADDY']]
#                         ,axis=1)
# for b in ['AMAZON','OB','NEGRO','IRRAWADDY']:
# for b in ['AMAZON','OB','IRRAWADDY']:
for b in ['SANTA_CRUZ']:
    if b in normdif_means.keys():
        normdif_means.pop(b)

# normdif_means.columns = OF_df.GF99.values
# normdif_means = normdif_means.reindex(sorted(normdif_means.columns),axis=1)

import matplotlib.colors
# cmap = plt.cm.gnuplot2_r
# cmap = plt.cm.viridis_r
# cmap = plt.cm.plasma_r
# cmap   = plt.cm.PuBu

pubu = plt.get_cmap('PuBu',200)
cmap = ListedColormap(pubu(range(200))[40:])


# cmap = newcmp
# cmap   = plt.cm.Blues
vmin =0
norm = matplotlib.colors.Normalize(vmin=vmin, vmax=1)
# cmap = cmap(norm(OF_df.GF99.values))
# OF_sorted_ascending= OF_df.sort_values(by=['GF99'],axis=0,ascending=True)
# rcParams['axes.prop_cycle'] = cycler(color=cmap(norm(OF_sorted_ascending.GF99.values)))
# rcParams['axes.prop_cycle'] = cycler(color=cmap(norm(OF_df.GF99.values)))
rcParams['axes.prop_cycle'] = cycler(color=cmap(norm(OF_sorted.GF99.values[2:])))

xticks = normdif_means.plot(linewidth=1.2,alpha =0,legend=False).get_xticks()

f1,ax1 = plt.subplots(figsize=(12,8))
# im = ax1.plot(normdif_means.index,normdif_means.values,alpha=0)
im = ax1.plot(normdif_means.index,normdif_means.values[:,20])
#plusminplus ND : 4,22, 20SKEENA

# cbar = f1.colorbar(cmap)
quantile_alpha = 0
normdif_means.quantile(0.75,axis=1).plot(ax=ax1,color='red',
                                        linewidth=1.5,label='75',
                                        linestyle='--',alpha=quantile_alpha)
normdif_means.mean(axis=1).plot(ax=ax1,color='black',linewidth=2,
                                label='50',alpha=quantile_alpha)
normdif_means.quantile(0.25,axis=1).plot(ax=ax1,color='red',
                                          linewidth=1.5,label='25',
                                          linestyle='--',alpha=quantile_alpha)

ax1.set_ylabel('Coupled model - Benchmark \n Normalized difference  [-]')
# ax1.set_title('Comparison between Benchmark and Coupled')
# ax1.legend(loc='upper left',title='Quantiles')
ax1.axhline(0,color='k',linestyle='--')
monthcombis = ['January / July',
               'February / August',
               'March / September',
               'April / October',
               'May / November',
               'June / December',
               'July / January',
               'August / February',
               'September / March',
               'October / April',
               'November / May',
               'December / June']
# xticks = ax1.get_xticks()
plt.xticks(xticks[:-1],monthcombis,rotation = 45)
ax1.grid()
ax1.set_xlim(normdif_means.index[0],normdif_means.index[-1])
ax2 = f1.add_axes([0,0,0.93,1])
ax2.set_visible(False)
im = ax2.imshow(np.array([[vmin,1]]),cmap=cmap)
bar =plt.colorbar(im,fraction=0.015,label=r'$Q_{99}$'+' glacier contribution [-]')
ax1.set_ylim(-0.8125,0.9957)
# xlabels=ax1.get_xticklabels()
# ax1.set_xticks(xticks-5,)

# ax1.set_title('Thjorsa')
f1.savefig(join(figures,'ND_oneplots.svg'),format='svg',bbox_inches = 'tight')
# 
# bar.ax.set_yticklabels(['0','1','2','3','4','5'])
# plt.xticks(n)
# ax1.set_xlabel('Months')
# ax1.set_xticks()



#%%
FD_concat = pd.concat(FD_list,axis=1)
FD_means  = FD_concat.quantile(0.5,axis=1)
FD_75     = FD_concat.quantile(0.75,axis=1)
FD_25     = FD_concat.quantile(0.25,axis=1)

# pubu = plt.get_cmap('PuBu',200)
# cmap = ListedColormap(pubu(range(200))[40:])
# # cmap = newcmp
# # cmap   = plt.cm.Blues
# vmin =0
# norm = matplotlib.colors.Normalize(vmin=vmin, vmax=1)
xx,yy = np.meshgrid(np.arange(25),np.linspace(0,100,1001))
f1,ax1 = plt.subplots(figsize=(12,8))
# FD_concat.plot(ax=ax1,legend=False)
ax1.plot(yy,FD_concat)
FD_75.plot(ax=ax1,color='red',linewidth = 1.5,label='75',linestyle='--')
FD_means.plot(ax=ax1,color='black',linewidth=2,label='50')
FD_25.plot(ax=ax1,color='red',linewidth = 1.5,label='25',linestyle='--')
ax1.grid()
ax1.legend(title='Quantiles',loc='upper left')
ax1.set_ylabel('Flow-duration ratio [-]')
ax1.set_xlabel('Exceedance probability [%]')
ax2 = f1.add_axes([0,0,0.93,1])
ax2.set_visible(False)
im = ax2.imshow(np.array([[vmin,1]]),cmap=cmap)
bar =plt.colorbar(im,fraction=0.015,label=r'$Q_{99}$'+' glacier contribution [-]')
ax1.set_xlim(0,100)
# ax1.semilogy()
# ax1.set_ylim(-1,1)
# f1.savefig(join(figures,'FDR_vectorplot_99.svg'),format='svg',bbox_inches = 'tight')

#%%
NRD_concat = pd.concat(NRD2_list,axis=1)
# NRD_concat['AMAZON'].shift(periods=6,freq='MS').plot()
for b in ['NEGRO','SANTA_CRUZ','AMAZON','CLUTHA']:
    if b in NRD_concat.keys():
        # NRD_concat[b]=np.roll(NRD_concat[b].values,6)
        NRD_concat[b] = NRD_concat[b].shift(periods=6,freq='MS')
    
NRD_means = NRD_concat.groupby([NRD_concat.index.month]).mean()
# NRD_Q50   = NRD_concat.groupby([NRD_concat.index.month]).quantile(0.5)
NRD_means.to_csv(join(run_dir,'Output','RRD_means.csv'))


f1,ax1 = plt.subplots(figsize=(20,15))
NRD_means.plot(ax=ax1,alpha=0.8,legend=False)
# NRD_concat.plot(ax=ax1,alpha=0.5)
# NRD_concat.quantile(0.5,axis=1).plot(ax=ax1,color='black',label='Q50')
# NRD_concat.quantile(0.75,axis=1).plot(ax=ax1,color='red',label='Q75')
# NRD_concat.quantile(0.25,axis=1).plot(ax=ax1,color='red',label='Q25')


NRD_means.quantile(0.5,axis=1).plot(ax=ax1,color='black',label='Q50',linewidth=3)
NRD_means.quantile(0.75,axis=1).plot(ax=ax1,color='red',label='Q75',linewidth=3)
NRD_means.quantile(0.25,axis=1).plot(ax=ax1,color='red',label='Q25',linewidth=3)
ax1.legend(loc=(1,0))
ax1.axhline(0,color='black',linestyle='--')
ax1.grid()
ax1.set_ylabel('NRD [-]')
ax1.set_title('Mean normalized RMSE difference')
ax1.set_xlabel('Months')
# ax1.axvline(5,color='black',linestyle='--')
# ax1.axvline(9,color='black',linestyle='--')
plt.xticks(range(1,13),monthcombis)
# ax1.set_ylim(-2,1)
#%%
OF27_path = join(run_dir,'Output','OF_BE25.csv')
OF27 = pd.read_csv(OF27_path,index_col=0)
NRD_path = join(run_dir,'Output','RRD_means.csv')
NRD = pd.read_csv(NRD_path,index_col=0).transpose()
NRD = NRD[NRD.index!='SANTA_CRUZ']
# NRD = NRD[NRD.index!='RHONE']
NRD['GF99']=OF27['GF99']
NRD = NRD.sort_values(by='GF99',ascending=False)
NRD.index = [NRD.index[i].title() for i in range(len(NRD.index))]

NRD_summermean = NRD[list(range(5,10))].mean(axis=1)
a,b,r,p,std = stats.linregress(range(len(NRD_summermean)),NRD_summermean.values)
x = np.arange(26)
NRD_lin = a*x+b
#%%
f1,ax1 = plt.subplots(figsize=(10,6))
colors = ['tab:cyan','tab:orange','tab:green','tab:red','tab:purple']
# colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple']
rscores = []
NRD_stack =np.zeros(2*len(NRD))
GF99_stack = np.zeros(2*len(NRD))
for i in range(5,10):
    if (i==7)or (i==8):
        NRD_stack[(i-7)*len(NRD):(i-6)*len(NRD)]= NRD[i].values
        GF99_stack[(i-7)*len(NRD):(i-6)*len(NRD)] = NRD['GF99'].values
    NRD.plot(kind='scatter',x=i,y='GF99',ax=ax1,color=colors[i-5])
    a,b,r,p,std = stats.linregress(x=NRD[i].values,
                                    y=NRD['GF99'].values)
    x = np.linspace(-1,1,20)
    ax1.plot(x,a*x+b,label = 'Month '+str(i)+'\n r2='+str(round(r**2,2))+ '\n',
             color=colors[i-5])
    rscores.append(r**2)
a,b,r,p,std = stats.linregress(x=NRD_stack,
                                y=GF99_stack)
x = np.linspace(-1,1,20)
ax1.plot(x,a*x+b,label = 'July/August \n r2= '+str(round(r**2,2)),
         color='black',linestyle='--')
ax1.set_xlim(-1,1)
ax1.set_xlabel('NRD')
ax1.grid()
ax1.legend(loc=(1,0))
#%%
# Without cax and R2 
# f1,ax1 = plt.subplots(figsize=(6,10))
# for month in range(5,10):
#     ax1.plot(NRD[month],NRD.index,color=colors[month-5],
#               linestyle='-',marker='o',label=str(month))
    
# markers=['o','v','^','s','D']
# f1,ax1 = plt.subplots(figsize=(6,10))
# for month in range(5,10):
#     ax1.plot(NRD[month],NRD.index,color=colors[month-5],
#               linestyle='-',marker=markers[month-5],label=str(month),
#               markeredgecolor='black',markeredgewidth=0.5,markersize=7)
# plt.xticks(rotation = 'vertical')

# ax1.grid()
# ax1.set_xlim(right=1,left=-1)
# ax1.set_xlabel('BE [-]')
# # ax1.plot(NRD_lin,NRD.index)
# plt.legend(title='Month')
# xmax = 1
# zz  = np.array([NRD['GF99'].values]).transpose()
# palette = plt.cm.Blues
# norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
# im = ax1.imshow(zz,cmap=palette,aspect='auto',vmin=0,vmax=1,
#                 extent = (*ax1.get_xlim(),len(NRD)-0.5,-0.5))
# cbar = f1.colorbar(im,ticks=np.linspace(0,1,6),fraction=0.05
#                     ,label='Maximum glacier contribution [-]')
# # cbar.ax.set_yticklabels(['0','0.2','0.4','0.6','0.8','1'])
# ax1.set_title('Monthly BE score')
# ax1.axvline(0,linestyle='--',color='black',alpha=0.6)


# With cax and R2
# f1,ax1 = plt.subplots(figsize=(6,10))
# for month in range(5,10):
#     ax1.plot(NRD[month],NRD.index,color=colors[month-5],
#               linestyle='-',marker='o',label=str(month))
markers=['o','v','^','s','D']
f1,ax1 = plt.subplots(figsize=(6,10))
alpha = [0,0,0,0,1]
for month in range(5,10):
    ax1.plot(NRD[month],NRD.index,color=colors[month-5],
              linestyle='-.',marker=markers[month-5],
              label=monthcombis[month-1]+'\n'+r'$R^2$ = '+str(round(rscores[month-5],2))+'\n',
              markeredgecolor='black',markeredgewidth=0.5,markersize=7,
              alpha=alpha[month-5])
plt.xticks(rotation = 0)

ax1.grid()
ax1.set_xlim(right=1,left=-1)
ax1.set_xlabel('RRD [-]')
# ax1.plot(NRD_lin,NRD.index)
plt.legend(loc = (1.01,0))
xmax = 1
zz  = np.array([NRD['GF99'].values]).transpose()
palette = plt.cm.Blues
import matplotlib.colors
norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
im = ax1.imshow(zz,cmap=palette,aspect='auto',vmin=0,vmax=1,
                extent = (*ax1.get_xlim(),len(NRD)-0.5,-0.5))
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
                   ,label=r'$Q_{99}$'+' glacier contribution [-]')
# cbar.ax.set_yticklabels(['0','0.2','0.4','0.6','0.8','1'])
# ax1.set_title('Monthly NRBE score')
ax1.axvline(0,linestyle='--',color='black',alpha=0.6)
# plt.xticks(rotation='horizontal')

f1.savefig(join(figures,'RRD_00001.svg'),format='svg',bbox_inches = 'tight')

# f1.savefig(join(figures,'NRBE_99.svg'),format='svg',bbox_inches = 'tight')
# f1.savefig(join(figures,'BE_99.svg'),format='svg',bbox_inches = 'tight')

# #%%
# xx,yy = np.meshgrid(NRD.keys()[5:10],NRD.GFmax)
# zz    = NRD[list(range(5,10))].values
# plt.pcolor(xx,yy,zz)

#%%
HG_dic = {Basin_names[i]:HG_list[i] for i in range(len(Basin_names))}
Qobs_dic = {Basin_names[i]:Qobs_list[i] for i in range(len(Basin_names))}
glacier_sum_dic = {Basin_names[i]:glacier_sum_list[i] for i in range(len(Basin_names))}
#%%
Basins = ['ALSEK','COLUMBIA','OELFUSA','MACKENZIE','RHINE','RHONE']
# Basins = ['JOEKULSA','SANTA_CRUZ','NELSON','LULE']
NB = len(Basins)

# year = 2005
for year in range(2001,2012):
    
    f1,axes = plt.subplots(NB,1,figsize=(10,1.9*NB),sharex=True)
    # f1.subplots_adjust(top=0.8)
    plt.subplots_adjust(hspace=0.31)
    # f1.suptitle(Basin_name[0]+Basin_name[1:].lower(),
    #             size='xx-large',y=(figheigth-0.15*topmargin)/figheigth)
    colors = ['tab:orange','tab:red','tab:green']
    for i in range(NB):
        hg = HG_dic[Basins[i]]
        Qobs = Qobs_dic[Basins[i]]
        glacier_sum = glacier_sum_dic[Basins[i]]
        
        ax = axes[i]
        months = months_slice(year,hemisphere,seasonal)
        Qo = Qobs.sel(time=months)
    
        c=0
        zorders = [2,1,2]
        for key,item in hg.items():
            if key in main_plot:
                Qm = item[months]
                ax.plot(Qm.index,Qm,
                        linewidth=1.5,
                        label=model_names[key[1]],
                        color=colors[c],
                        zorder = zorders[c])
                c+=1
        ax.plot(Qo.time,Qo,linewidth=1.5,alpha=1,label='Observed',
                color='tab:blue')
                    # CDB1 = CDB[months]
        # ax.plot(CDB1.index,CDB1,label='CDB')
        # ax.legend() 
        # ax.set_ylim(0,np.nanpercentile(Qobs.data,99.9))
        ax.set_ylim(0,ax.get_yticks()[-1])
        xticks = ax.get_xticks()
        ax.set_xlim(xticks[0]-30,xticks[-1]+29)
        ax.tick_params(axis='both', which='major', labelsize=8)
        # ax.tick_params(axis='both', which='minor', labelsize=8)
        
        ax.set_ylabel(r'$Q_{basin}\/[m^3/s]$')
        ax.grid()
        # gd = OF_sorted.loc[Basins[i].title(),'GF99']
        gd = basin_info.loc[Basins[i],'glac_degree']
        if i==0:
            ax.set_title(Basins[i].title()+' (Glaciation degree ='+str(round(gd,2))+'%)')
        else:
            ax.set_title(Basins[i].title()+' ('+str(round(gd,2))+'%)')
        twin = ax.twinx()
        glacR = glacier_sum.sel(time=months)
        twin.bar(glacR.time,glacR,alpha=0.4,
                  width = 1.5, label='GloGEM',color='tab:blue')
        twin.set_ylabel(r'$Q_{glacier}\/ [m^3/s]$')
        # twin.set_ylim(0,twin.get_ylim()[1]*2.5)
        twin.set_ylim(0,glacier_sum.data.max()*1.5)

        # twin.invert_yaxis()
        
        axticks = ax.get_yticks()
        twinticks = twin.get_yticks()
        while len(twinticks)<len(axticks):
            twinticks = np.append(twinticks,twinticks[-1]+twinticks[1])
        twin.set_ylim(0,twinticks[len(axticks)-1])
        twin.invert_yaxis()
        
        if i==0:
            lines, labels = ax.get_legend_handles_labels()
            lines2, labels2 = twin.get_legend_handles_labels()
            twin.legend(
                lines + lines2, labels + labels2, loc=(1.105, 0.1), prop={"size": 8}) #loc=(0.8,1.1)
           
        if i==NB-1:
            binitials = ''.join([B[0] for B in Basins ])
            # f1.savefig(join(figures,'Ensembleplots','Ensembleplot'+binitials+str(year)+'small.svg'),format='svg',
            # bbox_inches = 'tight')
    # f1.legend(lines + lines2, labels + labels2, loc=(0.8, 
    #                         (figheigth-0.8*topmargin)/figheigth), prop={"size": 8})
    # f1.add_subplot(111,frameon=False)
    # plt.tick_params(labelcolor='none',top=False,bottom=False,left=False,right=False)
    # plt.ylabel('Common Y')
    # make_space_above(axes,topmargin=topmargin)
    # binitials = ''.join([B[0] for B in Basins ])
    # f1.savefig(join(figures,'Ensembleplots','Ensembleplot'+binitials+str(year)+'small.svg'),format='svg',
    #             bbox_inches = 'tight')

#%%
# for i in range(1):#len(Basin_names)):
#     Basin_name = Basin_names[i]
#     hg = HG_dic[Basin_name]
#     Qobs = Qobs_dic[Basin_name]
#     glacier_sum = glacier_sum_dic[Basin_name]
#     figheigth = 2.5*N
#     topmargin = 1
    
    
#     f1,(axes) = plt.subplots(N,1,figsize=(20,figheigth),sharey=True)
#     f1.subplots_adjust(top=0.8)
#     plt.subplots_adjust(hspace=0.21)
#     # axes[0].set_title(Basin_name)
#     f1.suptitle(Basin_name[0]+Basin_name[1:].lower(),
#                 size='xx-large',y=(figheigth-0.15*topmargin)/figheigth)
#     colors = ['tab:orange','tab:red','tab:green']
#     # f1.text(0.1,(figheigth-0.8*topmargin)/figheigth,OF_plot)
#     for i in range(N):
#         ax = axes[i]
#         months = months_slice(years[i],hemisphere,seasonal)
#         Qo = Qobs.sel(time=months)
    
#         c=0
#         for key,item in hg.items():
#             if key in main_plot:
#                 Qm = item[months]
#                 ax.plot(Qm.index,Qm,
#                         linewidth=1.5,
#                         label=model_names[key[1]],
#                         color=colors[c])
#                 c+=1
#         ax.plot(Qo.time,Qo,linewidth=1.5,alpha=1,label='GRDC',
#                 color='tab:blue')
#                     # CDB1 = CDB[months]
#         # ax.plot(CDB1.index,CDB1,label='CDB')
#         # ax.legend() 
#         # ax.set_ylim(0,np.nanpercentile(Qobs.data,99.9))
#         ax.set_ylim(0,ax.get_yticks()[-1])
#         xticks = ax.get_xticks()
#         ax.set_xlim(xticks[0]-30,xticks[-1]+29)
#         ax.set_ylabel(r'$Q_{basin}\/[m^3/s]$')
#         ax.grid()
    
#         twin = ax.twinx()
#         glacR = glacier_sum.sel(time=months)
#         twin.bar(glacR.time,glacR,alpha=0.4,
#                   width = 1.5, label='GloGEM',color='tab:blue')
#         twin.set_ylabel(r'$Q_{glacier}\/ [m^3/s]$')
#         # twin.set_ylim(0,twin.get_ylim()[1]*2.5)
#         twin.set_ylim(0,glacier_sum.data.max()*1.5)
#         # twin.invert_yaxis()
        
#         axticks = ax.get_yticks()
#         twinticks = twin.get_yticks()
#         while len(twinticks)<len(axticks):
#             twinticks = np.append(twinticks,twinticks[-1]+twinticks[1])
#         twin.set_ylim(0,twinticks[len(axticks)-1])
#         twin.invert_yaxis()
        
#         if i==0:
#             lines, labels = ax.get_legend_handles_labels()
#             lines2, labels2 = twin.get_legend_handles_labels()
#             twin.legend(
#                 lines + lines2, labels + labels2, loc=(0.8, 1.1), prop={"size": 8},
#             )
#             # ax.text(0,1.1,OF_plot,transform = ax.transAxes)
#     # f1.legend(lines + lines2, labels + labels2, loc=(0.8, 
#     #                         (figheigth-0.8*topmargin)/figheigth), prop={"size": 8})
#     # f1.add_subplot(111,frameon=False)
#     # plt.tick_params(labelcolor='none',top=False,bottom=False,left=False,right=False)
#     # plt.ylabel('Common Y')
#     make_space_above(axes,topmargin=topmargin)
#     # if save_figs==True:
#     f1.savefig(join(figures,'Main_plots_plusbare_svg',
#                     '_'.join([Basin_name,'main'])+'.svg'),format='svg',
#                 bbox_inches = 'tight')

#%% Calculate percentages normdif_stack = pd.concat(normdiflist,axis=1)
normdif_stack = pd.concat(normdiflist,axis=1)
for b in ['SANTA_CRUZ','NEGRO','AMAZON','CLUTHA']: #NEGRO, AMAZON, SANTA CRUZ, CLUTHA
    if b in normdif_stack.keys():    
        normdif_stack[b]=normdif_stack[b].shift(periods=182,freq='D')

normdif_means = normdif_stack.groupby([normdif_stack.index.month, normdif_stack.index.day]).mean()
normdif_means.index=normdif_stack['2004'].index


import matplotlib.colors

pubu = plt.get_cmap('PuBu',200)
cmap = ListedColormap(pubu(range(200))[40:])


vmin =0
norm = matplotlib.colors.Normalize(vmin=vmin, vmax=1)
rcParams['axes.prop_cycle'] = cycler(color=cmap(norm(OF_df.GF99.values)))

xticks = normdif_means.plot(linewidth=1.2,alpha =0,legend=False).get_xticks()

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
# ax1.set_title('Comparison between Benchmark and Coupled')
ax1.legend(loc='upper left',title='Quantiles')
ax1.axhline(0,color='k',linestyle='--')
monthcombis = ['January / July',
               'February / August',
               'March / September',
               'April / October',
               'May / November',
               'June / December',
               'July / January',
               'August / February',
               'September / March',
               'October / April',
               'November / May',
               'December / June']
# plt.xticks(xticks[:-1],monthcombis,rotation = 45)
ax1.grid()
ax1.set_xlim(normdif_means.index[0],normdif_means.index[-1])
ax1.set_ylim(-0.8125,0.9957)
ax1.text(xticks[-1]-10,0.85,'a)',size=15)
xtickss = ax1.get_xticks()
xtickss = np.append(xtickss,xtickss[-1]+30)
summer ={}
for Basin_name in Basin_names:
    ratio = HG_dic[Basin_name]['s2']/HG_dic[Basin_name]['s0']
    # for b in ['SANTA_CRUZ','NEGRO','AMAZON','CLUTHA']: #NEGRO, AMAZON, SANTA CRUZ, CLUTHA
    if Basin_name in['SANTA_CRUZ','NEGRO','AMAZON','CLUTHA']:    
        ratio=ratio.shift(periods=182,freq='D')

    
    ratio_monthly = ratio.groupby([ratio.index.month]).mean()
    ratio_daily = ratio.groupby([ratio.index.month, ratio.index.day]).mean()
    # ratio_daily.plot(ax=ax1,legend=False)
    # ratio_monthly.plot(x=xticks,ax=ax2,legend=False)
    ratio_monthly = np.append(ratio_monthly.values,ratio_monthly.values[0])
    ax2.plot(xtickss,ratio_monthly,alpha=0)
    summer[Basin_name]=(np.mean(ratio_monthly[6:8]))
ax2.axhline(1,linestyle='--',color='black')
plt.xticks(xticks[:-1],monthcombis,rotation=30)
ax2.grid()
ax2.set_ylabel('Coupled model / Benchmark \n Ratio [-]')
ax2.text(xticks[-1]-10,3.2,'b)',size=15)
# ax2.set_xticks(ax2.get_xticks(),monthcombis)
# newticks=ax2.get_xticks()-10
# ax2.set_xticks(newticks)
# ax2.set_xticklabels(monthcombis,rotation=20)



ax3 = f1.add_axes([0,0,0.93,1])
ax3.set_visible(False)
im = ax3.imshow(np.array([[vmin,1]]),cmap=cmap)
bar =plt.colorbar(im,fraction=0.015,label=r'$Q_{99}$'+' glacier contribution [-]')



f1.savefig(join(figures,'ND+0Ratio_vectorplot.svg'),format='svg',bbox_inches = 'tight')


# xlabels=ax1.get_xticklabels()
# ax1.set_xticks(xticks-5,)

# ax1.set_title('Thjorsa')
# 
# bar.ax.set_yticklabels(['0','1','2','3','4','5'])
# plt.xticks(n)
# ax1.set_xlabel('Months')
# ax1.set_xticks()











# f1,ax1  = plt.subplots(figsize=(10,4))
# for Basin_name in Basin_names:
#     ratio = HG_dic[Basin_name]['s2']/HG_dic[Basin_name]['s0']
#     # for b in ['SANTA_CRUZ','NEGRO','AMAZON','CLUTHA']: #NEGRO, AMAZON, SANTA CRUZ, CLUTHA
#     if Basin_name in['SANTA_CRUZ','NEGRO','AMAZON','CLUTHA']:    
#         ratio=ratio.shift(periods=182,freq='D')

    
#     ratio_monthly = ratio.groupby([ratio.index.month]).mean()
#     ratio_daily = ratio.groupby([ratio.index.month, ratio.index.day]).mean()
#     # ratio_daily.plot(ax=ax1,legend=False)
#     ratio_monthly.plot(ax=ax1,legend=False)
# ax1.grid()
# ax1.set_ylabel('Ratio [-]')
# ax1.set_xtick


