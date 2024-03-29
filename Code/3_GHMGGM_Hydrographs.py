#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 10:00:50 2020

@author: Pau Wiersma

Scripts accompanying PCR-GLOBWB 2 and GloGEM coupling
(3/5): Hydrograph plotting and objective function (OF) calculation

This script:
    Loads the modelled hydrographs and plots them
        Single basin all years
        One year with a selection of basins
    Calculates the objective functions which can be plotted in the next script
        ND
        Ratio
        RRD
        Benchmark efficiency
        Flow duration curve efficiency
        Mass balance efficiency
        Standalone GHM evaluation
            NSE
            Calendar benchmark efficiency (Schaefli&Gupta2007)

Files needed:
    GHMGGM_basin_info.csv
    glaciers_nc files (including discharge observations)
    Model hydrographs 


Directories needed in run_dir:
    Files
        glacier_dailybasinsum
        glaciers_nc
    Figures
    Output (where hydrographs are )

Output: 
    OF_sorted.cv
    ratio.csv
    ND.csv
    RRD_means.csv


"""
import os
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import xarray as xr
import hydroeval as he
from matplotlib import cm




#%% Functions
def months_slice(year,hemisphere):
    """ Create a slice of all months in one hydrological year
    Southern hemispere is  months earlier
    """
    if hemisphere =='North':
        return slice(str(year-1)+'-10',str(year)+'-09')
    elif hemisphere =='South':
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

def FD_curve(data):
    """ Flow duration curve calculation"""
    y = data.sort_values(ascending=False).values
    x= np.cumsum(y)/np.nansum(y) *100
    return x,y

def load_hg2(GRDC_no, full_daterange):
    hg2_path = join(RUN_DIR,'Files','Q_obs',GRDC_no+'_Q_Day.Cmd.txt')
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

def RMSE(a,b):
    return np.sqrt(((a-b)**2).mean())

def BE(mod,obs,benchmark):
    """Calculation of the Benchmark efficiency (comparison
    between the coupled model and the benchmark against observations)"""
    return 1-(np.nansum((mod-obs)**2)/np.nansum((obs-benchmark)**2))
#%% INPUTS
RUN_DIR = r'd:\Documents\Master vakken\Thesis\Code'
os.chdir(RUN_DIR)
NC_DIR =r'd:\Documents\Master vakken\Thesis\Code\Files\glaciers_nc'
FIG_DIR             =join(RUN_DIR,'Figures')


basin_info = pd.read_csv(join(
    RUN_DIR,'Files','GHMGGM_basin_info.csv'),index_col = 0)

### 25 basins used in the paper
BASIN_NAMES = basin_info[basin_info['suitable']=='y'].index

MODEL_SETUPS        = ['2','0','1']
MODEL_NAMES         ={'0':'Modelled (Benchmark)',
                      '1':'Modelled (Bare)',
                      '2':'Modelled (Coupled)'}


MAIN_PLOT           =['s0','s1','s2'] #Which model settings to plot
ONLY_OBSYEARS       =True  #Consider only years for which GRDC is available
PLOT_BOOL           = False
SAVE_FIGS           = False
CALENDAR_DAY        =True #Calculate Caldendar day bencmark (Schaefli&Gupta2007)

GHM_NAME          ='PCRG'
GG_GCM              ='HadGEM2-ES'
GG_rcp              ='rcp26'
Fromyear            =2001
Untilyear           =2012

OF_list = []
normdiflist=[]
FD_list = []
RRD_list     = []
MBE_list = []
HG_list = []
glacier_sum_list =[]
Qobs_list = []
ratio_list =[]

#%%

for Basin_name in BASIN_NAMES:
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

    #Load GloGEM glacier runoff for analysis
    glacier_sum_path = join(RUN_DIR,'Files','glacier_dailybasinsum',Basin_name+'_glacsum.nc') #Write to this path if it doesn't exist
    nc_path = join(NC_DIR,'_'.join([Basin_name,
                                GG_GCM,
                                GG_rcp,
                                '2000',
                                '2016',
                                'R.nc']))
    nc_obs          =xr.open_dataset(nc_path).sel(time=daterange)

    if not os.path.isfile(glacier_sum_path):
        glacier_sum = nc_obs.R.sum(axis=(1,2))/(24*60*60)
        glacier_sum.to_netcdf(glacier_sum_path) #Write to file for later use
    else:
        glacier_sum = xr.open_dataarray(glacier_sum_path)
    glacier_sum_list.append(glacier_sum)

    #Load GRDC observations from .nc files
    if Basin_name in ['CLUTHA','COLUMBIA','RHINE','SUSITNA','DANUBE']:
        #NC files still has observations at basin mouth instead of more upstream
        full_date_range= pd.to_datetime(nc_obs.time.data)
        Qobs = load_hg2(str(int(basin_info.loc[Basin_name,'grdc_no'])),
                                full_date_range)
    else:
        Qobs       = nc_obs.hg



    Qobs = Qobs.where(Qobs!=-999,np.nan,drop=False)
    Qobs_list.append(Qobs)
    Q_time     = Qobs.where(~xr.ufuncs.isnan(Qobs),
                            drop=True).time
    Q_dates     = pd.to_datetime(Q_time.data)
    Q_years     = Q_dates.year.unique()

    #Load model runoff output and store in dictionary
    hg = {}
    for setup in MODEL_SETUPS:
        hgdic = {'Model_name'          :GHM_NAME,
           'Basin_name'          : Basin_name,
            'GG_GCM'              : GG_GCM ,
            'GG_rcp'              : GG_rcp,
           'Setting_no'          :setup,
           'Fromyear'           :'2001',
           'Untilyear'          :'2016',
           'station_no'         :'hg1'}
        hg_name = '_'.join(hgdic.values()) +'.txt'
        hg['s'+setup]=pd.read_csv(
                join(RUN_DIR,'Output','HG_2001_2016_',
                hg_name)
                ,index_col=0,
                parse_dates=True).loc[daterange]
    hg_years = daterange.year.unique()[1:]

    #-------------------------------------------------------------


    ###Calculate Objective functions and store everything in OF dictionary
    #make calendar day benchmark CDB
    if CALENDAR_DAY == True:
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

    OF = {}
    nanmask = ~xr.ufuncs.isnan(Qobs).data
    Obs = Qobs.data[nanmask]
    for setup in ['0','1','2']:
        Mod = hg['s'+setup].hg.values[nanmask]
        if CALENDAR_DAY ==True:
            OF['BE'+setup]  =BE(Mod,Obs,CDB.CDB.values[nanmask]) #Caldendar day bencmark efficiency
        OF['RMSE'+setup]=he.rmse(Mod,Obs)
        OF['NSE'+setup] =he.nse(Mod,Obs) #Nash sutcliffe
        # OF['KGE'+setup] =he.kge(Mod,Obs)[0][0]
        OF['MB'+setup]  =Mod.sum() - Obs.sum() #Mass balance error
    OF['BE']   = BE(hg['s2'].hg.values[nanmask],
                            Obs,hg['s0'].hg.values[nanmask]) #Overall benchmark efficiency

    if ('s0' in hg.keys())&('s2' in hg.keys()):
        Q02dif = hg['s2']-hg['s0']
        Qnormdif = Q02dif/Q02dif.quantile(0.99)
        # Qnormdif = Q02dif/Q02dif.max()
        normdiflist.append(Qnormdif.rename(columns={'hg':Basin_name}))

        Qratio = hg['s2']/hg['s0']
        ratio_list.append(Qratio.rename(columns={'hg':Basin_name}))



        #Calculate RRD
    if ('s0' in hg.keys())&('s1' in hg.keys())&('s2' in hg.keys()):
        #Make Dataframe with all hydrographs
        DF = pd.concat([hg['s0'],hg['s1'],hg['s2']],axis=1)
        DF.columns =  ['s0','s1','s2']
        DF['Obs'] = Qobs

        RRD = DF.groupby([DF.index.year,DF.index.month]).apply(
            lambda df,s0,s1,s2,obs:
                (RMSE(df[s0],df[obs])-RMSE(df[s2],df[obs]))/
                np.where(df[s2].sum()>df[s0].sum(),RMSE(df[s2],df[s1]),RMSE(df[s0],df[s1])),
                    's0','s1','s2','Obs')

        RRD.index = [pd.datetime.datetime(*RRD.index[i],1) for i in np.arange(len(RRD))]
        #Fix RRD-values that are just above 1 or below -1
        RRD[RRD>1] = 1
        RRD[RRD<-1]=-1

        #Mask all months that have on average less than 1 percent of annual glacier runoff
        glacsumDF = glacier_sum.to_dataframe()['R']
        pctmask = (glacsumDF.groupby([glacsumDF.index.month]).sum()/glacsumDF.sum())>0.01
        mask = [pctmask.loc[RRD.index[i].month] for i in range(len(RRD))]
        RRD=RRD.where(mask,np.nan)
        RRD_list.append(pd.DataFrame({Basin_name:RRD}))

    if ('s0' in hg.keys())&('s2' in hg.keys()):
        FD_bm = FD_curve(hg['s0'].hg)
        FD_cpl= FD_curve(hg['s2'].hg)
        FD_obs=FD_curve(Qobs.to_dataframe().hg)
        p = np.linspace(0,100,1001)
        FD_bm_int = np.interp(p,FD_bm[0],FD_bm[1])
        FD_cpl_int = np.interp(p,FD_cpl[0],FD_cpl[1])
        FD_obs_int = np.interp(p,FD_obs[0],FD_obs[1])
        FD_diff   =1- FD_bm_int/FD_cpl_int
        # FD_diff = (FD_bm_int-FD_cpl_int)/np.max(np.abs(FD_bm_int-FD_cpl_int))
        FD_list.append(pd.DataFrame({Basin_name:FD_diff},index=p))
        OF['FDBE'] = BE(FD_cpl_int,FD_obs_int,FD_bm_int)

    if 's1' in hg.keys():
        routed_R    = hg['s2']-hg['s1'] #routed glacier runoff
        glac_frac = routed_R.hg / hg['s2'].hg
        glac_frac_obs = routed_R.hg / Qobs.to_dataframe().hg
        max_frac  = glac_frac.max() #.max() or .quantile()
        OF['GFmax'] = max_frac
        OF['GC99']  = glac_frac.quantile(0.99)
        OF['GFmean'] = glac_frac.mean()
        basin_info.loc[Basin_name,'GC99'] = OF['GC99']

    OF_list.append(pd.DataFrame(OF,index=[Basin_name]))
    HG_list.append(hg)



    # Rsum = (nc_obs.R.sum(axis=(1,2))/(24*60*60))
    # frac = Rsum[0:-2]/Q2.hg[1:-1]

    # PCRG_glacR = Q0-Q1
    # routed_R   = Q2-Q1
    # Q02dif      =Q2-Q0

##%% Plotting
    if PLOT_BOOL:
            # Plotting parameters
        if ONLY_OBSYEARS==True:
            years = Q_years
        else:
            years=hg_years
        if years[0]==2000:
            years = years[1:]
        N = len(years)


        #% Regular plot
        multiplier = 1.7
        figheigth = multiplier*N #was 2.5
        topmargin = 1



        ### Main hydrograps plot
        f1,(axes) = plt.subplots(N,1,figsize=(10,figheigth),sharey=True)
        f1.subplots_adjust(top=0.8)
        plt.subplots_adjust(hspace=0.22) #was 0.21
        f1.suptitle(Basin_name.title(),
                    size='xx-large',y=(figheigth-0.15*topmargin)/figheigth)
        # colors = ['#d95f02','#7570b3','#1b9e77'] #colorblind safe
        #colors = ['tab:orange','tab:'red','tab:'green']

        turbo = cm.get_cmap('turbo',10)
        colors = turbo((9,1,5))


        for i in range(N): #One subplot per year
            ax = axes[i]
            months = months_slice(years[i],hemisphere)
            Qo = Qobs.sel(time=months)

            c=0
            zorders  =[2,1,0]
            for key,item in hg.items():
                if key in MAIN_PLOT:
                    Qm = item[months]
                    ax.plot(Qm.index,Qm,
                            linewidth=1.5,
                            label=MODEL_NAMES[key[1]],
                            color=colors[c],
                            zorder = zorders[c])
                    c+=1
            ax.plot(Qo.time,Qo,linewidth=1.5,alpha=0.8,label='Observed',
                    color='k',linestyle=(0,(5,1)))
            ax.set_ylim(0,ax.get_yticks()[-1])
            xticks = ax.get_xticks()
            ax.tick_params(axis='both', which='major', labelsize=8)
            ax.set_xlim(xticks[0]-30,xticks[-1]+29)
            ax.set_ylabel(r'$Q_{basin}\/[m^3/s]$')
            ax.grid()

            #Twin reverse axis for the glacier runoff barplot
            twin = ax.twinx()
            glacR = glacier_sum.sel(time=months)
            twin.bar(glacR.time,glacR,alpha=0.4,
                      width = 1.5, label='GloGEM',color='tab:blue')
            twin.set_ylabel(r'$Q_{glacier}\/ [m^3/s]$')
            twin.set_ylim(0,glacier_sum.data.max()*1.5)

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
        ##Make space above needed for title and legend
        make_space_above(axes,topmargin=topmargin)
        if SAVE_FIGS==True:
            f1.savefig(join(FIG_DIR,'Main_plots_colorblind',
                            '_'.join([Basin_name,'main'])+'.svg'),format='svg',
                        bbox_inches = 'tight')


        ###additional plot to visualize the metrics over the months
        # f1,(axes) = plt.subplots(N,1,figsize=(10,figheigth),sharey=True)
        # f1.subplots_adjust(top=0.8)
        # f1.suptitle(Basin_name.title()+' - Normalized difference by maximum',
        #             size='xx-large',y=(figheigth-0.15*topmargin)/figheigth)
        # for i in range(N):
        #     ax = axes[i]
        #     months = months_slice(years[i],hemisphere)
        #     DF[months].plot(ax=ax,
        #                     color=['tab:orange','tab:red','tab:green','tab:blue'],
        #                     legend=False)
        #     ax2=ax.twinx()
        #     RRD[months].plot(ax=ax2,color='black')
        #     Qnormdif[months].plot(ax=ax2,color='grey')
        #     ax.grid()
        #     ax2.axhline(0,color='black',linestyle='--')
        #     ax.set_ylim(0,ax.get_yticks()[-1])
        # make_space_above(axes,topmargin=0.4)

### Include OF results in figures :
    # OF_plot = ''
    # for i in MAIN_PLOT:
    #     F=''
    #     for key,item in OF.items():
    #         if i[1] in key:
    #             F+= key+' = '+'{:8.3f}'.format(item)+'  '
    #     F+= '\n'
    #     OF_plot +=F
    #     # OF_plot+= key+'{8.3f}'.format(item)
    # if 's1' in hg.keys():
    #     routed_R    = hg['s2']-hg['s1'] #routed glacier runoff
    #     glac_frac = routed_R.hg / hg['s2'].hg
    #     max_frac  = glac_frac.max() #.max() or .quantile()
    #     OF['GFmax'] = max_frac
    #     OF['GC99']  = glac_frac.quantile(0.99)
    #     OF['GFmean'] = glac_frac.mean()
    #+ this part of code for the plotting
        # #     if i==0:
    # #         lines, labels = ax.get_legend_handles_labels()
    # #         ax.legend(
    # #             lines , labels, loc=(0.8, 1.1), prop={"size": 8},
    # #         )
    # #         ax.text(0,1.1,OF_plot,transform = ax.transAxes)


#%% Organize OFs and save to files for next script
# Concat and sort OF's

basin_info.to_csv(join(RUN_DIR,'Files','GHMGGM_basin_info.csv'))

###OFs
OF_df = pd.concat(OF_list)
OF_df = OF_df.reindex(sorted(OF_df.columns,reverse=False),axis=1)

OF_sorted= OF_df.sort_values(by=['GC99'],axis=0,ascending=False)
OF_sorted.to_csv(join(RUN_DIR,'Output','OF_sorted'+str(len(OF_list))+'.csv'))



###ND
normdif_stack = pd.concat(normdiflist,axis=1)
#Shift southern hemisphere by6 months
for b in ['SANTA_CRUZ','NEGRO','AMAZON','CLUTHA']: #NEGRO, AMAZON, SANTA CRUZ, CLUTHA
    if b in normdif_stack.keys():
        normdif_stack[b]=normdif_stack[b].shift(periods=182,freq='D')

normdif_means = normdif_stack.groupby([normdif_stack.index.month, normdif_stack.index.day]).mean()
normdif_means.index=normdif_stack['2004'].index #Random year, it's only about the calendar days


normdif_means.columns = OF_df.GC99.values
normdif_means = normdif_means.reindex(sorted(normdif_means.columns),axis=1)
normdif_means.to_csv(join(RUN_DIR,'Output','ND_'+str(len(OF_list))+'.csv'))

### Ratio
ratio_stack = pd.concat(ratio_list,axis=1)
if Basin_name in['SANTA_CRUZ','NEGRO','AMAZON','CLUTHA']:
    if b in ratio_stack.keys():
        ratio=ratio_stack.shift(periods=182,freq='D')

ratio_monthly = ratio_stack.groupby([ratio_stack.index.month]).mean()
ratio_monthly = (ratio_monthly -1)*100
ratio_monthly.to_csv(join(RUN_DIR,'Output','ratio_'+str(len(OF_list))+'.csv'))

###RRD
RRD_concat = pd.concat(RRD_list,axis=1)
#Shift Southern hemisphere  months
for b in ['NEGRO','SANTA_CRUZ','AMAZON','CLUTHA']:
    if b in RRD_concat.keys():
        RRD_concat[b] = RRD_concat[b].shift(periods=6,freq='MS')

RRD_means = RRD_concat.groupby([RRD_concat.index.month]).mean()
RRD_means.to_csv(join(RUN_DIR,'Output','RRD_means.csv'))



#%%# Plot one year of multiple basin in one figure

HG_dic = {BASIN_NAMES[i]:HG_list[i] for i in range(len(BASIN_NAMES))}
Qobs_dic = {BASIN_NAMES[i]:Qobs_list[i] for i in range(len(BASIN_NAMES))}
glacier_sum_dic = {BASIN_NAMES[i]:glacier_sum_list[i] for i in range(len(BASIN_NAMES))}

# Basins = ['OELFUSA','ALSEK','RHONE','COLUMBIA','RHINE','MACKENZIE']
Basins = ['JOEKULSA','SANTA_CRUZ','NELSON','LULE']
NB = len(Basins)

# year = 2005
for year in range(2010,2011):

    f1,axes = plt.subplots(NB,1,figsize=(10,1.9*NB),sharex=True)
    # f1.subplots_adjust(top=0.8)
    plt.subplots_adjust(hspace=0.31)
    # f1.suptitle(Basin_name[0]+Basin_name[1:].lower(),
    #             size='xx-large',y=(figheigth-0.15*topmargin)/figheigth)
    # colors = ['tab:orange','tab:red','tab:green']
    colors = ['#d95f02','#7570b3','#1b9e77'] #colorblind safe
    colors = cm.get_cmap('turbo',10)((9,1,5))
    for i in range(NB):
        hg = HG_dic[Basins[i]]
        Qobs = Qobs_dic[Basins[i]]
        glacier_sum = glacier_sum_dic[Basins[i]]

        ax = axes[i]
        months = months_slice(year,hemisphere)
        Qo = Qobs.sel(time=months)


        c=0
        zorders = [2,1,0]
        for key,item in hg.items():
            if key in MAIN_PLOT:
                Qm = item[months]
                ax.plot(Qm.index,Qm,
                        linewidth=1.5,
                        label=MODEL_NAMES[key[1]],
                        color=colors[c],
                        zorder = zorders[c])
                c+=1
        ax.plot(Qo.time,Qo,linewidth=1.5,alpha=0.8,label='Observed',
                color='k',linestyle=(0,(5,1)))
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
        ax.grid(alpha=0.6)
        gd = OF_sorted.loc[Basins[i],'GC99']
        # gd = basin_info.loc[Basins[i],'glac_degree']
        if i==0:
            # ax.set_title(Basins[i].title()+' (Glaciation degree ='+str(round(gd,2))+'%)')
            # ax.set_title(Basins[i].title()+r' ($fQ_{99}$'+' = '+str(round(gd,2))+')')
            ax.set_title(Basins[i].title()+r' (GC99 '+' = '+str(round(gd,2))+')')
        else:
            ax.set_title(Basins[i].title()+' ('+str(round(gd,2))+')')
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

        if Basins[i]=='ALSEK':
            lines, labels = ax.get_legend_handles_labels()
            lines2, labels2 = twin.get_legend_handles_labels()
            twin.legend(
                lines + lines2, labels + labels2,
                loc=(0.19, 0.102), prop={"size": 8}) #loc=(0.8,1.1)

        if Basins[i]=='LULE':
            lines, labels = ax.get_legend_handles_labels()
            lines2, labels2 = twin.get_legend_handles_labels()
            twin.legend(
                lines + lines2, labels + labels2,
                loc=(0.35, 4.13), prop={"size": 7.5}) #loc=(0.8,1.1)

        if i==NB-1:
            binitials = ''.join([B[0] for B in Basins ])
            if SAVE_FIGS==True:
                f1.savefig(join(FIG_DIR,'Ensembleplots','Ensembleplot'+binitials+str(year)+'_nowhitespace_GC99.svg'),format='svg',
                bbox_inches = 'tight')







