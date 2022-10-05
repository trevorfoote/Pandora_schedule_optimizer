import numpy as np
import pandas as pd
from astropy.time import Time

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['font.size'] = 18
mpl.rcParams['axes.linewidth'] = 6
plt.rcParams['font.size'] = '18'

### List of Graphs worth adding plotting function for
# Per Transit Coverage for GJ 3470b
# Visibility calendar



def Visibility_Calendar(fdir, target_name, output_dir=None, save=False):
    """ Plot Visibility Calendar for target
        
    Parameters
    ----------
    fdir:           string
                    directory containing input files 
    target_file:    string
                    name of csv file with list of targets
    output_dir:       string
                    directory where to save output plots if save_fig = True
    save_fig:       bool
                    Argument whether to save a copy of the figure or not
                                    
    Returns
    -------
    plot
        plot of transits schedule for year
    """
### Read in Visibility Data
    data = pd.read_csv(fdir + target_name + '/' + \
                       'Visibility for ' + target_name + '.csv', sep=',')
    t_mjd_utc = data['Time(MJD_UTC)']
    Earth_req = data['Earth_Clear']
    Moon_req = data['Moon_Clear']
    Sun_req = data['Sun_Clear']
    Visible = data['Visible']
    
    # Convert time to ISO_UTC
    T_mjd_utc = Time(t_mjd_utc, format='mjd', scale='utc')
    dt_iso_utc = T_mjd_utc.to_value('datetime')

    # Space out the individual requirements for plotting
    Sun_r = Sun_req * 0.25
    Moon_r = Moon_req * 0.50
    Earth_r = Earth_req * 0.75

    
### Plot entire observing year
    fig, axs = plt.subplots(1,1, figsize=(20,10))
    axs.plot(dt_iso_utc, Sun_r  , 'y|')
    axs.plot(dt_iso_utc, Moon_r , 'C7|')
    axs.plot(dt_iso_utc, Earth_r, 'g|')
    axs.plot(dt_iso_utc, Visible, 'b|')

    axs.set_title('Observing Window for %s' %target_name)
    axs.xaxis.set_major_locator(mpl.dates.MonthLocator())
    axs.xaxis.set_minor_locator(mpl.dates.DayLocator())
    axs.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%b-%d'))
    axs.set_yticks([0.25, 0.5, 0.75, 1])
    axs.set_yticklabels(['Sun Clear', 'Moon Clear', 'Earth Clear', 'Visible'])
    for label in axs.get_xticklabels(which='major'):
        label.set(rotation=30, horizontalalignment='right')
    axs.set_ylim(0.2, 1.05)
    axs.set_xlim(dt_iso_utc[0], dt_iso_utc[-1])

### Could add argument with specified start and stop times to plot
### But should major and minor axes change too?
### e.g.
### dt_iso_utc[np.where((datetime(2025, 5, 1, 0, 0)<=dt_iso_utc) & 
###                       (dt_iso_utc<=datetime(2025, 7, 1, 0, 0)))]
    
### Save observing year plot
    if save == True:
        plt.savefig(output_dir+'Observing Calendar for '
                +target_name+'.png')
    plt.show()
    plt.close()

def Transit_Plots(fdir, target_name, output_dir=None, save=False):
    data = pd.read_csv(fdir + target_name + '/' + \
                       'Visibility for ' + target_name + '.csv', sep=',')
    dt_iso_utc = Time(data['Time(MJD_UTC)'], format='mjd', scale='utc').to_value('datetime')
    Visible = data['Visible']

### Transit data
    tdata = data[~data['All_Transits'].isnull()]
    start_obs = Time(tdata['Observation_Start'], format='mjd', scale='utc').to_value('datetime')
    end_obs = Time(tdata['Observation_Stop'], format='mjd', scale='utc').to_value('datetime')
    start_transits = Time(tdata['Transit_Start'], format='mjd', scale='utc').to_value('datetime')
    end_transits = Time(tdata['Transit_Stop'], format='mjd', scale='utc').to_value('datetime')

    observation_coverage = tdata['Observation_Coverage']
    transit_coverage = tdata['Transit_Coverage']

    fig, axs = plt.subplots(1,1, figsize=(20,5))
    axs.plot(dt_iso_utc, Visible * 0.5, 'b|')
    axs.plot((start_obs, end_obs), (0.4, 0.4), '|', 
            color='orange', ms=20, mew=3)
    axs.plot((start_transits, end_transits), (0.4, 0.4), 'r-')

    axs.set_title('Observing Window for %s' %target_name)
    axs.xaxis.set_major_locator(mpl.dates.MonthLocator())
    axs.xaxis.set_minor_locator(mpl.dates.DayLocator())
    axs.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%b-%d'))
    axs.set_yticks([0.4, 0.5])
    axs.set_yticklabels(['Transit', 'Visible'])
    for label in axs.get_xticklabels(which='major'):
        label.set(rotation=30, horizontalalignment='right')
    axs.set_ylim(0.35, 0.55)
    axs.set_xlim(dt_iso_utc[0], dt_iso_utc[-1])

### Save calendar of transits
    if save == True:
        plt.savefig(output_dir+'Transit Calendar for '
                +target_name+'.png')
    plt.show()
    plt.close()


### Plot Transit Coverage
    plt.plot(transit_coverage*100, 'ro')
    # plt.plot(observation_coverage*100, 'bo') #option for plotting observation coverage too
    plt.ylim(0, 105)
    plt.xlabel('Transit Number')
    plt.ylabel('Coverage (%)')
    plt.title('Per Transit Coverage for %s' %target_name)

### Save calendar of transits
    if save == True:
        plt.savefig(output_dir+'Transit Coverage Plot for '
                +target_name+'.png')
    plt.show()
    plt.close()