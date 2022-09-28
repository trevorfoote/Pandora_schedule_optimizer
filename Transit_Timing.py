from mpl_toolkits import mplot3d
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from datetime import datetime, timedelta
import barycorr #.py file stored in same location as jupyter notebook

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

def transit_timing(input_dir, target_file, output_dir, obs_start, obs_stop):
    """ Determine primary transits for target(s) during Pandora's science 
    observation lifetime.
        
    Parameters
    ----------
    input_dir:      string
                    directory containing input files 
    target_file:    string
                    name of csv file with list of targets
    gmat_file:      string
                    name of txt file from GMAT containing time and positions
                    of Pandora, Sun, Moon, and Earth
    output_dir:     string
                    directory where to save output csv and plots
    obs_start:      string
                    Date and time of start of Pandora science observing 
                    ex. '2025-04-25 00:00:00.000'
    obs_stop:      string
                    Date and time of end of Pandora science observing 
                    ex. '2026-04-25 00:00:00.000'
                                    
    Returns
    -------
    csv file
        file containing target's transits during Pandora's lifetime
    plot
        plot of transits schedule for year
    """

### Create an array of Pandora's science observation year with exactly 1 min spacings
    sci_start = Time(obs_start, format='iso', scale='utc') #Pandora science start in Greg_UTC
    sci_stop  = Time(obs_stop, format='iso', scale='utc') #Pandora science stop in Greg_UTC
    t_iso_utc = np.arange(sci_start, sci_stop+(1*u.minute), 1*u.minute)
    t_mjd_utc = t_iso_utc
    for i in range(len(t_iso_utc)):
        t_mjd_utc[i] = t_iso_utc[i].mjd

    # Create a time array in Gregarian that can be used in matplotlib
    T_mjd_utc = np.array(list(t_mjd_utc[:]))
    T_mjd_utc = Time(T_mjd_utc, format='mjd', scale= 'utc')

    T_iso_utc =[]
    for i in range(len(T_mjd_utc)):
        T_iso_utc.append(datetime.strptime(T_mjd_utc.iso[i], '%Y-%m-%d %H:%M:%S.%f'))
