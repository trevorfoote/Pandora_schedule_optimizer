import os
from progressbar import ProgressBar
import numpy as np
import pandas as pd
import spectres
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.visualization import time_support
from datetime import datetime

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['font.size'] = 18
mpl.rcParams['axes.linewidth'] = 6
plt.rcParams['font.size'] = '18'


def target_vis(input_dir, target_file, gmat_file, output_dir, \
    sun_block, moon_block, earth_block, obs_start, obs_stop):
    """ Determine visibility for target(s) with Pandora given avoidance angles
    for Sun, Moon, and Earth limb.
        
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
    sun_block:      float
                    Avoidance angle for the Sun
    moon_block:     float
                    Avoidance angle for the Moon
    earth_block:    float
                    Avoidance angle for the Earth limb
    obs_start:      string
                    Date and time of start of Pandora science observing 
                    ex. '2025-04-25 00:00:00.000'
    obs_stop:      string
                    Date and time of end of Pandora science observing 
                    ex. '2026-04-25 00:00:00.000'
                                    
    Returns
    -------
    csv file
        file containing target's visibility by Pandora 
    plot
        plot of visibility schedule for year
    """

    
### Read in GMAT results
    print('Importing GMAT data')
    df = pd.read_csv(input_dir+gmat_file, sep='\t')


### Extract time and positions for objects of interest from GMAT file
    #Time
    GMAT_mjd_tdb = np.array(df['Earth.TDBModJulian'])

    #Earth
    GMAT_ex = np.array(df['Earth.EarthMJ2000Eq.X'])
    GMAT_ey = np.array(df['Earth.EarthMJ2000Eq.Y'])
    GMAT_ez = np.array(df['Earth.EarthMJ2000Eq.Z'])

    #Pandora
    GMAT_px = np.array(df['Pandora.EarthMJ2000Eq.X'])
    GMAT_py = np.array(df['Pandora.EarthMJ2000Eq.Y'])
    GMAT_pz = np.array(df['Pandora.EarthMJ2000Eq.Z'])

    #Sun
    GMAT_sx = np.array(df['Sun.EarthMJ2000Eq.X'])
    GMAT_sy = np.array(df['Sun.EarthMJ2000Eq.Y'])
    GMAT_sz = np.array(df['Sun.EarthMJ2000Eq.Z'])

    #Moon
    GMAT_mx = np.array(df['Luna.EarthMJ2000Eq.X'])
    GMAT_my = np.array(df['Luna.EarthMJ2000Eq.Y'])
    GMAT_mz = np.array(df['Luna.EarthMJ2000Eq.Z'])


### Convert GMAT times from MJD_TDB to MJD_UTC
 ## Note: GMAT uses different offset for it's MJD (GMAT uses 2,430,000.0 rather than 2,400,000.5)

    print('Converting GMAT times using astropy.Time')
    GMAT_jd_tdb = GMAT_mjd_tdb+2430000.0
    GMAT_jd_tdb = Time(GMAT_jd_tdb, format='jd', scale='tdb')
    GMAT_mjd_tdb = GMAT_jd_tdb.mjd
    GMAT_mjd_tdb = Time(GMAT_mjd_tdb, format='mjd', scale='tdb')
    GMAT_mjd_utc = GMAT_mjd_tdb.utc
    GMAT_mjd_utc = Time(GMAT_mjd_utc, format='mjd', scale='utc')
    GMAT_mjd_utc_np = GMAT_mjd_utc.to_value('mjd')

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

    # Fix end points
    if GMAT_mjd_utc_np[0] != t_mjd_utc[0]:
        GMAT_mjd_utc_np[0] = t_mjd_utc[0]
    if GMAT_mjd_utc_np[-1] != t_mjd_utc[-1]:
        GMAT_mjd_utc_np[-1] = t_mjd_utc[-1]


### Interpolate all positions from GMAT to map to 1 minute intervals
    print('Interpolating locations at 1 minute intervals')
    #Earth
    ex = spectres.spectres(t_mjd_utc, GMAT_mjd_utc_np, GMAT_ex, fill=GMAT_ex[-1], verbose=False)
    ey = spectres.spectres(t_mjd_utc, GMAT_mjd_utc_np, GMAT_ey, fill=GMAT_ey[-1], verbose=False)
    ez = spectres.spectres(t_mjd_utc, GMAT_mjd_utc_np, GMAT_ez, fill=GMAT_ez[-1], verbose=False)

    #Pandora
    px = spectres.spectres(t_mjd_utc, GMAT_mjd_utc_np, GMAT_px, fill=GMAT_px[-1], verbose=False)
    py = spectres.spectres(t_mjd_utc, GMAT_mjd_utc_np, GMAT_py, fill=GMAT_py[-1], verbose=False)
    pz = spectres.spectres(t_mjd_utc, GMAT_mjd_utc_np, GMAT_pz, fill=GMAT_pz[-1], verbose=False)

    #Moon
    mx = spectres.spectres(t_mjd_utc, GMAT_mjd_utc_np, GMAT_mx, fill=GMAT_mx[-1], verbose=False)
    my = spectres.spectres(t_mjd_utc, GMAT_mjd_utc_np, GMAT_my, fill=GMAT_my[-1], verbose=False)
    mz = spectres.spectres(t_mjd_utc, GMAT_mjd_utc_np, GMAT_mz, fill=GMAT_mz[-1], verbose=False)

    #Sun
    sx = spectres.spectres(t_mjd_utc, GMAT_mjd_utc_np, GMAT_sx, fill=GMAT_sx[-1], verbose=False)
    sy = spectres.spectres(t_mjd_utc, GMAT_mjd_utc_np, GMAT_sy, fill=GMAT_sy[-1], verbose=False)
    sz = spectres.spectres(t_mjd_utc, GMAT_mjd_utc_np, GMAT_sz, fill=GMAT_sz[-1], verbose=False)

### Coordinate shift to Pandora centric
    #Earth 
    exx = ex-px
    eyy = ey-py
    ezz = ez-pz

    #Pandora
    pxx = px-px
    pyy = py-py
    pzz = pz-pz

    #Sun
    sxx = sx-px
    syy = sy-py
    szz = sz-pz

    #Moon
    mxx = mx-px
    myy = my-py
    mzz = mz-pz


### Create SkyCoord for angular seperation calculations
    ss = SkyCoord(x=sxx, y=syy, z=szz, unit='km', representation_type='cartesian')
    ee = SkyCoord(x=exx, y=eyy, z=ezz, unit='km', representation_type='cartesian')
    mm = SkyCoord(x=mxx, y=myy, z=mzz, unit='km', representation_type='cartesian')
    pp = SkyCoord(x=pxx, y=pyy, z=pzz, unit='km', representation_type='cartesian')

### Define Constraints for each Solar system body
    Sun_constraint = sun_block * u.deg
    Moon_constraint = moon_block * u.deg
    Earth_constraint = earth_block *u.deg #based on worst case orbital alt of 450km should be 63 deg

    #Could be more robust to pull altitude from each time step but need to investigate
    # Pandora_alt = 450
    # Earth_constraint = np.arctan((1.*u.earthRad)/(1.*u.earthRad+Pandora_alt)).to(u.deg)

### Import Target list
    t_list = pd.read_csv(input_dir + target_file, sep=',')

    #Cycle through targets
    for i in range(len(t_list['Simbad Name'])):
        target_name = t_list['Simbad Name'][i]
        target_sc = SkyCoord.from_name(target_name)
        print('Analyzing constraints for:', target_name)

        #Evaluate at each time step whether target is blocked by each contraints
        Sun_req = np.zeros(len(pp))
        Moon_req = np.zeros(len(pp))
        Earth_req = np.zeros(len(pp))

        print('Calculating angular seperation requirements')
        pbar = ProgressBar()
        for i in pbar(range(len(pp))):
            Sun_req[i] = ss[i].separation(target_sc).deg * u.deg > Sun_constraint
            Moon_req[i] = mm[i].separation(target_sc).deg * u.deg > Moon_constraint
            Earth_req[i] = ee[i].separation(target_sc).deg * u.deg > Earth_constraint
        all_req = Sun_req * Moon_req * Earth_req

        #Check if folder exists for planet and if not create new folder for 
        #output products
        save_dir = output_dir + target_name + '/'
        if os.path.exists(save_dir) != True:
            os.makedirs(save_dir)
        
        #Save results for each planet to csv file
        data = np.vstack((np.array(list(t_mjd_utc[:])), Earth_req, Moon_req, Sun_req, all_req))
        data = data.T.reshape(-1,5)
        
        output_file_name = 'Visibility for %s.csv' %target_name
        print('Saving results to ', save_dir)
        np.savetxt((save_dir + output_file_name), data, delimiter=',',
                header='Time(MJD_UTC),Earth_Clear,Moon_Clear,Sun_Clear,Visible', comments='')
        

###Plot visibility calendar for each planet and save

        # Space out the individual requirements
        Sun_r = Sun_req * 0.25
        Moon_r = Moon_req * 0.50
        Earth_r = Earth_req * 0.75

        # Plot entire observing year
        fig, axs = plt.subplots(1,1, figsize=(20,10))
        axs.plot(T_iso_utc, Sun_r  , 'y|')
        axs.plot(T_iso_utc, Moon_r , 'C7|')
        axs.plot(T_iso_utc, Earth_r, 'g|')
        axs.plot(T_iso_utc, all_req, 'b|')

        axs.set_title('Observing Window for %s' %target_name)
        axs.xaxis.set_major_locator(mpl.dates.MonthLocator())
        axs.xaxis.set_minor_locator(mpl.dates.DayLocator())
        axs.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%b-%d'))
        axs.set_yticks([0,0.25, 0.5, 0.75, 1])
        axs.set_yticklabels(['Blocked', 'Sun Clear', 'Moon Clear', 'Earth Clear', 'Visible'])
        for label in axs.get_xticklabels(which='major'):
            label.set(rotation=30, horizontalalignment='right')
        axs.set_ylim(0.2, 1.05)
        axs.set_xlim(T_iso_utc[0], T_iso_utc[-1])

        # Save observing year plot
        plt.savefig(save_dir+'Observing Calendar for '
                    +target_name+'.png')
        plt.close()