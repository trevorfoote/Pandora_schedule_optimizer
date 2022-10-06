import os
from progressbar import ProgressBar
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

def target_vis(fdir, target_list, gmat_file, output_dir, \
    sun_block, moon_block, earth_block, obs_start, obs_stop):
    """ Determine visibility for target(s) with Pandora given avoidance angles
    for Sun, Moon, and Earth limb.
        
    Parameters
    ----------
    fdir:      string
                    directory containing input files 
    target_list:    string
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
    """

### Create an array of Pandora's science observation year with exactly 1 min spacings
    # datetime 
    Dt_iso_utc = pd.date_range(obs_start, obs_stop, freq='min')
    Dt_iso_utc = Dt_iso_utc.to_pydatetime()
    Dt_iso_utc = Time(Dt_iso_utc, format='datetime', scale='utc')

    # iso 
    t_iso_utc  = Dt_iso_utc.to_value('iso')
    T_iso_utc  = Time(t_iso_utc, format='iso', scale='utc')

    # mjd 
    t_mjd_utc = T_iso_utc.to_value('mjd')
    T_mjd_utc = Time(t_mjd_utc, format='mjd', scale='utc')


### Read in GMAT results
    print('Importing GMAT data')
    df = pd.read_csv(fdir+gmat_file, sep='\t')

    # Trim dataframe to slightly larger than date range of 
    # Pandora science lifetime defined as obs_start and obs_stop
    df = df[(df['Earth.UTCModJulian']>=(T_mjd_utc.jd[0]-2430000.0)-1) & 
            (df['Earth.UTCModJulian']<=(T_mjd_utc.jd[-1]-2430000.0)+1)]


### Convert GMAT times into standard MJD_UTC
    # Note: GMAT uses different offset for it's MJD (uses 2,430,000.0 rather than 2,400,000.5)
    # Time
    GMAT_jd_utc = np.array(df['Earth.UTCModJulian']) + 2430000.0
    GMAT_jd_utc = Time(GMAT_jd_utc, format='jd', scale='utc')
    GMAT_mjd_utc = GMAT_jd_utc.mjd
    GMAT_mjd_utc = Time(GMAT_mjd_utc, format='mjd', scale='utc')
    gmat_mjd_utc = GMAT_mjd_utc.to_value('mjd')


### Extract GMAT positions in MJ2000 Earth fixed cartesian coordinates
    # Earth
    GMAT_ex = np.array(df['Earth.EarthMJ2000Eq.X'])
    GMAT_ey = np.array(df['Earth.EarthMJ2000Eq.Y'])
    GMAT_ez = np.array(df['Earth.EarthMJ2000Eq.Z'])

    # Pandora
    GMAT_px = np.array(df['Pandora.EarthMJ2000Eq.X'])
    GMAT_py = np.array(df['Pandora.EarthMJ2000Eq.Y'])
    GMAT_pz = np.array(df['Pandora.EarthMJ2000Eq.Z'])

    # Sun
    GMAT_sx = np.array(df['Sun.EarthMJ2000Eq.X'])
    GMAT_sy = np.array(df['Sun.EarthMJ2000Eq.Y'])
    GMAT_sz = np.array(df['Sun.EarthMJ2000Eq.Z'])

    # Moon
    GMAT_mx = np.array(df['Luna.EarthMJ2000Eq.X'])
    GMAT_my = np.array(df['Luna.EarthMJ2000Eq.Y'])
    GMAT_mz = np.array(df['Luna.EarthMJ2000Eq.Z'])


### Interpolate all positions from GMAT to map to 1 minute intervals
    print('Interpolating locations at 1 minute intervals')
    # Earth
    ex = np.interp(t_mjd_utc, gmat_mjd_utc, GMAT_ex)
    ey = np.interp(t_mjd_utc, gmat_mjd_utc, GMAT_ey)
    ez = np.interp(t_mjd_utc, gmat_mjd_utc, GMAT_ez)

    # Pandora
    px = np.interp(t_mjd_utc, gmat_mjd_utc, GMAT_px)
    py = np.interp(t_mjd_utc, gmat_mjd_utc, GMAT_py)
    pz = np.interp(t_mjd_utc, gmat_mjd_utc, GMAT_pz)

    # Moon
    mx = np.interp(t_mjd_utc, gmat_mjd_utc, GMAT_mx)
    my = np.interp(t_mjd_utc, gmat_mjd_utc, GMAT_my)
    mz = np.interp(t_mjd_utc, gmat_mjd_utc, GMAT_mz)

    # Sun
    sx = np.interp(t_mjd_utc, gmat_mjd_utc, GMAT_sx)
    sy = np.interp(t_mjd_utc, gmat_mjd_utc, GMAT_sy)
    sz = np.interp(t_mjd_utc, gmat_mjd_utc, GMAT_sz)

### Coordinate shift to Pandora reference frame
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
    target_data = pd.read_csv(fdir + target_list, sep=',')
    
    #Cycle through targets
    for i in range(len(target_data['Simbad Name'])):
        target_name = target_data['Planet Name'][i]
        target_name_sc = target_data['Simbad Name'][i]
        target_sc = SkyCoord.from_name(target_name_sc)
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
        data = np.vstack((t_mjd_utc, Earth_req, Moon_req, Sun_req, all_req))
        data = data.T.reshape(-1,5)
        df1 = pd.DataFrame(data, columns = ['Time(MJD_UTC)','Earth_Clear','Moon_Clear','Sun_Clear','Visible'])

        output_file_name = 'Visibility for %s.csv' %target_name
        print('Saving results to ', save_dir)
        df1.to_csv((save_dir + output_file_name), sep=',', index=False)
