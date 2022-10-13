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
    t_mjd_utc  = T_iso_utc.to_value('mjd')
    T_mjd_utc  = Time(t_mjd_utc, format='mjd', scale='utc')


### Read in GMAT results
    print('Importing GMAT data')
    gmat_data = pd.read_csv(fdir+gmat_file, sep='\t')

    # Trim dataframe to slightly larger than date range of 
    # Pandora science lifetime defined as obs_start and obs_stop
    gmat_data = gmat_data[(gmat_data['Earth.UTCModJulian']>=(T_mjd_utc.jd[0]-2430000.0)-0.0007) & 
            (gmat_data['Earth.UTCModJulian']<=(T_mjd_utc.jd[-1]-2430000.0)+0.0007)]


### Convert GMAT times into standard MJD_UTC
    # Note: GMAT uses different offset for it's MJD (uses 2,430,000.0 rather than 2,400,000.5)
    # Time
    GMAT_jd_utc  = np.array(gmat_data['Earth.UTCModJulian']) + 2430000.0
    GMAT_jd_utc  = Time(GMAT_jd_utc, format='jd', scale='utc')
    GMAT_mjd_utc = Time(GMAT_jd_utc.mjd, format='mjd', scale='utc')
    gmat_mjd_utc = GMAT_mjd_utc.to_value('mjd')


### Extract GMAT positions in MJ2000 Earth fixed cartesian coordinates
    # Earth
    gmat_ex = np.array(gmat_data['Earth.EarthMJ2000Eq.X'])
    gmat_ey = np.array(gmat_data['Earth.EarthMJ2000Eq.Y'])
    gmat_ez = np.array(gmat_data['Earth.EarthMJ2000Eq.Z'])

    # Pandora
    gmat_px = np.array(gmat_data['Pandora.EarthMJ2000Eq.X'])
    gmat_py = np.array(gmat_data['Pandora.EarthMJ2000Eq.Y'])
    gmat_pz = np.array(gmat_data['Pandora.EarthMJ2000Eq.Z'])

    # Sun
    gmat_sx = np.array(gmat_data['Sun.EarthMJ2000Eq.X'])
    gmat_sy = np.array(gmat_data['Sun.EarthMJ2000Eq.Y'])
    gmat_sz = np.array(gmat_data['Sun.EarthMJ2000Eq.Z'])

    # Moon
    gmat_mx = np.array(gmat_data['Luna.EarthMJ2000Eq.X'])
    gmat_my = np.array(gmat_data['Luna.EarthMJ2000Eq.Y'])
    gmat_mz = np.array(gmat_data['Luna.EarthMJ2000Eq.Z'])

    # Pandora's Lat and Lon
    gmat_lat = np.array(gmat_data['Pandora.Earth.Latitude'])
    gmat_lon = np.array(gmat_data['Pandora.Earth.Longitude'])


### Interpolate all positions from GMAT to map to 1 minute intervals
    print('Interpolating locations at 1 minute intervals')
    # Earth
    ex = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_ex)
    ey = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_ey)
    ez = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_ez)

    # Pandora
    px = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_px)
    py = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_py)
    pz = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_pz)

    # Moon
    mx = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_mx)
    my = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_my)
    mz = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_mz)

    # Sun
    sx = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_sx)
    sy = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_sy)
    sz = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_sz)

    # Pandora Latitude & Longitude
    gmat_lon_cont = np.copy(gmat_lon)
    for i in range(len(gmat_lon_cont)-1):
        if np.abs(gmat_lon_cont[i+1]-gmat_lon_cont[i]) > 100:
            gmat_lon_cont[i+1:] = gmat_lon_cont[i+1:] - 360
    
    p_lon_cont = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_lon_cont)
    p_lon = np.copy(p_lon_cont)
    for i in range(len(p_lon)):
        while p_lon[0] < -180:
            p_lon = p_lon + 360        
        if p_lon[i] < -180:
            p_lon[i:] = p_lon[i:] + 360
    p_lon = p_lon * u.deg
    p_lat = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_lat) * u.deg

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

### Evaluate at each time step whether Pandora is crossing SAA
    # SAA coordinates at altitude of ~500km
    saa_lat_max = 0. * u.deg
    saa_lat_min = -50. * u.deg
    saa_lon_max = 40. * u.deg
    saa_lon_min = -90. * u.deg
    saa_cross   = np.zeros(len(p_lat))
    for i in range(len(p_lat)):
        if (saa_lat_min <= p_lat[i]) and (p_lat[i] <= saa_lat_max) and \
        (saa_lon_min <= p_lon[i]) and (p_lon[i] <= saa_lon_max):
            saa_cross[i] = 1.


### Import Target list
    target_data = pd.read_csv(fdir + target_list, sep=',')
    
    #Cycle through targets
    for i in range(len(target_data['Simbad Name'])):
        target_name    = target_data['Planet Name'][i]
        target_name_sc = target_data['Simbad Name'][i]
        target_sc      = SkyCoord.from_name(target_name_sc)
        print('Analyzing constraints for:', target_name)

        #Evaluate at each time step whether target is blocked by each contraints
        Sun_sep   = np.zeros(len(pp))
        Moon_sep  = np.zeros(len(pp))
        Earth_sep = np.zeros(len(pp))
        Sun_req   = np.zeros(len(pp))
        Moon_req  = np.zeros(len(pp))
        Earth_req = np.zeros(len(pp))

        print('Calculating angular seperation requirements')
        pbar = ProgressBar()
        for i in pbar(range(len(pp))):
            Sun_sep[i]   = ss[i].separation(target_sc).deg
            Sun_req[i]   = Sun_sep[i]* u.deg > Sun_constraint
            Moon_sep[i]  = mm[i].separation(target_sc).deg
            Moon_req[i]  = Moon_sep[i]* u.deg > Moon_constraint
            Earth_sep[i] = ee[i].separation(target_sc).deg
            Earth_req[i] = Earth_sep[i]* u.deg > Earth_constraint
        all_req = Sun_req * Moon_req * Earth_req

        #Check if folder exists for planet and if not create new folder for 
        #output products
        save_dir = output_dir + target_name + '/'
        if os.path.exists(save_dir) != True:
            os.makedirs(save_dir)
        
        #Save results for each planet to csv file
        data = np.vstack((t_mjd_utc, saa_cross, Earth_req, Moon_req, Sun_req, \
            all_req, Earth_sep, Moon_sep, Sun_sep))
        data = data.T.reshape(-1,9)
        df1  = pd.DataFrame(data, columns = ['Time(MJD_UTC)','SAA_Crossing','Earth_Clear','Moon_Clear','Sun_Clear',\
            'Visible','Earth_Sep','Moon_Sep','Sun_Sep'])

        output_file_name = 'Visibility for %s.csv' %target_name
        print('Saving results to ', save_dir)
        df1.to_csv((save_dir + output_file_name), sep=',', index=False)
