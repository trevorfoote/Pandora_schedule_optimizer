import os
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from datetime import timedelta
import barycorr #.py file stored in same location as jupyter notebook

def transit_timing(fdir, target_list, planet_name, star_name):
    """ Determine primary transits for target(s) during Pandora's science 
    observation lifetime.
        
    Parameters
    ----------
    fdir:           string
                    high level directory 
    target_list:    string
                    name of csv file with list of targets
    planet_name:    string
                    name of target planet
    star_name:      string
                    name of target planet's host star
                                    
    Returns
    -------
    csv file
        file containing target's transits during Pandora's lifetime
    """

### Read in Visibility Data
    star_data = pd.read_csv(fdir + 'Targets/' + star_name + '/' + \
                    'Visibility for ' + star_name + '.csv', sep=',')
    t_mjd_utc = star_data['Time(MJD_UTC)']
    Visible = np.array(star_data['Visible'])

    # Convert time to datetime
    T_mjd_utc = Time(t_mjd_utc, format='mjd', scale='utc')
    T_iso_utc = Time(T_mjd_utc.iso, format='iso', scale='utc')
    dt_iso_utc = T_iso_utc.to_value('datetime')


### Extract planet specific parameters from target list
    target_data = pd.read_csv(fdir + target_list, sep=',')
    planet_name_sc = target_data.loc[target_data['Planet Name'] == planet_name,
                                'Planet Simbad Name'].iloc[0]
    planet_sc = SkyCoord.from_name(planet_name_sc)

    transit_dur = target_data.loc[target_data['Planet Name'] == planet_name, 
                            'Transit Duration (hrs)'].iloc[0] * u.hour
    period = target_data.loc[target_data['Planet Name'] == planet_name, 
                            'Period (day)'].iloc[0] *u.day
                            
    epoch_BJD_TDB = target_data.loc[target_data['Planet Name'] == planet_name, 
                'Transit Epoch (BJD_TDB-2400000.5)'].iloc[0]+2400000.5
    epoch_JD_UTC = barycorr.bjd2utc(epoch_BJD_TDB, planet_sc.ra.degree, planet_sc.dec.degree)
    epoch_JD_UTC = Time(epoch_JD_UTC, format='jd', scale= 'utc')
    epoch_MJD_UTC = epoch_JD_UTC.mjd
    epoch_MJD_UTC = Time(epoch_MJD_UTC, format='mjd', scale='utc')


### Calculate transit times
    # Calculate Pre mid-transit time on target
    half_obs_width = 0.75*u.hour + \
        np.maximum((1.*u.hour+transit_dur/2), transit_dur)

    # Determine minimum number of periods between Epoch and Pandora start of observing run
    min_start_epoch = epoch_MJD_UTC-half_obs_width
    min_pers_start = np.ceil((T_mjd_utc[0]-min_start_epoch)/period)

    # Calculate first transit within Pandora lifetime
    first_transit = epoch_MJD_UTC+(min_pers_start*period)

    # Calc transit times
    Mid_transits = np.arange(first_transit, T_mjd_utc[-1], period)
    for i in range(len(Mid_transits)):
        Mid_transits[i] = Mid_transits[i].mjd
    Mid_transits = (np.array(list(Mid_transits[:])))
    Mid_transits = Time(Mid_transits, format='mjd', scale= 'utc')

    Start_transits = Mid_transits-transit_dur/2
    End_transits   = Mid_transits+transit_dur/2

    start_transits = Start_transits.to_value('datetime')
    end_transits   = End_transits.to_value('datetime')

    # Truncate everything after the minutes place
    for i in range(len(start_transits)):
        start_transits[i] = start_transits[i] - timedelta(seconds=start_transits[i].second,
                                        microseconds=start_transits[i].microsecond)
        end_transits[i]   = end_transits[i] - timedelta(seconds=end_transits[i].second,
                                        microseconds=end_transits[i].microsecond)
    all_transits = np.arange(len(start_transits))


### Calculate which transits are visible to Pandora
    dt_vis_times = [] 
    for i in range(len(dt_iso_utc)):
        if Visible[i] == 1.0:
            dt_vis_times.append(dt_iso_utc[i])

    transit_coverage = np.zeros(len(start_transits))
    for i in range(len(start_transits)):
        tran_rng = pd.date_range(start_transits[i], end_transits[i], freq='min')
        tran_rng = tran_rng.to_pydatetime()
        
        tset = set(tran_rng)
        tran_test = tset.intersection(dt_vis_times)
        
        if len(tran_test)>0:
            transit_coverage[i] = len(tran_test)/len(tran_rng)
        else:
            continue
    
    #Check if folder exists for planet and if not create new folder for 
    #output products
    save_dir =  fdir + 'Targets/' + star_name + '/' + planet_name + '/'
    if os.path.exists(save_dir) != True:
        os.makedirs(save_dir)

### Save transit data to Visibility file
    transit_data = np.vstack((all_transits, Start_transits.value, Mid_transits.value, End_transits.value, transit_coverage))
    transit_data = transit_data.T.reshape(-1, 5)
    df = pd.DataFrame(transit_data, columns = ['All_Transits','Transit_Start','Transit_Mid','Transit_Stop','Transit_Coverage'])

    output_file_name = 'Visibility for ' + planet_name + '.csv'
    df.to_csv((save_dir + output_file_name), sep=',', index=False)