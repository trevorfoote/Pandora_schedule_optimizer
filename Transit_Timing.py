import os
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from datetime import timedelta
import barycorr #.py file stored in same location as jupyter notebook

def transit_timing(fdir, target_list, planet_name, star_name, output_dir):
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
    output_dir:     string
                    directory where to save output csv and plots
                                    
    Returns
    -------
    csv file
        file containing target's transits during Pandora's lifetime
    """

### Read in Visibility Data
    data = pd.read_csv(fdir + 'Results/' + star_name + '/' + \
                    'Visibility for ' + star_name + '.csv', sep=',')
    t_mjd_utc = data['Time(MJD_UTC)']
    Visible = np.array(data['Visible'])

    # Convert time to datetime
    T_mjd_utc = Time(t_mjd_utc, format='mjd', scale='utc')
    T_iso_utc = T_mjd_utc.iso
    T_iso_utc = Time(T_iso_utc, format='iso', scale='utc')
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
    End_transits = Mid_transits+transit_dur/2
    Start_obs = Mid_transits-half_obs_width
    End_obs = Mid_transits+half_obs_width

    start_transits = Start_transits.to_value('datetime')
    end_transits = End_transits.to_value('datetime')
    start_obs = Start_obs.to_value('datetime')
    end_obs = End_obs.to_value('datetime')

    # Truncate everything after the minutes place
    for i in range(len(Start_transits)):
        start_transits[i] = start_transits[i] - timedelta(seconds=start_transits[i].second,
                                        microseconds=start_transits[i].microsecond)
        end_transits[i] = end_transits[i] - timedelta(seconds=end_transits[i].second,
                                        microseconds=end_transits[i].microsecond)
        start_obs[i] = start_obs[i] - timedelta(seconds=start_obs[i].second,
                                        microseconds=start_obs[i].microsecond)
        end_obs[i] = end_obs[i] - timedelta(seconds=end_obs[i].second,
                                        microseconds=end_obs[i].microsecond)
    all_transits = np.arange(len(start_transits))


### Calculate which transits are visible to Pandora
    dt_vis_times = [] 
    for i in range(len(dt_iso_utc)):
        if Visible[i] == 1.0:
            dt_vis_times.append(dt_iso_utc[i])

    obs_coverage = np.zeros(len(start_obs))
    transit_coverage = np.zeros(len(start_transits))

    for i in range(len(start_obs)):
        temp = pd.date_range(start_obs[i], end_obs[i], freq='min')
        obs_rng = temp.to_pydatetime()

        temp = pd.date_range(start_transits[i], end_transits[i], freq='min')
        tran_rng = temp.to_pydatetime()
        
        oset = set(obs_rng)
        tset = set(tran_rng)
        obs_test = oset.intersection(dt_vis_times)
        tran_test = tset.intersection(dt_vis_times)
        
        if len(obs_test)>0 and len(tran_test)>0:
            obs_coverage[i] = len(obs_test)/len(obs_rng)
            transit_coverage[i] = len(tran_test)/len(tran_rng)
        else:
            continue
    
    #Check if folder exists for planet and if not create new folder for 
    #output products
    save_dir = output_dir + star_name + '/' + planet_name + '/'
    if os.path.exists(save_dir) != True:
        os.makedirs(save_dir)

### Save transit data to Visibility file
    transit_data = np.vstack((all_transits, Start_transits.value, Mid_transits.value, End_transits.value, Start_obs.value, End_obs.value, obs_coverage, transit_coverage))
    transit_data = transit_data.T.reshape(-1, 8)
    # df1 = data[['Time(MJD_UTC)','Earth_Clear','Moon_Clear','Sun_Clear','Visible']].copy()
    # df2 = pd.DataFrame(transit_data, columns = ['All_Transits','Transit_Start','Transit_Mid','Transit_Stop','Observation_Start','Observation_Stop','Observation_Coverage','Transit_Coverage'])
    # result = pd.concat([df1,df2], axis=1)

    output_file_name = 'Visibility for ' + planet_name + '.csv'
    # result.to_csv((save_dir + output_file_name), sep=',', index=False)
    df2.to_csv((save_dir + output_file_name), sep=',', index=False)