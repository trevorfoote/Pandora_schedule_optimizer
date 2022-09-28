from astropy import units as u

import Visibility_Calc as VC
import Transit_Timing as TT

### Establish Directories and files
input_dir   = '/home/trevor/OneDrive/Coding/Pandora/Schedule_Optimizer/'
target_file = 'Pandora_target_list.csv'
gmat_file   = 'GMAT_Pandora.txt'
output_dir  = '/home/trevor/OneDrive/Coding/Pandora/Schedule_Optimizer/Results/'

### Establish observation parameters
sun_block   = 90. *u.deg
moon_block  = 40. *u.deg
earth_block = 63. *u.deg
obs_start = '2025-04-25 00:00:00.000'
obs_stop  = '2026-04-25 00:00:00.000'

### Calculate blocking regions from Earth-limb, Moon, & Sun
VC.target_vis(input_dir, target_file, gmat_file, output_dir, \
    sun_block, moon_block, earth_block, obs_start, obs_stop)

### Calculate all the transits that will occur for each Target Planet during Pandora's lifetime
TT.transit_timing(input_dir, target_file, output_dir, obs_start, obs_stop)