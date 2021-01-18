import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
import time
import pickle
import argparse
import sys
import pandas as pd
import os
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

def cli(args = sys.argv[0]):
    usage = "{} [options]".format(args)
    description = """Get observing progress (cumulative count)."""
    parser = argparse.ArgumentParser(prog = 'star-lists', 
        usage = usage, description = description) 
    help_p_dir = """Specify directory of pre-processed pointing files (each 
    containing duration, date and list of stars)"""
    help_t_obs = """Desired observation duration for each star"""
    help_n_beams = """Specified number of beams"""
    help_d_min = """Minimum acceptable primary observation duration"""
    help_cat = """Main catalogue of stars which was drawn from"""
    parser.add_argument('p_dir', type = str, help = help_p_dir)
    parser.add_argument('t_obs', type = int, help = help_t_obs)
    parser.add_argument('n_beams', type = int, help = help_n_beams)
    parser.add_argument('d_min', type = int, help = help_d_min)
    parser.add_argument('catalogue', type = str, help = help_cat)
    if(len(sys.argv[1:]) == 0): 
        print('Missing args')
        parser.exit()
    args = parser.parse_args()
    main(p_dir = args.p_dir, t_obs = args.t_obs, 
            n_beams = args.n_beams, d_min = args.d_min, catalogue = args.catalogue)

def accumulate_dates(dates, stars):
    """Accumulate observed stars on the same dates.
    
    Args:
        dates: list of datetime objects
        stars: list of associated numbers of observed stars

    Returns:
        a_dates: list of unique datetime objects
        a_stars: list of associated accumulated numbers of stars
    """
    start = min(dates)
    stop = max(dates)
    t_range = (stop - start).days
    a_dates = [start + timedelta(days = n) for n in range(t_range + 1)]
    a_stars = [0 for n in range(t_range + 1)]
    for i in range(len(dates)):
        idx = (dates[i] - start).days
        a_stars[idx] = a_stars[idx] + stars[i]
    return a_dates, a_stars

def main(p_dir, t_obs, n_beams, d_min, catalogue):
    VERBOSE = False
    # Set main index based on 32M database:
    observed = np.zeros(34000000)
    # Progress list:
    dates = []
    stars = []
    # Walk through files in directory
    p_files = os.listdir(p_dir)
    # Sort earliest first:
    p_files = sorted(p_files)
    # For each file, step through pointings
    for p_file in p_files:
        p_path = os.path.join(p_dir, p_file)
        with open(p_path, 'rb') as f:
            p_list = pickle.load(f)
        for pointing in p_list:
            # Limit duration if desired:
            if(pointing[1] < d_min):
                continue
            # Calculate number of beam slots available
            # Note: Partial slots are discarded here (partial being < d_min)
            p_duration = (pointing[1]//d_min)*d_min
            n_slots = (n_beams*p_duration)//t_obs # t_obs in seconds
            # check which stars still to be observed
            unobserved = np.where(observed[pointing[3]] == 0)[0]
            # mark as observed
            new_idxs = unobserved[0:int(n_slots)]
            observed[pointing[3][new_idxs]] = observed[pointing[3][new_idxs]] + 1
            if(VERBOSE):
                print("    Duration: {} Slots: {}".format(pointing[1], n_slots))
                print("    New stars: {}".format(len(new_idxs)))
            # format date/time for accumulation later
            ep_time = datetime.strptime(pointing[2], '%Y-%m-%d')
            # increment novel observed
            dates.append(ep_time) 
            stars.append(len(new_idxs))
    print("Saving...")
    a_dates, a_stars = accumulate_dates(dates, stars)
    progress = [a_dates, a_stars]
    with open('progress_{}_{}_a.pkl'.format(n_beams, d_min), 'wb') as f:
        pickle.dump(progress, f)
    # Access statistics on observed stars
    # Load catalogue
    catalogue = pd.read_csv(catalogue, delimiter = ',',
            dtype={'source_id':str,
                'ra':float,
                'dec':float,
                'dist_c':float})
    maindb = catalogue.to_numpy()
    maindb = maindb[1:, :] # remove column headers
    # Use lines below if using full unrestricted DB:
    # maindb_ = np.zeros((maindb.shape[0], 3))
    # maindb_[:, 1:2] = maindb[:, 3:4]
    # maindb_[:, 2:3] = maindb[:, 5:6]
    # Sort by distance:
    dist_idx = np.argsort(maindb[:, -1])
    db_full = maindb[dist_idx, :]
    print("Catalogue loaded")
    # Stars to consider:
    star_idx = np.where(observed > 0)[0]
    distances = maindb[:, 3][star_idx]
    print(np.min(distances))
    print(np.max(distances))

if(__name__=="__main__"):
    cli()
