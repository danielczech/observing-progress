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
    parser.add_argument('p_dir', type = str, help = help_p_dir)
    parser.add_argument('t_obs', type = int, help = help_t_obs)
    parser.add_argument('n_beams', type = int, help = help_n_beams)
    parser.add_argument('d_min', type = int, help = help_d_min)
    if(len(sys.argv[1:]) == 0): 
        print('Missing args')
        parser.exit()
    args = parser.parse_args()
    main(p_dir = args.p_dir, t_obs = args.t_obs, 
            n_beams = args.n_beams)

def main(p_dir, t_obs, n_beams):
    # Set main index based on 32M database:
    observed = np.zeros(34000000)
    # Progress list:
    progress = []
    # Walk through files in directory
    p_files = os.listdir(p_dir)
    # For each file, step through pointings
    for p_file in p_files:
        p_path = os.path.join(p_dir, p_file)
        with open(p_path, 'rb') as f:
            p_list = pickle.load(f)
        for pointing in p_list:
            # Limit duration if desired (uncomment below):
            # if(pointing[1] < 300):
            #    continue
            # calculate number of beam slots available
            n_slots = n_beams*pointing[1]//t_obs # in seconds
            # check which stars still to be observed
            unobserved = np.where(observed[pointing[3]] == 0)[0]
            # mark as observed
            new_idxs = unobserved[0:int(n_slots)]
            observed[pointing[3][new_idxs]] = observed[pointing[3][new_idxs]] + 1
            # increment novel observed
            progress.append([pointing[2], len(new_idxs)])
    with open('progress.pkl', 'wb') as f:
        pickle.dump(progress, f)

if(__name__=="__main__"):
    cli()
