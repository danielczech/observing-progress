import numpy as np
import time
import pickle
import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt

"""
Script to calculate some basic statistics from a list of
completed pointings.
"""

def cli(args = sys.argv[0]):
    usage = "{} [options]".format(args)
    description = """Calculate some basic stats from pointing list"""
    parser = argparse.ArgumentParser(prog = 'star-lists',
        usage = usage, description = description)
    help_info = """Specify catalogue, starting position in file, 
                number of lines, and minimum obs duration."""
    parser.add_argument('pointings', type = str, help = help_info)
    if(len(sys.argv[1:]) == 0):
        print('Missing args')
        parser.exit()
    args = parser.parse_args()
    main(pointings = args.pointings)

def main(pointings):
    # Load file of VLITE pointings
    # For now, start with RA, Dec, duration and date.
    pd.options.display.max_colwidth = 2000 # Surprising this is necessary...
    pointings = pd.read_csv(pointings,
                            usecols = [0,1,2,3],
                            delimiter = '  ',
                            dtype={'RA_[deg]':float,
                                'Dec_[deg]':float,
                                'Duration_[s]':float,
                                'YYYY-MM-DD':str})
    pointings = pointings.to_numpy()
    with open('durations.pkl', 'wb') as f:
        pickle.dump(pointings[:, 2], f)
    plt.hist(pointings[:, 2])
    plt.show()

if(__name__ == '__main__'):
    cli()
