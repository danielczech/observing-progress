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

def convert_time(MJD_array):
    """
    Convert an array of MJDs into UNIX timestamps.
    """
    ts_array = (MJD_array + 2400000.5 - 2440587.5)*86400.0
    return ts_array

def gap_durations(given, MJDstart, MJDend):
    """
    Calculate pointing gap durations.
    Convert to unix timestamp in seconds.
    """
    MJDstart_i = convert_time(MJDstart[0:10000])
    MJDend_i = convert_time(MJDend[0:10000])
    gaps = np.subtract(MJDstart_i[1:], MJDend_i[0:len(MJDstart_i) - 1])
    #for i in range(2000):
    #    print('{}    {}    {}'.format(MJDstart[i], MJDend[i], gaps[i])) 
    durations = np.subtract(MJDend_i, MJDstart_i)
    plt.plot(durations)
    plt.show()
    plt.plot(gaps)
    plt.ylabel('Gap duration (s)')
    plt.xlabel('Pointing number')
    plt.show()
    plt.hist(gaps, bins = 500, facecolor = 'w', edgecolor = 'k', hatch = '/')
    plt.ylabel('Number of pointings')
    plt.xlabel('Gap Duration (s)')
    plt.yscale('log')
    plt.show()

def fov(band, D):
    """
    Estimate field of view from frequency band (str).
    Expects naming format: 1.5GHz
    Returns radius of FoV in degrees.
    """
    if('M' in band):
        band = band.split('M')[0]
        band = float(band)*10**6
    elif('G' in band):
        band = band.split('G')[0]
        band = float(band)*10**9
    else:
        print('Unable to read band: {}'.format(band))
    wavelength = 299792458.0/band
    field_of_view = 1.02*wavelength/D*180.0/np.pi/2.0
    return field_of_view

def main(pointings):
    # Load file of VLITE pointings
    # For now, start with RA, Dec, duration and date.
    pd.options.display.max_colwidth = 2000 # Surprising this is necessary...
    pointings = pd.read_csv(pointings,
                            usecols = [0,1,2,3,4,5,10],
                            delimiter = '  ',
                            dtype={'RA_[deg]':float,
                                'Dec_[deg]':float,
                                'Duration_[s]':float,
                                'YYYY-MM-DD':str,
                                'MJDstart':float,
                                'MJDend':float,
                                'priband':str
                                })
    pointings = pointings.to_numpy()
    for i in range(1, 1000):
        print(fov(pointings[i, -1], 25))
    gap_durations(pointings[:, 2], pointings[:, 4], pointings[:, 5])
    with open('durations.pkl', 'wb') as f:
        pickle.dump([pointings[:, 3], pointings[:, 2]], f)
    plt.hist(pointings[:, 2], bins = 100, facecolor = 'w', edgecolor = 'k', hatch = '/')
    plt.xlabel('Duration (s)')
    plt.ylabel('Number of pointings (s)')
    plt.show()

if(__name__ == '__main__'):
    cli()
