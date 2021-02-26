import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import time
import pickle
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt

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

if(__name__ == '__main__'):
   
    print(fov('3GHz', 25))

    pd.options.display.max_colwidth = 2000 # Surprising this is necessary...
    pointings = pd.read_csv('tmp_coords.clean',
                            usecols = [0,1,2,3,9,10],
                            delimiter = '  ',
                            dtype={'RA_[deg]':float,
                                'Dec_[deg]':float,
                                'Duration_[s]':float,
                                'YYYY-MM-DD':str,
                                'project_code':str,
                                'priband':str
                                })
    pointings = pointings.to_numpy()
    print("Pointings loaded")

    vlass_count = 0
    other_count = 0
    vlass_durations = []
    other_durations = []

    for i in range(pointings.shape[0]):
        # Date range:
        date = datetime.strptime(pointings[i, 3], '%Y-%m-%d')
        start_date = datetime(2020, 10, 23)
        end_date = datetime(2020, 11, 9)
        if(date < start_date):
            #print('Date too old')
            continue
        elif(date > end_date):
            #print('Date too new')
            continue
        # Check project_code
        if('VLASS' not in pointings[i, 4]):
            #print(pointings[i, 4])
            other_count += 1
            other_durations.append(pointings[i, 2])
        else:
            vlass_count += 1
            vlass_durations.append(pointings[i, 2])
 
 # print('OTHER:')
   # print(other_count)
   # print(np.histogram(other_durations, bins = 100))
   # print(np.sum(other_durations))
   # print(np.histogram(other_durations, bins = [0, 30, 1000000]))
   # plt.hist(other_durations, bins = 100)

   # plt.title('Other')
   # plt.show()
   # print('VLASS:')
   # print(vlass_count)
   # print(np.histogram(vlass_durations, bins = 100))
   # print(np.sum(vlass_durations))
   # print(np.histogram(vlass_durations, bins = [0, 30, 100000000]))
   # plt.hist(vlass_durations, bins = 100)
   # plt.title('VLASS')
   # plt.show()
