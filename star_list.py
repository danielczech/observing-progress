import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
import time
import pickle
import argparse
import sys
import pandas as pd

def region_tf(points, ra_1, dec_1):
    """Transform to make the local region flat
    """
    # Points into SkyCoords
    sc_points = SkyCoord(points[:, 1], points[:, 2], unit='deg')
    new_origin = SkyCoord(ra_1, dec_1, unit='deg')
    # Transform points so locally flat
    new_frame = new_origin.skyoffset_frame()
    sc_points_tf = sc_points.transform_to(new_frame)
    sc_dec_tf = sc_points_tf.lat.value
    sc_ra_tf = sc_points_tf.lon.value
    return sc_ra_tf, sc_dec_tf

def gen_star_lists(pointing, catalogue, radius):
    """Calculate the stars available as drawn from the catalogue
    that fall within the field of view surrounding the pointing. 
    """
    # Transform catalogue:
    tf_ra, tf_dec = region_tf(catalogue, pointing[0], pointing[1])
    # Extract region:
    idxs = np.where(tf_ra**2 + tf_dec**2 < radius**2)[0]
    return idxs

def cli(args = sys.argv[0]):
    usage = "{} [options]".format(args)
    description = """Get index list of stars per FoV."""
    parser = argparse.ArgumentParser(prog = 'star-lists', 
        usage = usage, description = description) 
    help_info = """Specify catalogue, starting position in file 
                and number of lines."""
    parser.add_argument('catalogue', type = str, help = help_info)
    parser.add_argument('pointings', type = str, help = help_info)
    parser.add_argument('lstart', type = int, help = help_info)
    parser.add_argument('n', type = int, help = help_info)
    if(len(sys.argv[1:]) == 0): 
        print('Missing args')
        parser.exit()
    args = parser.parse_args()
    main(catalogue = args.catalogue, pointings = args.pointings, 
            lstart = args.lstart, n = args.n)

def main(catalogue, pointings, lstart, n):

    RADIUS = 0.004284 # Estimate: 1.02*0.21/25.0/2.0
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
    print("Pointings loaded")

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

    pointing_list = []
    # Step through each pointing from pointing lstart to pointing n:
    for i in range(lstart, lstart + n):
        start = time.time()
        idx_i = gen_star_lists(pointings[i, 0:2], db_full, RADIUS)
        if(len(idx_i) > 0):
            # field 1: pointing index
            # field 2: duration
            # field 3: date
            # field 4: indices of stars observed in main catalogue
            star_list = [i, pointings[i, 2], pointings[i, 3], idx_i]
            pointing_list.append(star_list)
        print('{} of {} for {} stars in {:.3} s'.format(i - lstart + 1, 
            n, len(idx_i), time.time() - start))
    print(pointing_list)
    print("Saving source lists: ")
    with open('stars_per_pointing_{}.pkl'.format(RADIUS), 'wb') as f:
        pickle.dump(pointing_list, f)

if(__name__=="__main__"):
    cli()
