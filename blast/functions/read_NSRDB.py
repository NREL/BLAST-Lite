# Extract time-series data for a single site from NSRDB, format for battery life
import h5pyd
import pandas as pd
from scipy.spatial import cKDTree
from geopy.geocoders import Nominatim
import numpy as np

def get_nsrdb_temperature_data(location):
    # Get latitude and longitude from the location string
    # Get temperature data from the nearest coordinates in the NSRDB

    def nearest_site(tree, lat_coord, lon_coord):
        lat_lon = np.array([lat_coord, lon_coord])
        dist, pos = tree.query(lat_lon)  
        return pos
    
    loc = Nominatim(user_agent="Geopy Library")
    # entering the location name
    getLoc = loc.geocode(location)
    # printing address
    print("Found location: ", getLoc.address)
    # printing latitude and longitude
    lat = getLoc.latitude
    lon = getLoc.longitude
    coords = (lat, lon)

    # Open NSRDB .h5 file
    f = h5pyd.File("/nrel/nsrdb/v3/nsrdb_2018.h5", 'r')
    # Get spatial coordinates
    dset_coords = f['coordinates'][...]
    tree = cKDTree(dset_coords)
    pos_idx = nearest_site(tree, coords[0], coords[1])
    print("Input coordiantes: \t {}".format(coords))
    print("Coordinates of nearest point in NSRDB: \t {}".format(dset_coords[pos_idx])) 
    # Extract time_index and convert to datetime
    # NOTE: time_index is saved as byte-strings and must be decoded
    time_index = pd.to_datetime(f['time_index'][...].astype(str))
    # Initialize DataFrame to store time-series data
    time_series = pd.DataFrame(index=time_index)
    # Extract variables needed
    for var in ['air_temperature']:
        # Get dataset
        ds = f[var]
        # Extract scale factor
        scale_factor = ds.attrs['psm_scale_factor']
        # Extract site 100 and add to DataFrame
        time_series[var] = ds[:, pos_idx] / scale_factor
    
    time_series = time_series.reset_index()
    time_series['dt'] = time_series['index'] - time_series['index'][0]
    time_series["Time_s"] = time_series['dt'].dt.seconds
    time_series = time_series.drop(columns=["dt","index"])

    return time_series