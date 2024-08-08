"""Helper functions for fetching and formatting simulation input data."""

import h5pyd
import pandas as pd
from scipy.spatial import cKDTree
from geopy.geocoders import Nominatim
import numpy as np
from numpy import interp
import importlib
from blast.models._available_models import available_models

def simulate_all_models(**kwargs):
    models = available_models()
    cells = []
    for model in models:
        This_Battery_Model = getattr(importlib.import_module("blast.models"), model)
        cell = This_Battery_Model()
        cell.simulate_battery_life(**kwargs)
        cells.append(cell)
    
    return cells

def get_nsrdb_temperature_data(location: str = "Honolulu, Hawaii") -> pd.DataFrame:
    """
    Get temperature time-series data from NSRDB at the nearest coordinates to the
    provided location, and format for battery life simulation.

    Args:
        location (str):      Descriptive location string.

    Returns:
        pd.DataFrame:  Temperature time series at the specified location.
    """

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
    f = h5pyd.File("/nrel/nsrdb/v3/nsrdb_2018.h5", "r")
    # Get spatial coordinates
    dset_coords = f["coordinates"][...]
    tree = cKDTree(dset_coords)
    pos_idx = nearest_site(tree, coords[0], coords[1])
    print("Input coordiantes: \t {}".format(coords))
    print("Coordinates of nearest point in NSRDB: \t {}".format(dset_coords[pos_idx]))
    # Extract time_index and convert to datetime
    # NOTE: time_index is saved as byte-strings and must be decoded
    time_index = pd.to_datetime(f["time_index"][...].astype(str))
    # Initialize DataFrame to store time-series data
    time_series = pd.DataFrame(index=time_index)
    # Extract variables needed
    for var in ["air_temperature"]:
        # Get dataset
        ds = f[var]
        # Extract scale factor
        scale_factor = ds.attrs["psm_scale_factor"]
        # Extract site 100 and add to DataFrame
        time_series[var] = ds[:, pos_idx] / scale_factor

    time_series = time_series.reset_index()
    time_series["dt"] = time_series["index"] - time_series["index"][0]
    time_series["Time_s"] = time_series["dt"].values / np.timedelta64(1, "s")
    time_series["Temperature_C"] = time_series["air_temperature"]
    time_series = time_series.drop(columns=["dt", "index", "air_temperature"])

    return time_series


def make_inputs_periodic(
    input_timeseries: dict, interp_time_window_hours: float
) -> pd.DataFrame:
    """
    Linearly interpolate the last 'interp_time_window_hours' of the input to ensure periodicity.

    Args:
        input_timeseries (dict):     Dictionary with keys ['Temperature_C', 'SOC', 'Time_s']
                                    and timeseries values.
        interp_time_window_hours (float):    Number of hours to interpolate, from end of timeseries

    Returns:
        dict:    Dictionary with keys ['Temperature_C', 'SOC', 'Time_s'] after interpolation.
    """

    # Check if required keys are in the input_timeseries dictionary
    required_keys = ["Temperature_C", "SOC", "Time_s"]
    for key in required_keys:
        if key not in input_timeseries:
            raise ValueError(
                f"Required key '{key}' is missing in the input_timeseries dictionary."
            )
    # Check if all values are numpy arrays and of the same length
    lengths = [len(input_timeseries[key]) for key in required_keys]
    if len(set(lengths)) != 1:
        raise ValueError(
            "All numpy arrays in input_timeseries must have the same length."
        )

    # Extract arrays from input_timeseries
    t_secs = input_timeseries["Time_s"]
    temperature = input_timeseries["Temperature_C"]
    soc = input_timeseries["SOC"]

    # Convert interp_time_window_hours to seconds
    interp_time_window_secs = interp_time_window_hours * 3600.0

    # Identify the range of indices to interpolate
    end_time = t_secs[-1]
    start_time = end_time - interp_time_window_secs
    idx = np.argmin(np.abs(t_secs - start_time))

    # Interpolate SOC and Temperature_C
    soc_interp = interp(t_secs[idx:], t_secs[idx:], [soc[idx], soc[0]])
    temperature_interp = interp(
        t_secs[idx:], t_secs[idx:], [temperature[idx], temperature[0]]
    )
    # Replace the last part of soc and temperature arrays with interpolated values
    soc[idx:] = soc_interp
    temperature[idx:] = temperature_interp

    # Update the input_timeseries dictionary with interpolated values
    input_timeseries["SOC"] = soc
    input_timeseries["Temperature_C"] = temperature

    return input_timeseries


def assemble_one_year_input(
    soc: pd.DataFrame, climate: pd.DataFrame, time_interval: int = 3600
) -> pd.DataFrame:
    """
    Construct SOC and climate data into a single dataframe,
    linked by timestamp.

    Args:
        soc (pd.DataFrame):      Dataframe with SOC profile and 'Time_s' timestamp
        climate (pd.DataFrame):  Dataframe with temperature data and 'Time_s' timestamp column
        time_interval (int):     Interval at which the final data is resampled to

    Return:
        pd.DataFrame:  Combined DataFrame with time, soc, and temperature columns.
    """
    # Check if 'Time_s' column exists in both data frames
    if "Time_s" not in soc.columns:
        raise ValueError("soc dataframe does not contain 'Time_s' column.")
    if "Time_s" not in climate.columns:
        raise ValueError("climate dataframe does not contain 'Time_s' column.")
    # Make each data frame one year long
    soc = tile_to_one_year(soc)
    climate = tile_to_one_year(climate)
    # Merge data frames on 'Time_s' column
    combined_df = pd.merge(soc, climate, on="Time_s", how="outer")
    # Sort values by 'Time_s'
    combined_df = combined_df.sort_values("Time_s")
    # Interpolate to fill missing values
    combined_df = combined_df.interpolate(method="linear")
    # Resample to at the input interval (default 1 hour) from 0 to 31104000 (one year)
    index = pd.Index(np.arange(0, 31104000, time_interval), name="Time_s")
    combined_df = combined_df.set_index("Time_s").reindex(index)

    return combined_df.reset_index()


def tile_to_one_year(df: pd.DataFrame) -> pd.DataFrame:
    """
    Take a dataframe with timestamp column called 'Time_s'
    and resample to 1 year.

    Args:
        df (pd.DataFrame):   Dataframe to cut off to 1 year.

    Returns:
        pd.DataFrame:   Dataframe after cutting off to 1 year.
    """
    t_secs = df["Time_s"].values
    if t_secs[-1] < 365 * 24 * 3600:
        df = df.drop(columns="Time_s")
        df_values = df.to_numpy()
        # Tile the inputs. When tiling time, assume the timestep
        # between repeats is the most common timestep.
        n_tile = np.ceil((365 * 24 * 3600) / t_secs[-1]).astype(int)
        delta_t_secs = np.ediff1d(t_secs, to_begin=0)
        delta_t_secs[0] = np.median(delta_t_secs)
        delta_t_secs = np.tile(delta_t_secs, n_tile)
        delta_t_secs[0] = 0
        t_secs = np.cumsum(delta_t_secs)
        df_values = np.tile(df_values, (n_tile, 1))

        # Cutoff once simulation is at simulation_years
        idx_end = np.argwhere(t_secs > 365 * 24 * 3600)
        idx_end = idx_end[0][0]
        t_secs = t_secs[:idx_end]
        df_values = df_values[:idx_end]

        # remake into dataframe
        df = pd.DataFrame(df_values, columns=df.columns)
        df["Time_s"] = t_secs

    return df
