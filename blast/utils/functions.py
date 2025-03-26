"""Helper functions for fetching and formatting simulation input data."""
import h5pyd
import pandas as pd
from scipy.spatial import cKDTree
from geopy.geocoders import Nominatim
import numpy as np
from numpy import interp
import importlib
from blast.models._available_models import available_models
from blast.utils.rainflow import reversals

def simulate_all_models(*args, **kwargs):
    models = available_models()
    cells = []
    for model in models:
        This_Battery_Model = getattr(importlib.import_module("blast.models"), model)
        cell = This_Battery_Model()
        cell.simulate_battery_life(*args, **kwargs)
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

    loc = Nominatim(user_agent="BLAST-Lite")
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

# TODO (Paul 2/3/2025): something is not exactly perfect but I think it's close enough
def scale_vehicle_profile_to_annual_efcs(profile, desired_efcs_per_year, show_efcs=False):
    profile['dSOC'] = profile['SOC'].diff().fillna(0)
    efcs = 0.5 * profile['dSOC'].abs().sum()
    duration = profile['Time_s'].iloc[-1] - profile['Time_s'].iloc[0]
    efcs_per_year = efcs / (duration / (24*3600*365))
    # Add high throughput/dod events on to the end of the profile until the desired efcs per year is exceeded
    while desired_efcs_per_year > efcs_per_year:
        # Make a high dod event
        high_dod = pd.DataFrame()
        high_dod['Time_s'] = np.arange(0, 4*3600, 300)
        high_dod['SOC'] = np.linspace(0.95, 0.1, int(4*3600/300))
        high_dod['dSOC'] = high_dod['SOC'].diff().fillna(0)
        high_dod['Time_s'] = high_dod['Time_s'] + profile['Time_s'].iloc[-1] + 300
        high_dod2 = pd.DataFrame()
        high_dod2['Time_s'] = np.arange(0, 3600, 300)
        high_dod2['SOC'] = np.linspace(0.1, 0.7, int(3600/300))
        high_dod2['dSOC'] = high_dod2['SOC'].diff().fillna(0)
        high_dod2['Time_s'] = high_dod2['Time_s'] + high_dod['Time_s'].iloc[-1] + 300
        high_dod3 = pd.DataFrame()
        high_dod3['Time_s'] = np.arange(0, 2*3600, 300)
        high_dod3['SOC'] = np.linspace(0.7, 0.4, int(2*3600/300))
        high_dod3['dSOC'] = high_dod3['SOC'].diff().fillna(0)
        high_dod3['Time_s'] = high_dod3['Time_s'] + high_dod2['Time_s'].iloc[-1] + 300
        high_dod4 = pd.DataFrame()
        high_dod4['Time_s'] = np.arange(0, 4*3600, 300)
        high_dod4['SOC'] = np.linspace(0.4, 0.95, int(4*3600/300))
        high_dod4['dSOC'] = high_dod4['SOC'].diff().fillna(0)
        high_dod4['Time_s'] = high_dod4['Time_s'] + high_dod3['Time_s'].iloc[-1] + 300
        # concat the high DOD event to the end of the profile
        profile = pd.concat([profile, high_dod, high_dod2, high_dod3, high_dod4], axis=0).reset_index(drop=True)
        # recalculate efcs
        profile['dSOC'] = profile['SOC'].diff().fillna(0)
        efcs = 0.5 * profile['dSOC'].abs().sum()  
        duration = profile['Time_s'].iloc[-1] - profile['Time_s'].iloc[0]
        efcs_per_year = efcs / (duration / (24*3600*365))
        print('added a high throughput event')
    
    # Identify reversals
    soc = profile['SOC'].to_numpy()
    idx_reversals = []
    for reverse in reversals(soc):
        idx_reversals.append(reverse[0])
    
    # Calculate the total amount of rest time that would need to be added to the profile to make the actual EFCs per year equal to the desired EFCs per year
    total_rest_time = duration * ((efcs_per_year / desired_efcs_per_year) - 1)
    rest_time_per_interval = np.ceil(total_rest_time / len(idx_reversals))
    # print(rest_time_per_interval)

    # Add a rest at constant SOC (equal to the SOC at the start of the rest) to the end of each interval
    new_row_count = 0
    for i in range(len(idx_reversals)):
        # print(idx_reversals[i], new_row_count, idx_reversals[i] + new_row_count, profile.shape[0])
        if i == 0:
            before_rest = profile.iloc[:1].copy()
        else:
            before_rest = profile.iloc[:idx_reversals[i] + new_row_count].copy()
        new_rows_here = int(rest_time_per_interval / 300)
        rest_period = pd.DataFrame({
            'Time_s': [profile.loc[idx_reversals[i] + new_row_count, 'Time_s'] + j * 300 for j in range(1, new_rows_here + 1)],
            'SOC': profile.loc[idx_reversals[i] + new_row_count, 'SOC'],
            'dSOC': 0
        })
        after_rest = profile.iloc[idx_reversals[i] + new_row_count:]
        new_time = profile.loc[idx_reversals[i] + new_row_count:, 'Time_s'] + rest_time_per_interval
        # cast new_time to same dtype as 'Time_s' column
        new_time = new_time.astype(profile['Time_s'].dtype)
        after_rest.loc[:, 'Time_s'] = new_time
        
        # print(before_rest['Time_s'].to_numpy()[-1], rest_period['Time_s'].to_numpy()[0], rest_period['Time_s'].to_numpy()[-1], after_rest['Time_s'].to_numpy()[0])
        # print(before_rest['SOC'].to_numpy()[-1], rest_period['SOC'].to_numpy()[0], rest_period['SOC'].to_numpy()[-1], after_rest['SOC'].to_numpy()[0])

        new_row_count += new_rows_here
        profile = pd.concat([before_rest, rest_period, after_rest], axis=0).reset_index(drop=True)
        # plt.plot(profile['Time_s'] / (24*3600), profile['SOC'])
    
    profile['dSOC'] = profile['SOC'].diff().fillna(0)
    if show_efcs:
        efcs = 0.5 * profile['dSOC'].abs().sum()  
        duration = profile['Time_s'].iloc[-1] - profile['Time_s'].iloc[0]
        efcs_per_year = efcs / (duration / (24*3600*365))
        print(efcs_per_year)
    return profile

def decimate_and_rescale_profile(profile, decimation_factor, tol=1e-1, show_efcs=False):
    dSOC = profile['SOC'].diff().fillna(0)
    efcs_original = 0.5 * dSOC.abs().sum()
    if show_efcs:
        print(f"Original EFCs: {efcs_original}")
    profile = profile.iloc[::decimation_factor, :].reset_index(drop=True)
    dSOC = profile['SOC'].diff().fillna(0)
    efcs = 0.5 * dSOC.abs().sum()
    while np.abs(efcs_original - efcs) > tol:
        dSOC = dSOC * (efcs_original / efcs)
        profile['SOC'] = np.cumsum(dSOC) + profile['SOC'].iloc[0]
        # Ceiling SOC at 1 and floor at 0
        profile['SOC'] = np.maximum(0, np.minimum(1, profile['SOC']))
        dSOC = profile['SOC'].diff().fillna(0)
        efcs = 0.5 * dSOC.abs().sum()
    if show_efcs:
        print(f"Decimated and rescaled EFCs: {efcs}")
    return profile

def derate_profile(profile, derating_factor, max_soc=0.9):
    # Derate a profile by scaling dSOC by the derating factor.
    # If feasible, reduce max_soc as well.
    dSOC = profile['SOC'].diff().fillna(0)
    dSOC = dSOC * (1-derating_factor)
    profile['SOC'] = np.cumsum(dSOC) + profile['SOC'].iloc[0]
    if profile['SOC'].min() > (1 - max_soc) and profile['SOC'].max() > max_soc:
        profile['SOC'] = profile['SOC'] - (profile['SOC'].max() - max_soc)
    return profile

def rescale_profile(profile, rescaling_factor):
    # Rescale a profile by scaling dSOC by the rescaling factor.
    dSOC = profile['SOC'].diff().fillna(0)
    dSOC = dSOC * rescaling_factor
    profile['SOC'] = np.cumsum(dSOC) + profile['SOC'].iloc[0]
    return profile

def rescale_soc(soc, rescaling_factor):
    # Rescale an soc vector by scaling dSOC by the rescaling factor.
    dSOC = np.diff(soc, prepend=soc[0])
    dSOC = dSOC * rescaling_factor
    soc = np.cumsum(dSOC) + soc[0]
    if np.max(soc) > 1 or np.min(soc) < 0:
        # print('Warning: rescaled SOC is outside of the range [0, 1], modified SOC will not directly reflect inputs.')
        soc = np.maximum(0, np.minimum(1, soc)) 
    return soc

def simulate_battery_life_distribution(Battery, input, degradation_scalar_std=0.15, degradation_scalar_mean=1., **kwargs):
    # Simulate the distribution of battery life from normally distributed degradation scalars.
    # Simulations are done at -3, -2, -1, 0, 1, 2, 3 sigma. The resulting outputs from the 
    # simulations at any non-zero sigma are resized to be the same size as the nominal simulation
    # to allow for easy plotting.
    sigmas = [0, -3, -2, -1, 1, 2, 3]
    batteries_sigma = {}
    for sigma in sigmas:
        degradation_scalar = degradation_scalar_mean + sigma*degradation_scalar_std
        batt = Battery(degradation_scalar=degradation_scalar)
        batt.simulate_battery_life(input, **kwargs)
        if sigma == 0:
            output_len = len(batt.stressors['t_days'])
        else:
            outputs = batt.outputs
            stressors = batt.stressors
            for var in outputs:
                outputs[var] = np.interp(np.linspace(0, len(outputs[var]), output_len), np.arange(len(outputs[var])), outputs[var])
            for var in stressors:
                stressors[var] = np.interp(np.linspace(0, len(stressors[var]), output_len), np.arange(len(stressors[var])), stressors[var])
            batt.outputs = outputs
            batt.stressors = stressors

        batteries_sigma[sigma] = batt
    
    return batteries_sigma
    