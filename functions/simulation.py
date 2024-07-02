import numpy as np
from numpy import interp
from .rainflow import extract_cycles

def simulate_battery_life(battery, input_timeseries: dict, simulation_years: float = None,
                            is_constant_input: bool = False, breakpoints_max_time_diff_s: float = 86400, breakpoints_max_EFC_diff: float = 1):
    # run battery life simulation over the input, or repeat for the number of years specified
    #   check input timeseries
    #       needs Temperature_C, SOC, t_secs keys
    #       values all need to be the same length
    #   if is_constant_input
    #       run life sim repeating until simulation is longer than simulation_years
    #   else
    #       check if we need to tile
    #           check that the input is periodic (SOC and temperatures don't change absurdly)
    #           tile input timeseries until it is 'simulation_years' long
    #       step through (maybe) tiled input_timeseries

    # Check if required keys are in the input_timeseries dictionary
    required_keys = ['Temperature_C', 'SOC', 'Time_s']
    for key in required_keys:
        if key not in input_timeseries:
            raise ValueError(f"Required key '{key}' is missing in the input_timeseries dictionary.")
    # Check if all values are numpy arrays and of the same length
    lengths = [len(input_timeseries[key]) for key in required_keys]
    if len(set(lengths)) != 1:
        raise ValueError("All numpy arrays in input_timeseries must have the same length.")
    
    if is_constant_input:
        # Only calculate stressors / degradation rates on the first timestep to dramatically accelerate the simulation
        battery.update_battery_state(input_timeseries['Time_s'], input_timeseries['SOC'], input_timeseries['Temperature_C'])
        years_simulated = battery.stressors['t_days'][-1]/365
        while years_simulated < simulation_years:
            battery.update_battery_state_repeating()
            years_simulated += battery.stressors['t_days'][-1]/365
        return battery
    else:
        # Unpack the inputs, calculate equivalent full cycles (EFCs)
        t_secs = input_timeseries['Time_s']
        soc = input_timeseries['SOC']
        temperature = input_timeseries['Temperature_C']

        # Check if we need to tile the inputs, do it if we do
        if simulation_years is not None and t_secs[-1] < simulation_years*365*24*3600:
            # Check that inputs are periodic before tiling
            if soc[-1] - soc[0] > 0.2:
                raise UserWarning(f"Inputs are being tiled to simulate {simulation_years} years, big jumps between repeats may lead to unrealistic simulations. Consider using the 'make_inputs_periodic' helper function.")
            if temperature[-1] - temperature[0] > 10:
                raise UserWarning(f"Inputs are being tiled to simulate {simulation_years} years, big jumps between repeats may lead to unrealistic simulations. Consider using the 'make_inputs_periodic' helper function.")
            
            # Tile the inputs. When tiling time, assume the timestep between repeats is the most common timestep.
            n_tile = np.ceil((simulation_years*365*24*3600)/t_secs[-1]).astype(int)
            delta_t_secs = np.ediff1d(t_secs, to_begin=0)
            delta_t_secs[0] = np.median(delta_t_secs)
            delta_t_secs = np.tile(delta_t_secs, n_tile)
            delta_t_secs[0] = 0
            t_secs = np.cumsum(delta_t_secs)
            soc = np.tile(soc, n_tile)
            temperature = np.tile(temperature, n_tile)

            # Cutoff once simulation is at simulation_years
            idx_end = np.argwhere(t_secs > simulation_years*365*24*3600)
            idx_end = idx_end[0][0]
            t_secs = t_secs[:idx_end]
            soc = soc[:idx_end]
            temperature = temperature[:idx_end]

        # Chunk the input timeseries into cycles (SOC turning points) using the rainflow algorithm
        cycles_generator = (
                [span, mean, count, idx_start, idx_end] for span, mean, count, idx_start, idx_end in extract_cycles(soc)
            )
        idx_turning_point = []
        for cycle in cycles_generator:
            idx_turning_point.append(cycle[-1])
        idx_turning_point = np.array(idx_turning_point)

        # Low DOD cycles can be extremely short, or the cycle time during long periods of storage excessively long.
        # Find breakpoints for simulating degradation, either simulating degradation when the number of EFCs between
        # turning points is greater than 1, or every day if EFCs accumlate slower than that.        
        EFCs = np.cumsum(np.abs(np.ediff1d(soc, to_begin=0)))/2
        simulation_breakpoints = find_breakpoints(t_secs, EFCs, idx_turning_point, max_time_diff_s=breakpoints_max_time_diff_s, max_EFC_diff=breakpoints_max_EFC_diff)
        # Add the final point manually if it's not there by coinkydink
        if simulation_breakpoints[-1] != len(t_secs)-1:
            simulation_breakpoints.append(len(t_secs)-1)

        # Simulate battery life between breakpoints until we reach the end
        prior_breakpoint = 0
        for breakpoint in simulation_breakpoints:
            battery.update_battery_state(t_secs[prior_breakpoint:breakpoint], soc[prior_breakpoint:breakpoint], temperature[prior_breakpoint:breakpoint])
            prior_breakpoint = breakpoint
        return battery
        
def find_breakpoints(time_seconds, EFCs, idx_turning_point, max_time_diff_s = 86400, max_EFC_diff = 1):
    # Find breakpoints for simulating battery degradation, where degradation is calculated once
    # a certain number of equivalent full cycles has passed, otherwise after a certain time has passed.
    # Defaults are to calculate degradation every 1 EFC or every day.
    breakpoints = []
    last_breakpoint = 0
    for i in idx_turning_point:
        # Check if change in time_seconds greater than max_time_diff_s
        if time_seconds[i] - time_seconds[last_breakpoint] > max_time_diff_s:
            # Add a breakpoint at each day
            for i_time in range(last_breakpoint, i):
                if time_seconds[i_time] - time_seconds[last_breakpoint] > max_time_diff_s:
                    breakpoints.append(i_time)
                    last_breakpoint = i_time
            continue

        # Check change in EFCs greater than max_EFC_diff
        if EFCs[i] - EFCs[last_breakpoint] > max_EFC_diff:
            breakpoints.append(i)
            last_breakpoint = i

    return breakpoints

def make_inputs_periodic(input_timeseries: dict, interp_time_window_hours: float):
    # linear interpolate the last 'interp_time_window_hours' of the input to the first value to ensure periodicity
    #   check input timeseries
    #       needs Temperature_C, SOC, Time_s keys
    #       values of each key all need to be numpy arrays of the same length
    #   interpolate the last 'interp_time_window_hours' of SOC and Temperature_C values from the point at the start of 'interp_time_window_hours' to the first value of each series using the t_secs array, which is time in seconds

    # Check if required keys are in the input_timeseries dictionary
    required_keys = ['Temperature_C', 'SOC', 'Time_s']
    for key in required_keys:
        if key not in input_timeseries:
            raise ValueError(f"Required key '{key}' is missing in the input_timeseries dictionary.")
    # Check if all values are numpy arrays and of the same length
    lengths = [len(input_timeseries[key]) for key in required_keys]
    if len(set(lengths)) != 1:
        raise ValueError("All numpy arrays in input_timeseries must have the same length.")

    # Extract arrays from input_timeseries
    t_secs = input_timeseries['Time_s']
    temperature = input_timeseries['Temperature_C']
    soc = input_timeseries['SOC']

    # Convert interp_time_window_hours to seconds
    interp_time_window_secs = interp_time_window_hours * 3600.0

    # Identify the range of indices to interpolate
    end_time = t_secs[-1]
    start_time = end_time - interp_time_window_secs
    idx = np.argmin(np.abs(t_secs - start_time))

    # Interpolate SOC and Temperature_C
    soc_interp = interp(t_secs[idx:], t_secs[idx:], [soc[idx], soc[0]])
    temperature_interp = interp(t_secs[idx:], t_secs[idx:], [temperature[idx], temperature[0]])
    # Replace the last part of soc and temperature arrays with interpolated values
    soc[idx:] = soc_interp
    temperature[idx:] = temperature_interp

    # Update the input_timeseries dictionary with interpolated values
    input_timeseries['SOC'] = soc
    input_timeseries['Temperature_C'] = temperature

    return input_timeseries

"""
# Test of 'find_breakpoints'
from .rainflow import extract_cycles

hours = 6 * 24
t_hours = np.arange(hours + 1)
t_secs = t_hours * 3600
soc = np.tile(
    np.array(
        [1, 0.6, 0.4, 0.4, 0.6, 1,
        1, 0.5, 0, 0.5, 1, 1,
        1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1,
        0.98, 1, 0.98, 1, 0.98, 1,
        0.98, 1, 0.98, 1, 0.98, 1,
        0.98, 1, 0.98, 1, 0.98, 1,
        0.98, 1, 0.98, 1, 0.98, 1,
        0.98, 1, 0.98, 1, 0.98, 1,
        0.98, 1, 0.98, 1, 0.98, 1,
        0.98, 1, 0.98, 1, 0.98, 1,
        0.98, 1, 0.98, 1, 0.98, 1,]
        ),
    2
)
soc = np.append(soc, 1)
EFCs = np.cumsum(np.abs(np.ediff1d(soc, to_begin=0)))/2

# Chunk the input timeseries into cycles (SOC turning points) using the rainflow algorithm
cycles_generator = (
        [span, mean, count, idx_start, idx_end] for span, mean, count, idx_start, idx_end in extract_cycles(soc)
    )
idx_turning_point = []
for cycle in cycles_generator:
    idx_turning_point.append(cycle[-1])
idx_turning_point = np.array(idx_turning_point)

# Find the breakpoints for battery life simulation
idx_simulation_breakpoints = find_breakpoints(t_secs, EFCs, idx_turning_point)

plt.plot(t_hours/24, soc)
plt.plot(t_hours[idx_turning_point]/24, soc[idx_turning_point], '.g')
plt.plot(t_hours[idx_simulation_breakpoints]/24, soc[idx_simulation_breakpoints], 'or')
"""