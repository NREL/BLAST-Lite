# Paul Gasper, NREL
# This model is fit to SECOND LIFE data on Nissan Leaf half-modules (2p cells) by Braco et al.
# https://doi.org/10.1109/EEEIC/ICPSEUROPE54979.2022.9854784 (calendar aging data)
# https://doi.org/10.1016/j.est.2020.101695 (cycle aging data)
# Note that these cells are already hugely degraded, starting out at an average relative capacity
# of 70%. So the model reports q and qNew, where qNew is relative to initial

import numpy as np
from functions.extract_stressors import extract_stressors
from functions.state_functions import update_power_state

# EXPERIMENTAL AGING DATA SUMMARY:
# Calendar aging widely varied SOC and temperature.
# Cycle aging is only at a single condition (25 Celsius, 100% DOD, 1C-1C).

# MODEL SENSITIVITY
# The model predicts degradation rate versus time as a function of temperature and average
# state-of-charge and degradation rate is only a function equivalent full cycles.

# MODEL LIMITATIONS
# Cycling degradation IS ONLY A FUNCTION OF CHARGE THROUGHPUT due to limited aging data.
# Cycling degradation predictions ARE ONLY VALID NEAR 25 CELSIUS, 100% DOD, 1 C CHARGE/DISCHARGE RATE.

class Lmo_Gr_NissanLeaf66Ah_2ndLife_Battery:

    def __init__(self):
        # States: Internal states of the battery model
        self.states = {
            'qLoss_t': np.array([0]),
            'qLoss_EFC': np.array([0]),
        }

        # Outputs: Battery properties derived from state values
        self.outputs = {
            'q': np.array([1]),
            'q_t': np.array([1]),
            'q_EFC': np.array([1]),
            'qNew': np.array([0.7]),
        }

        # Stressors: History of stressors on the battery
        self.stressors = {
            'delta_t_days': np.array([np.nan]), 
            't_days': np.array([0]),
            'delta_efc': np.array([np.nan]), 
            'efc': np.array([0]),
            'TdegK': np.array([np.nan]),
            'soc': np.array([np.nan]),
        }

        # Rates: History of stressor-dependent degradation rates
        self.rates = {
            'k_cal': np.array([np.nan]),
        }

    # Nominal capacity
    @property
    def _cap_2ndLife(self):
        return 46
    
    @property
    def _cap(self):
        return 66

    # Define life model parameters
    @property
    def _params_life(self):
        return {
            # Capacity fade parameters
            'qcal_A': 3.25e+08,
            'qcal_B': -7.58e+03,
            'qcal_C': 162,
            'qcal_p': 0.464,
            'qcyc_A': 7.58e-05,
            'qcyc_p': 1.08,
        }
        
    # Battery model
    def update_battery_state(self, t_secs, soc, T_celsius):
        # Update the battery states, based both on the degradation state as well as the battery performance
        # at the ambient temperature, T_celsius
        # Inputs:
            #   t_secs (ndarry): vector of the time in seconds since beginning of life for the soc_timeseries data points
            #   soc (ndarry): vector of the state-of-charge of the battery at each t_sec
            #   T_celsius (ndarray): the temperature of the battery during this time period, in Celsius units.
            
        # Check some input types:
        if not isinstance(t_secs, np.ndarray):
            raise TypeError('Input "t_secs" must be a numpy.ndarray')
        if not isinstance(soc, np.ndarray):
            raise TypeError('Input "soc" must be a numpy.ndarray')
        if not isinstance(T_celsius, np.ndarray):
            raise TypeError('Input "T_celsius" must be a numpy.ndarray')
        if not (len(t_secs) == len(soc) and len(t_secs) == len(T_celsius)):
            raise ValueError('All input timeseries must be the same length')

        self.__update_states(t_secs, soc, T_celsius)
        self.__update_outputs()
    
    def __update_states(self, t_secs, soc, T_celsius):
        # Update the battery states, based both on the degradation state as well as the battery performance
        # at the ambient temperature, T_celsius
        # Inputs:
            #   t_secs (ndarry): vector of the time in seconds since beginning of life for the soc_timeseries data points
            #   soc (ndarry): vector of the state-of-charge of the battery at each t_sec
            #   T_celsius (ndarray): the temperature of the battery during this time period, in Celsius units.
            
        # Extract stressors
        delta_t_secs = t_secs[-1] - t_secs[0]
        delta_t_days, delta_efc, TdegK, soc, Ua, dod, Crate, cycles = extract_stressors(t_secs, soc, T_celsius)

        # Grab parameters
        p = self._params_life

        # Calculate the degradation coefficients
        k_cal = p['qcal_A'] * np.exp(p['qcal_B']/TdegK) * np.exp(p['qcal_C']*soc/TdegK)

        # Calculate time based average of each rate
        k_cal = np.trapz(k_cal, x=t_secs) / delta_t_secs

        # Calculate incremental state changes
        states = self.states
        # Capacity
        dq_t = update_power_state(states['qLoss_t'][-1], delta_t_days, k_cal, p['qcal_p'])
        dq_EFC = update_power_state(states['qLoss_EFC'][-1], delta_efc, p['qcyc_A'], p['qcyc_p'])

        # Accumulate and store states
        dx = np.array([dq_t, dq_EFC])
        for k, v in zip(states.keys(), dx):
            x = self.states[k][-1] + v
            self.states[k] = np.append(self.states[k], x)

        # Store stressors
        t_days = self.stressors['t_days'][-1] + delta_t_days
        efc = self.stressors['efc'][-1] + delta_efc
        stressors = np.array([delta_t_days, t_days, delta_efc, efc, np.mean(TdegK), np.mean(soc)])
        for k, v in zip(self.stressors.keys(), stressors):
            self.stressors[k] = np.append(self.stressors[k], v)

        # Store rates
        rates = np.array([k_cal])
        for k, v in zip(self.rates.keys(), rates):
            self.rates[k] = np.append(self.rates[k], v)
    
    def __update_outputs(self):
        # Calculate outputs, based on current battery state
        states = self.states

        # Capacity
        q_t = 1 - states['qLoss_t'][-1]
        q_EFC = 1 - states['qLoss_EFC'][-1]
        q = 1 - states['qLoss_t'][-1] - states['qLoss_EFC'][-1]
        qNew = 0.7 * q

        # Assemble output
        out = np.array([q, q_t, q_EFC, qNew])
        # Store results
        for k, v in zip(list(self.outputs.keys()), out):
            self.outputs[k] = np.append(self.outputs[k], v)