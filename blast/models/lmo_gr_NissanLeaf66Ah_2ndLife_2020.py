# Paul Gasper, NREL
# This model is fit to SECOND LIFE data on Nissan Leaf half-modules (2p cells) by Braco et al.
# https://doi.org/10.1109/EEEIC/ICPSEUROPE54979.2022.9854784 (calendar aging data)
# https://doi.org/10.1016/j.est.2020.101695 (cycle aging data)
# Note that these cells are already hugely degraded, starting out at an average relative capacity
# of 70%. So the model reports q and qNew, where qNew is relative to initial

import numpy as np
from blast.utils.state_functions import update_power_state
from blast.models.degradation_model import BatteryDegradationModel

# EXPERIMENTAL AGING DATA SUMMARY:
# Calendar aging widely varied SOC and temperature.
# Cycle aging is only at a single condition (25 Celsius, 100% DOD, 1C-1C).

# MODEL SENSITIVITY
# The model predicts degradation rate versus time as a function of temperature and average
# state-of-charge and degradation rate is only a function equivalent full cycles.

# MODEL LIMITATIONS
# Cycling degradation IS ONLY A FUNCTION OF CHARGE THROUGHPUT due to limited aging data.
# Cycling degradation predictions ARE ONLY VALID NEAR 25 CELSIUS, 100% DOD, 1 C CHARGE/DISCHARGE RATE.

class Lmo_Gr_NissanLeaf66Ah_2ndLife_Battery(BatteryDegradationModel):

    def __init__(self, degradation_scalar=1, label="LMO-Gr Nissan Leaf"):
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

        # Expermental range: details on the range of experimental conditions, i.e.,
        # the range we expect the model to be valid in
        self.experimental_range = {
            'cycling_temperature': [20, 30],
            'dod': [0.8, 1],
            'soc': [0, 1],
            'max_rate_charge': 1,
            'max_rate_discharge': 1,
        }

        # Degradation scalar - scales all state changes by a coefficient
        self._degradation_scalar = degradation_scalar
        # Label for plotting
        self._label = label

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
    def update_rates(self, stressors):
        # Calculate and update battery degradation rates based on stressor values
        # Inputs:
        #   stressors (dict): output from extract_stressors

        # Unpack stressors
        t_secs = stressors["t_secs"]
        delta_t_secs = t_secs[-1] - t_secs[0]
        TdegK = stressors["TdegK"]
        soc = stressors["soc"]
        Ua = stressors["Ua"]
        dod = stressors["dod"]
        Crate = stressors["Crate"]

        # Grab parameters
        p = self._params_life

        # Calculate the degradation coefficients
        k_cal = p['qcal_A'] * np.exp(p['qcal_B']/TdegK) * np.exp(p['qcal_C']*soc/TdegK)

        # Calculate time based average of each rate
        k_cal = np.trapz(k_cal, x=t_secs) / delta_t_secs

        # Store rates
        rates = np.array([k_cal])
        for k, v in zip(self.rates.keys(), rates):
            self.rates[k] = np.append(self.rates[k], v)

    def update_states(self, stressors):
        # Update the battery states, based both on the degradation state as well as the battery performance
        # at the ambient temperature, T_celsius
        # Inputs:
            #   stressors (dict): output from extract_stressors
            
        # Unpack stressors
        delta_t_days = stressors["delta_t_days"]
        delta_efc = stressors["delta_efc"]
        
        # Grab parameters
        p = self._params_life

        # Grab rates, only keep most recent value
        r = self.rates.copy()
        for k, v in zip(r.keys(), r.values()):
            r[k] = v[-1]

        # Calculate incremental state changes
        states = self.states
        # Capacity
        dq_t = self._degradation_scalar * update_power_state(states['qLoss_t'][-1], delta_t_days, r['k_cal'], p['qcal_p'])
        dq_EFC = self._degradation_scalar * update_power_state(states['qLoss_EFC'][-1], delta_efc, p['qcyc_A'], p['qcyc_p'])

        # Accumulate and store states
        dx = np.array([dq_t, dq_EFC])
        for k, v in zip(states.keys(), dx):
            x = self.states[k][-1] + v
            self.states[k] = np.append(self.states[k], x)
    
    def update_outputs(self, stressors):
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