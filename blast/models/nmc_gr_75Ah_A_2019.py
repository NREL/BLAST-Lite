# Paul Gasper, NREL
# Large format pouch NMC-Gr cell from a large manufacturer, with ~75 Ah capacity and an
# energy-to-power ratio of 16 h^-1 (relatively high power, suitable for fast charging).
# Experimental test data is reported in https://doi.org/10.1016/j.est.2023.109042.
# This is for the 'NMC-Gr A' cell reported in the paper.

import numpy as np
from blast.models.degradation_model import BatteryDegradationModel

# EXPERIMENTAL AGING DATA SUMMARY:
# Experimental test data is reported in https://doi.org/10.1016/j.est.2023.109042.
# Aging test matrix varied temperature and state-of-charge for calendar aging, and
# varied depth-of-discharge, average state-of-charge, and C-rates for cycle aging.

# MODEL SENSITIVITY
# The model predicts degradation rate versus time as a function of temperature and average
# state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as 
# a function of average state-of-charge during a cycle, depth-of-discharge, and average of the
# charge and discharge C-rates.

# MODEL LIMITATIONS
# Charging and discharging rates were conservative, following cycling protocols as suggested by the
# cell manufacturer. Degradation at higher charging or discharging rates will not be simulated accurately.


class NMC_Gr_75Ah_A(BatteryDegradationModel):
    # Model predicting degradation of 'NMC-Gr A' cells from https://doi.org/10.1016/j.est.2023.109042.

    def __init__(self, degradation_scalar=1, label="NMC-Gr A 75Ah"):
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
        }

        # Stressors: History of stressors on the battery
        self.stressors = {
            'delta_t_days': np.array([np.nan]), 
            't_days': np.array([0]),
            'delta_efc': np.array([np.nan]), 
            'efc': np.array([0]),
            'TdegK': np.array([np.nan]),
            'soc': np.array([np.nan]), 
            'Ua': np.array([np.nan]), 
            'dod': np.array([np.nan]), 
            'Crate': np.array([np.nan]),
        }

        # Rates: History of stressor-dependent degradation rates
        self.rates = {
            'kcal': np.array([np.nan]),
            'kcyc': np.array([np.nan]),
        }

        # Expermental range: details on the range of experimental conditions, i.e.,
        # the range we expect the model to be valid in
        self.experimental_range = {
            'cycling_temperature': [10, 45],
            'dod': [0.8, 1],
            'soc': [0, 1],
            'max_rate_charge': 2,
            'max_rate_discharge': 2,
        }

        # Degradation scalar - scales all state changes by a coefficient
        self._degradation_scalar = degradation_scalar
        # Label for plotting
        self._label = label

    # Nominal capacity
    @property
    def _cap(self):
        return 75

    # Define life model parameters
    @property
    def _params_life(self):
        return {
            # Capacity fade parameters
            'p1': 8.71e-07,
            'p2': 10,
            'pcal': 0.578,
            'p3': 0.0267,
            'p4': 1.64,
            'p5': 0.011,
            'p6': -0.44,
            'pcyc': 1.11,
        }
    
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
        TdegKN = TdegK / (293.15 + 35)
        UaN = Ua / 0.123 # Ua at about 50% SOC

        # Grab parameters
        p = self._params_life

        # Calculate the degradation coefficients
        kcal = (p['p1']
            * np.exp(p['p2']*(TdegKN**3) * (1/(UaN**(1/3)))))
        kcyc = np.abs(p['p3'] + p['p4']*(TdegKN**3)*(dod**0.5) 
        + p['p5']*np.exp((1/(TdegKN**0.5)) * (Crate**2))
        + p['p6']*np.exp((dod**0.5) * (TdegKN**2) * (Crate**0.5)))
        
        # Calculate time based average of each rate
        kcal = np.trapz(kcal, x=t_secs) / delta_t_secs
        kcyc = np.trapz(kcyc, x=t_secs) / delta_t_secs

        # Store rates
        rates = np.array([kcal, kcyc])
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
        dq_t = self._degradation_scalar * self.update_power_state(states['qLoss_t'][-1], delta_t_days/1e3, r['kcal'], p['pcal'])
        dq_EFC = self._degradation_scalar * self.update_power_state(states['qLoss_EFC'][-1], delta_efc/1e4, r['kcyc'], p['pcyc'])

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

        # Assemble output
        out = np.array([q, q_t, q_EFC])
        # Store results
        for k, v in zip(list(self.outputs.keys()), out):
            self.outputs[k] = np.append(self.outputs[k], v)
