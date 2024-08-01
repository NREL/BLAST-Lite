# Paul Gasper, NREL
# This model is fit to Panasonic 18650B NCA-Gr cells. High-ish energy density 18650 cells with adequate lifetime.
# Calendar data is reported by Keil et al (https://dx.doi.org/10.1149/2.0411609jes)
# Cycling data is reported by Preger et al (https://doi.org/10.1149/1945-7111/abae37) and
# is available at batteryarchive.com.
# I'm not aware of any study conducting both calendar aging and cycle aging of these cells.

import numpy as np
from blast.models.degradation_model import BatteryDegradationModel


class Nca_Gr_Panasonic3Ah_Battery(BatteryDegradationModel):

    """
    Model fit to Panasonic 18650B NCA-Gr cells. High-ish energy density 18650 cells with
    adequate lifetime.
    Calendar data is reported by Keil et al (https://dx.doi.org/10.1149/2.0411609jes)
    Cycling data is reported by Preger et al (https://doi.org/10.1149/1945-7111/abae37) and
    is available at batteryarchive.com.
    The authors of BLAST-Lite are not aware of any study conducting both calendar aging and
    cycle aging of these cells.

    .. note::
        EXPERIMENTAL AGING DATA SUMMARY:
            Calendar aging widely varied SOC at 25, 40, and 50 Celsius. 300 days max.
            Cycle aging varied temperature and C-rates, and DOD. Some accelerating fade is observed
            at room temperature and high DODs but isn't modeled well here. That's not a huge problem,
            because the modeled lifetime is quite short anyways.

        MODEL SENSITIVITY
            The model predicts degradation rate versus time as a function of temperature and average
            state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as 
            a function of C-rate, temperature, and depth-of-discharge (DOD dependence is assumed to be linear, no aging data)

        MODEL LIMITATIONS
            Cycle degradation predictions WILL NOT PREDICT KNEE-POINT due to limited data.
            Cycle aging is only modeled at 25, 35, and 45 Celsius, PREDICTIONS OUTSIDE THIS 
            TEMPERATURE RANGE MAY BE OPTIMISTIC.
    """

    def __init__(self, degradation_scalar: float = 1, label: str = "NCA-Gr Panasonic"):
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
            'dod': np.array([np.nan]),
            'Crate': np.array([np.nan]),
        }

        # Rates: History of stressor-dependent degradation rates
        self.rates = {
            'k_cal': np.array([np.nan]),
            'k_cyc': np.array([np.nan]),
        }

        # Expermental range: details on the range of experimental conditions, i.e.,
        # the range we expect the model to be valid in
        self.experimental_range = {
            'cycling_temperature': [15, 35],
            'dod': [0.8, 1],
            'soc': [0, 1],
            'max_rate_charge': 0.5,
            'max_rate_discharge': 2,
        }

        # Degradation scalar - scales all state changes by a coefficient
        self._degradation_scalar = degradation_scalar
        # Label for plotting
        self._label = label

    # Nominal capacity
    @property
    def cap(self):
        return 3.2

    # Define life model parameters
    @property
    def _params_life(self):
        return {
            # Capacity fade parameters
            'qcal_A': 75.4,
            'qcal_B': -3.34e+03,
            'qcal_C': 353,
            'qcal_p': 0.512,
            'qcyc_A': 1.86e-06,
            'qcyc_B': 4.74e-11,
            'qcyc_C': 0.000177,
            'qcyc_D': 3.34e-11,
            'qcyc_E': 2.81e-09,
            'qcyc_p': 0.699,
        }
        
    # Battery model
    def update_rates(self, stressors):
        # Calculate and update battery degradation rates based on stressor values
        # Inputs:
        #   stressors (dict): output from extract_stressors

        # Unpack stressors
        t_secs = stressors["t_secs"]
        delta_t_secs = t_secs[-1] - t_secs[0]
        delta_t_days = stressors["delta_t_days"]
        delta_efc = stressors["delta_efc"]
        TdegK = stressors["TdegK"]
        soc = stressors["soc"]
        Ua = stressors["Ua"]
        dod = stressors["dod"]
        Crate = stressors["Crate"]

        # Grab parameters
        p = self._params_life

        # Calculate the degradation coefficients
        k_cal = p['qcal_A'] * np.exp(p['qcal_B']/TdegK) * np.exp(p['qcal_C']*soc/TdegK)
        k_cyc  = (
            (p['qcyc_A'] + p['qcyc_B']*Crate + p['qcyc_C']*dod)
            * (np.exp(p['qcyc_D']/TdegK) + np.exp(-p['qcyc_E']/TdegK))
        )

        # Calculate time based average of each rate
        k_cal = np.trapz(k_cal, x=t_secs) / delta_t_secs
        k_cyc = np.trapz(k_cyc, x=t_secs) / delta_t_secs

        # Store rates
        rates = np.array([k_cal, k_cyc])
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
        dq_t = self._degradation_scalar * self._update_power_state(states['qLoss_t'][-1], delta_t_days, r['k_cal'], p['qcal_p'])
        dq_EFC = self._degradation_scalar * self._update_power_state(states['qLoss_EFC'][-1], delta_efc, r['k_cyc'], p['qcyc_p'])

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