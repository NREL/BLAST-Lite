# Paul Gasper, NREL
# This model is fit to LG MJ1 cell aging data reported as part of the EU EVERLASTING battery project, report D2.3
# https://everlasting-project.eu/wp-content/uploads/2020/03/EVERLASTING_D2.3_final_20200228.pdf
# Cell tests were reported in early 2020, so likely 2018 or 2019 LG MJ1 cells.
# High energy density 18650s but poor cycle life.

import numpy as np
from blast.models.degradation_model import BatteryDegradationModel


class Nmc811_GrSi_LGMJ1_4Ah_Battery(BatteryDegradationModel):
    """
    Model fit to LG MJ1 cell aging data reported as part of the EU EVERLASTING battery project, report D2.3
    https://everlasting-project.eu/wp-content/uploads/2020/03/EVERLASTING_D2.3_final_20200228.pdf
    Cell tests were reported in early 2020, so likely 2018 or 2019 LG MJ1 cells.
    High energy density 18650s but poor cycle life.

    .. note::
        EXPERIMENTAL AGING DATA SUMMARY:
            Calendar aging varied SOC (10%, 70%, 90%) and temperature.
            Cycle aging varied temperature and C-rates; all DOD is 80% (10%-90%). NO ACCELERATED FADE OBSERVED.
            Relative discharge capacity is reported from measurements recorded at 25 Celsius and C/20 rate.

        MODEL SENSITIVITY
            The model predicts degradation rate versus time as a function of temperature and average
            state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as 
            a function of C-rate, temperature, and depth-of-discharge (DOD dependence is assumed to be linear, no aging data)

        MODEL LIMITATIONS
            Cycle degradation predictions WILL NOT PREDICT KNEE-POINT due to limited data.
            OPERATION AT HIGH DOD PREDCTIONS ARE LIKELY INACCURATE (it is unclear what voltage window corresponds to SOCs defined in the test data).
            NMC811 is known to degrade quickly at voltages above 4.1 V.
    """

    def __init__(self, degradation_scalar: float = 1, label: str = "NMC811-GrSi LG MJ1"):
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
            'cycling_temperature': [0, 50],
            'dod': [0.2, 0.8],
            'soc': [0.1, 0.9],
            'max_rate_charge': 1,
            'max_rate_discharge': 3,
        }

        # Degradation scalar - scales all state changes by a coefficient
        self._degradation_scalar = degradation_scalar
        # Label for plotting
        self._label = label

    # Nominal capacity
    @property
    def _cap(self):
        return 3.5

    # Define life model parameters
    @property
    def _params_life(self):
        return {
            # Capacity fade parameters
            'qcal_A': 0.0353,
            'qcal_B': -1.03e+03,
            'qcal_C': 57.7,
            'qcal_p': 0.743,
            'qcyc_A': 1.77e-07,
            'qcyc_B': 8.08e-13,
            'qcyc_C': 2.21e-07,
            'qcyc_D': 2.25e+03,
            'qcyc_E': 1.14e+04,
            'qcyc_p': 0.695,
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
        dq_t = self._degradation_scalar * self.update_power_state(states['qLoss_t'][-1], delta_t_days, r['k_cal'], p['qcal_p'])
        dq_EFC = self._degradation_scalar * self.update_power_state(states['qLoss_EFC'][-1], delta_efc, r['k_cyc'], p['qcyc_p'])

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