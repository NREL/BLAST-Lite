# Paul Gasper, NREL
# This model is fit to data reported by Bank et al from commercial NMC-LTO cells.
# https://doi.org/10.1016/j.jpowsour.2020.228566

import numpy as np
from functions.extract_stressors import extract_stressors
from functions.state_functions import update_power_state

# EXPERIMENTAL AGING DATA SUMMARY:
# Calendar aging varies temperature and SOC. There is almost no calendar aging impact
# at all until 80 Celsius.
# Cycle aging varies temperature, C-rate, and depth-of-discharge.

# MODEL SENSITIVITY
# The model predicts degradation rate versus time as a function of temperature and average
# state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as 
# a function of C-rate, temperature, and depth-of-discharge (DOD dependence is assumed to be linear, no aging data)

# MODEL LIMITATIONS
# Calendar aging has competition between capacity gain and capacity loss. There is an experimental
# case (80 Celsius, 5% SOC) that has complex behavior not modeled here.
# Astonishingly enough, the cycling degradation model is actually _overestimating_ capacity fade for most cases.
# The exception here is at very high temperature (60+ Celsius), where the fade is high, but not quite as high as observed degradation.

class Nmc_Lto_10Ah_Battery:

    def __init__(self, degradation_scalar=1):
        # States: Internal states of the battery model
        self.states = {
            'qLoss_t': np.array([0]),
            'qGain_t': np.array([0]),
            'qLoss_EFC': np.array([0]),
        }

        # Outputs: Battery properties derived from state values
        self.outputs = {
            'q': np.array([1]),
            'q_t_loss': np.array([1]),
            'q_t_gain': np.array([1]),
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
            'alpha': np.array([np.nan]),
            'beta': np.array([np.nan]),
            'gamma': np.array([np.nan]),
        }

        # Expermental range: details on the range of experimental conditions, i.e.,
        # the range we expect the model to be valid in
        self.experimental_range = {
            'cycling_temperature': [30, 60],
            'dod': [0, 1],
            'soc': [0, 1],
            'max_rate_charge': 10,
            'max_rate_discharge': 10,
        }

        # Degradation scalar - scales all state changes by a coefficient
        self._degradation_scalar = degradation_scalar

    # Nominal capacity
    @property
    def _cap(self):
        return 10.2

    # Define life model parameters
    @property
    def _params_life(self):
        return {
            # Capacity fade parameters
            'alpha_0': 3.11e+11,
            'alpha_1': -34.8,
            'alpha_2': 1.07,
            'alpha_p': 0.473,
            'beta_0': 7.86e+10,
            'beta_1': -35.8,
            'beta_2': 3.94,
            'beta_p': -0.553,
            'gamma_0': 1.29,
            'gamma_1': 7.83e-05,
            'gamma_2': 4.02,
            'gamma_3': -8.33,
            'gamma_p': 0.526,
        }
        
    # Battery model
    def update_battery_state(self, t_secs, soc, T_celsius):
        # Update the battery states, based both on the degradation state as well as the battery performance
        # at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.
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
        
        stressors = extract_stressors(t_secs, soc, T_celsius)
        # Unpack and store some stressors for debugging or plotting
        delta_t_days = stressors["delta_t_days"]
        delta_efc = stressors["delta_efc"]
        TdegK = stressors["TdegK"]
        soc = stressors["soc"]
        Ua = stressors["Ua"]
        dod = stressors["dod"]
        Crate = stressors["Crate"]
        t_days = self.stressors['t_days'][-1] + delta_t_days
        efc = self.stressors['efc'][-1] + delta_efc
        stressors_norm = np.array([delta_t_days, t_days, delta_efc, efc, np.mean(TdegK), np.mean(soc), np.mean(Ua), dod, Crate])
        for k, v in zip(self.stressors.keys(), stressors_norm):
            self.stressors[k] = np.append(self.stressors[k], v)
            
        self.__update_rates(stressors)
        self.__update_states(stressors)
        self.__update_outputs()

    def update_battery_state_repeating(self):
        # Update the battery states, based both on the degradation state as well as the battery performance
        # at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
        # degradation rates are unchanging for every timestep, and don't need to be calculated again
        
        # Just accumulate the cumulative stressors (total time, charge throughput in units of EFCs)
        self.stressors['t_days'] = np.append(self.stressors['t_days'], self.stressors['t_days'][-1] + self.stressors['delta_t_days'][-1])
        self.stressors['efc'] = np.append(self.stressors['efc'], self.stressors['efc'][-1] + self.stressors['delta_efc'][-1])

        # copy end values of old stressors into a dict
        stressors = {}
        for k, v in zip(self.stressors.keys(), self.stressors.values()):
            stressors[k] = v[-1]

        self.__update_states(stressors)
        self.__update_outputs()
    
    def __update_rates(self, stressors):
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
        TdegKN = TdegK / (273.15 + 45)

        # Grab parameters
        p = self._params_life

        # Calculate the degradation coefficients
        alpha = p['alpha_0'] * np.exp(p['alpha_1']/TdegKN) * np.exp(p['alpha_2']*soc/TdegKN)
        beta = p['beta_0'] * np.exp(p['beta_1']/TdegKN) * np.exp(p['beta_2']*soc/TdegKN)
        gamma = (
            (p['gamma_0'] + p['gamma_1']*Crate + p['gamma_2']*(dod**3))
            * np.exp(p['gamma_3']/TdegKN)
        )

        # Calculate time based average of each rate
        alpha = np.trapz(alpha, x=t_secs) / delta_t_secs
        beta = np.trapz(beta, x=t_secs) / delta_t_secs
        gamma = np.trapz(gamma, x=t_secs) / delta_t_secs

        # Store rates
        rates = np.array([alpha, beta, gamma])
        for k, v in zip(self.rates.keys(), rates):
            self.rates[k] = np.append(self.rates[k], v)

    def __update_states(self, stressors):
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
        dq_t_gain = self._degradation_scalar * update_power_state(states['qGain_t'][-1], delta_t_days, r['alpha'], p['alpha_p'])
        dq_t_loss = self._degradation_scalar * update_power_state(states['qLoss_t'][-1], delta_t_days, r['beta'], p['beta_p'])
        dq_EFC = self._degradation_scalar * update_power_state(states['qLoss_EFC'][-1], delta_efc, r['gamma'], p['gamma_p'])

        # Accumulate and store states
        dx = np.array([dq_t_loss, dq_t_gain, dq_EFC])
        for k, v in zip(states.keys(), dx):
            x = self.states[k][-1] + v
            self.states[k] = np.append(self.states[k], x)
    
    def __update_outputs(self):
        # Calculate outputs, based on current battery state
        states = self.states

        # Capacity
        q_t_loss = 1 - states['qLoss_t'][-1]
        q_t_gain = 1 + states['qGain_t'][-1]
        q_EFC = 1 - states['qLoss_EFC'][-1]
        q = 1 - states['qLoss_t'][-1] + states['qGain_t'][-1] - states['qLoss_EFC'][-1]

        # Assemble output
        out = np.array([q, q_t_loss, q_t_gain, q_EFC])
        # Store results
        for k, v in zip(list(self.outputs.keys()), out):
            self.outputs[k] = np.append(self.outputs[k], v)