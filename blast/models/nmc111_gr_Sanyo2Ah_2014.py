# Paul Gasper, NREL
# This model is replicated as reported by Schmalsteig et al, J. Power Sources 257 (2014) 325-334
# http://dx.doi.org/10.1016/j.jpowsour.2014.02.012

import numpy as np
from blast.utils.state_functions import update_power_state
from blast.models.degradation_model import BatteryDegradationModel

# EXPERIMENTAL AGING DATA SUMMARY:
# Calendar aging varied SOC at 50 Celsius, and temperature at 50% state-of-charge.
# Cycle aging varied depth-of-discharge and average state-of-charge at 35 Celsius at
# charge and discharge rates of 1C. 
# Relative discharge capacity is reported from measurements recorded at 35 Celsius and 1C rate.
# Relative DC resistance is reported after fitting of 10s 1C discharge pulses near 50% state-of-charge.

# MODEL SENSITIVITY
# The model predicts degradation rate versus time as a function of temperature and average
# state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as 
# a function of average voltage and depth-of-discharge.

# MODEL LIMITATIONS
# Cycle degradation predictions are NOT SENSITIVE TO TEMPERATURE OR C-RATE. Cycling degradation predictions
# are ONLY ACCURATE NEAR 1C RATE AND 35 CELSIUS CELL TEMPERATURE. 

class Nmc111_Gr_Sanyo2Ah_Battery(BatteryDegradationModel):
    # Model predicting the degradation of Sanyo UR18650E cells, published by Schmalsteig et al:
    # http://dx.doi.org/10.1016/j.jpowsour.2014.02.012.
    # More detailed analysis of cell performance and voltage vs. state-of-charge data was copied from
    # Ecker et al: http://dx.doi.org/10.1016/j.jpowsour.2013.09.143 (KNEE POINTS OBSERVED IN ECKER ET AL
    # AT HIGH DEPTH OF DISCHARGE WERE SIMPLY NOT ADDRESSED DURING MODEL FITTING BY SCHMALSTEIG ET AL).
    # Voltage lookup table here use data from Ecker et al for 0/10% SOC, and other values were extracted
    # from Figure 1 in Schmalsteig et al using WebPlotDigitizer.

    def __init__(self, degradation_scalar=1, label="NMC111-Gr Sanyo"):
        # States: Internal states of the battery model
        self.states = {
            'qLoss_t': np.array([0]),
            'qLoss_EFC': np.array([0]),
            'rGain_t': np.array([0]),
            'rGain_EFC': np.array([0]),
        }

        # Outputs: Battery properties derived from state values
        self.outputs = {
            'q': np.array([1]),
            'q_t': np.array([1]),
            'q_EFC': np.array([1]),
            'r': np.array([1]),
            'r_t': np.array([1]),
            'r_EFC': np.array([1]),
        }

        # Stressors: History of stressors on the battery
        self.stressors = {
            'delta_t_days': np.array([np.nan]), 
            't_days': np.array([0]),
            'delta_efc': np.array([np.nan]), 
            'efc': np.array([0]),
            'TdegK': np.array([np.nan]),
            'soc': np.array([np.nan]),
            'Vrms': np.array([np.nan]),
            'dod': np.array([np.nan]),
        }

        # Rates: History of stressor-dependent degradation rates
        self.rates = {
            'alpha_cap': np.array([np.nan]),
            'beta_cap': np.array([np.nan]),
            'alpha_res': np.array([np.nan]),
            'beta_res': np.array([np.nan]),
        }

        # Expermental range: details on the range of experimental conditions, i.e.,
        # the range we expect the model to be valid in
        self.experimental_range = {
            'cycling_temperature': [20, 40],
            'dod': [0, 1],
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
    def _cap(self):
        return 2.15

    # SOC index
    @property
    def _soc_index(self):
        return np.array([0,0.008637153,0.026779514,0.044921875,0.063064236,0.081206597,0.099348958,0.117491319,0.135633681,0.153776042,0.171918403,0.190060764,0.208203125,0.226345486,0.244487847,0.262630208,0.280772569,0.298914931,0.317057292,0.335199653,0.353342014,0.371484375,0.389626736,0.407769097,0.425911458,0.444053819,0.462196181,0.480338542,0.498480903,0.516623264,0.534765625,0.552907986,0.571050347,0.589192708,0.607335069,0.625477431,0.643619792,0.661762153,0.679904514,0.698046875,0.716189236,0.734331597,0.752473958,0.770616319,0.788758681,0.806901042,0.825043403,0.843185764,0.861328125,0.879470486,0.897612847,0.915755208,0.933897569,0.952039931,0.970182292,0.988324653,0.998220486,1])
    
    # OCV
    @property
    def _ocv(self):
        return np.array([3.331,3.345014187,3.37917149,3.411603677,3.440585632,3.466289865,3.490268982,3.511315401,3.529946658,3.547197821,3.561688798,3.574972194,3.586357962,3.597053683,3.605506753,3.613442288,3.620342753,3.626380661,3.632073544,3.63690387,3.642079219,3.646909545,3.652084894,3.657605266,3.663470662,3.670198615,3.677444104,3.686759732,3.696420384,3.708323686,3.720572012,3.734372943,3.749553967,3.765252525,3.781123596,3.797857224,3.814590852,3.832532062,3.849783226,3.866861877,3.884458064,3.900846669,3.917752809,3.934831461,3.953462717,3.971921462,3.991415276,4.011771649,4.03195551,4.05196686,4.070770628,4.087849279,4.104237884,4.120108955,4.135980025,4.152541142,4.160649189,4.162])

    # Voltage lookup table
    def calc_voltage(self, soc):
        # calculate cell voltage from soc vector
        return np.interp(soc, self._soc_index, self._ocv, left=self._ocv[0], right=self._ocv[-1])

    # Define life model parameters
    @property
    def _params_life(self):
        return {
            # Capacity fade parameters
            'qcal_A_V': 7.543,
            'qcal_B_V': -23.75,
            'qcal_C_TdegK': -6976,
            'qcal_p': 0.75,
            'qcyc_A_V': 7.348e-3,
            'qcyc_B_V': 3.667,
            'qcyc_C': 7.6e-4,
            'qcyc_D_DOD': 4.081e-3,
            'qcyc_p': 0.5,

            # Resistance growth parameters
            'rcal_A_V': 5.270,
            'rcal_B_V': -16.32,
            'rcal_C_TdegK': -5986,
            'rcal_p': 0.75,
            'rcyc_A_V': 2.153e-4,
            'rcyc_B_V': 3.725,
            'rcyc_C': -1.521e-5,
            'rcyc_D_DOD': 2.798e-4,
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
        # Calculate RMS voltage, charge throughput
        V_rms = np.sqrt(np.mean(self.calc_voltage(soc)**2))

        # Grab parameters
        p = self._params_life

        # Calculate the degradation coefficients
        alpha_cap = (p['qcal_A_V'] * self.calc_voltage(soc) + p['qcal_B_V']) * 1e6 * np.exp(p['qcal_C_TdegK'] / TdegK)
        alpha_res = (p['rcal_A_V'] * self.calc_voltage(soc) + p['rcal_B_V']) * 1e5 * np.exp(p['rcal_C_TdegK'] / TdegK)
        beta_cap  = (
            p['qcyc_A_V'] * (V_rms - p['qcyc_B_V']) ** 2
            + p['qcyc_C']
            + p['qcyc_D_DOD'] * dod
        )
        beta_res  = (
            p['rcyc_A_V'] * (V_rms - p['rcyc_B_V']) ** 2
            + p['rcyc_C']
            + p['rcyc_D_DOD'] * dod
        )

        # Calculate time based average of each rate
        alpha_cap = np.trapz(alpha_cap, x=t_secs) / delta_t_secs
        alpha_res = np.trapz(alpha_res, x=t_secs) / delta_t_secs

        # Store rates
        rates = np.array([alpha_cap, beta_cap, alpha_res, beta_res])
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
        Ah_throughput = delta_efc * 2 * self._cap

        # Grab parameters
        p = self._params_life

        # Grab rates, only keep most recent value
        r = self.rates.copy()
        for k, v in zip(r.keys(), r.values()):
            r[k] = v[-1]

        # Calculate incremental state changes
        states = self.states
        # Capacity
        dq_t = self._degradation_scalar * update_power_state(states['qLoss_t'][-1], delta_t_days, r['alpha_cap'], p['qcal_p'])
        dq_EFC = self._degradation_scalar * update_power_state(states['qLoss_EFC'][-1], Ah_throughput, r['beta_cap'], p['qcyc_p'])
        # Resistance
        dr_t = self._degradation_scalar * update_power_state(states['rGain_t'][-1], delta_t_days, r['alpha_res'], p['rcal_p'])
        dr_EFC = self._degradation_scalar * r['beta_res'] * Ah_throughput

        # Accumulate and store states
        dx = np.array([dq_t, dq_EFC, dr_t, dr_EFC])
        for k, v in zip(states.keys(), dx):
            x = self.states[k][-1] + v
            self.states[k] = np.append(self.states[k], x)
    
    def update_outputs(self, stressors):
        # Calculate outputs, based on current battery state
        states = self.states
        p = self._params_life

        # Capacity
        q_t = 1 - states['qLoss_t'][-1]
        q_EFC = 1 - states['qLoss_EFC'][-1]
        q = 1 - states['qLoss_t'][-1] - states['qLoss_EFC'][-1]
        
        # Resistance
        r_t = 1 + states['rGain_t'][-1]
        r_EFC = 1 + states['rGain_EFC'][-1]
        r = 1 + states['rGain_t'][-1] + states['rGain_EFC'][-1]

        # Assemble output
        out = np.array([q, q_t, q_EFC, r, r_t, r_EFC])
        # Store results
        for k, v in zip(list(self.outputs.keys()), out):
            self.outputs[k] = np.append(self.outputs[k], v)


