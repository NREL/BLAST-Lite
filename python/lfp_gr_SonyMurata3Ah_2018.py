# Paul Gasper, NREL
import numpy as np
from functions.extract_stressors import extract_stressors
from functions.state_functions import update_power_B_state, update_sigmoid_state
import scipy.stats as stats

# EXPERIMENTAL AGING DATA SUMMARY:
# Aging test matrix varied temperature and state-of-charge for calendar aging, and
# varied depth-of-discharge, average state-of-charge, and C-rates for cycle aging.
# There is NO LOW TEMPERATURE cycling aging data, i.e., no lithium-plating induced by
# kinetic limitations on cell performance; CYCLING WAS ONLY DONE AT 25 CELSIUS AND 45 CELSIUS,
# so any model predictions at low temperature cannot incorporate low temperature degradation modes.
# Discharge capacity

# MODEL SENSITIVITY
# The model predicts degradation rate versus time as a function of temperature and average
# state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as 
# a function of average state-of-charge during a cycle, depth-of-discharge, and average of the
# charge and discharge C-rates.

# MODEL LIMITATIONS
# There is no influence of TEMPERATURE on CYCLING DEGRADATION RATE due to limited data. This is
# NOT PHYSICALLY REALISTIC AND IS BASED ON LIMITED DATA.

class Lfp_Gr_SonyMurata3Ah_Battery:
    # Model predicting the degradation of Sony-Murata 3 Ah LFP-Gr cylindrical cells.
    # Data is from Technical University of Munich, reported in studies led by Maik Naumann.
    # Capacity model identification was conducted at NREL. Resistance model is from Naumann et al.
    # Naumann et al used an interative fitting procedure, but it was found that lower model error could be
    # achieved by simply reoptimizing all resistance growth parameters to the entire data set.
    # Calendar aging data source: https://doi.org/10.1016/j.est.2018.01.019
    # Cycle aging data source: https://doi.org/10.1016/j.jpowsour.2019.227666
    # Model identification source: https://doi.org/10.1149/1945-7111/ac86a8
    # Degradation rate is a function of the aging stressors, i.e., ambient temperature and use.
    # The state of the battery is updated throughout the lifetime of the cell.
    # Performance metrics are capacity and DC resistance. These metrics change as a function of the
    # cell's current degradation state, as well as the ambient temperature. The model predicts time and 
    # cycling dependent degradation. Cycling dependent degradation includes a break-in mechanism as well
    # as long term cycling fade; the break-in mechanism strongly influenced results of the accelerated
    # aging test, but is not expected to have much influence on real-world applications.
    # Parameters to modify to change fade rates:
    #   q1_b0: rate of capacity loss due to calendar degradation
    #   q5_b0: rate of capacity loss due to cycling degradation
    #   k_ref_r_cal: rate of resistance growth due to calendar degradation
    #   A_r_cyc: rate of resistance growth due to cycling degradation

    def __init__(self):
        # States: Internal states of the battery model
        self.states = {
            'qLoss_LLI_t': np.array([0]),
            'qLoss_LLI_EFC': np.array([0]),
            'qLoss_BreakIn_EFC': np.array([1e-10]),
            'rGain_LLI_t': np.array([0]),
            'rGain_LLI_EFC': np.array([0]),
        }

        # Outputs: Battery properties derived from state values
        self.outputs = {
            'q': np.array([1]),
            'q_LLI_t': np.array([1]),
            'q_LLI_EFC': np.array([1]),
            'q_BreakIn_EFC': np.array([1]),
            'r': np.array([1]),
            'r_LLI_t': np.array([1]),
            'r_LLI_EFC': np.array([1]),
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
            'q1': np.array([np.nan]),
            'q3': np.array([np.nan]),
            'q5': np.array([np.nan]),
            'q7': np.array([np.nan]),
            'r_kcal': np.array([np.nan]),
            'r_kcyc': np.array([np.nan]),
        }

    # Nominal capacity
    @property
    def _cap(self):
        return 3

    # Define life model parameters
    @property
    def _params_life(self):
        return {
            # Capacity fade parameters
            'q2': 0.000130510034211874,
            'q1_b0': 0.989687151293590,             # CHANGE to modify calendar degradation rate
            'q1_b1': -2881067.56019324,
            'q1_b2': 8742.06309157261,
            'q3_b0': 0.000332850281062177,          
            'q3_b1': 734553185711.369,
            'q3_b2': -2.82161575620780e-06,
            'q3_b3': -3284991315.45121,
            'q3_b4': 0.00127227593657290,
            'q8': 0.00303553871631028,
            'q9': 1.43752162947637,
            'q7_b0': 0.582258029148225,
            'q7_soc_skew': 0.0583128906965484,
            'q7_soc_width': 0.208738181522897,
            'q7_dod_skew': -3.80744333129564,
            'q7_dod_width': 1.16126260428210,
            'q7_dod_growth': 25.4130804598602,
            'q6': 1.12847759334355,
            'q5_b0': -6.81260579372875e-06,         # CHANGE to modify cycling degradation rate
            'q5_b1': 2.59615973160844e-05,
            'q5_b2': 2.11559710307295e-06,
            # Resistance growth parameters
            'k_ref_r_cal': 3.4194e-10,    # CHANGE to modify calendar degradation rate
            'Ea_r_cal': 71827,
            'C_r_cal': -3.3903,
            'D_r_cal': 1.5604,
            'A_r_cyc': -0.002,        # CHANGE to modify cycling degradation rate
            'B_r_cyc': 0.0021,
            'C_r_cyc': 6.8477,
            'D_r_cyc': 0.91882
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
        q1 = np.abs(
        p['q1_b0']
        * np.exp(p['q1_b1']*(1/(TdegK**2))*(Ua**0.5))
        * np.exp(p['q1_b2']*(1/TdegK)*(Ua**0.5))
        )
        q3 = np.abs(
            p['q3_b0']
            * np.exp(p['q3_b1']*(1/(TdegK**4))*(Ua**(1/3)))
            * np.exp(p['q3_b2']*(TdegK**3)*(Ua**(1/4)))
            * np.exp(p['q3_b3']*(1/(TdegK**3))*(Ua**(1/3)))
            * np.exp(p['q3_b4']*(TdegK**2)*(Ua**(1/4)))
            )
        q5 = np.abs(
            p['q5_b0']
            + p['q5_b1']*dod
            + p['q5_b2']*np.exp((dod**2)*(Crate**3))
            )
        q7 = np.abs(
            p['q7_b0']
            * skewnormpdf(soc, p['q7_soc_skew'], p['q7_soc_width'])
            * skewnormpdf(dod, p['q7_dod_skew'], p['q7_dod_width'])
            * sigmoid(dod, 1, p['q7_dod_growth'], 1)
            )
        k_temp_r_cal = (
            p['k_ref_r_cal']
            * np.exp((-p['Ea_r_cal'] / 8.3144) * (1/TdegK - 1/298.15))
        )
        k_soc_r_cal = p['C_r_cal'] * (soc - 0.5)**3 + p['D_r_cal']
        k_Crate_r_cyc = p['A_r_cyc'] * Crate + p['B_r_cyc']
        k_dod_r_cyc = p['C_r_cyc']* (dod - 0.5)**3 + p['D_r_cyc']

        # Calculate time based average of each rate
        q1 = np.trapz(q1, x=t_secs) / delta_t_secs
        q3 = np.trapz(q3, x=t_secs) / delta_t_secs
        #q5 = np.trapz(q5, x=t_secs) / delta_t_secs # no time varying inputs
        q7 = np.trapz(q7, x=t_secs) / delta_t_secs # no time varying inputs
        k_temp_r_cal = np.trapz(k_temp_r_cal, x=t_secs) / delta_t_secs
        k_soc_r_cal = np.trapz(k_soc_r_cal, x=t_secs) / delta_t_secs # no time varying inputs
        #k_Crate_r_cyc = np.trapz(k_Crate_r_cyc, x=t_secs) / delta_t_secs # no time varying inputs
        #k_dod_r_cyc = np.trapz(k_dod_r_cyc, x=t_secs) / delta_t_secs # no time varying inputs
 
        # Calculate incremental state changes
        states = self.states
        # Capacity
        dq_LLI_t = update_sigmoid_state(states['qLoss_LLI_t'][-1], delta_t_days, q1, p['q2'], q3)
        dq_LLI_EFC = update_power_B_state(states['qLoss_LLI_EFC'][-1], delta_efc, q5, p['q6'])
        if delta_efc / delta_t_days > 2: # only evalaute if more than 2 full cycles per day
            dq_BreakIn_EFC = update_sigmoid_state(states['qLoss_BreakIn_EFC'][-1], delta_efc, q7, p['q8'], p['q9'])
        else:
            dq_BreakIn_EFC = 0

        # Resistance
        dr_LLI_t = k_temp_r_cal * k_soc_r_cal * delta_t_secs
        dr_LLI_EFC = k_Crate_r_cyc * k_dod_r_cyc * delta_efc / 100

        # Accumulate and store states
        dx = np.array([dq_LLI_t, dq_LLI_EFC, dq_BreakIn_EFC, dr_LLI_t, dr_LLI_EFC])
        for k, v in zip(states.keys(), dx):
            x = self.states[k][-1] + v
            self.states[k] = np.append(self.states[k], x)
        
        # Store stressors
        t_days = self.stressors['t_days'][-1] + delta_t_days
        efc = self.stressors['efc'][-1] + delta_efc
        stressors = np.array([delta_t_days, t_days, delta_efc, efc, np.mean(TdegK), np.mean(soc), np.mean(Ua), dod, Crate])
        for k, v in zip(self.stressors.keys(), stressors):
            self.stressors[k] = np.append(self.stressors[k], v)

        # Store rates
        rates = np.array([q1, q3, q5, q7, k_temp_r_cal * k_soc_r_cal, k_Crate_r_cyc * k_dod_r_cyc])
        for k, v in zip(self.rates.keys(), rates):
            self.rates[k] = np.append(self.rates[k], v)
    
    def __update_outputs(self):
        # Calculate outputs, based on current battery state
        states = self.states
        p = self._params_life

        # Capacity
        q_LLI_t = 1 - states['qLoss_LLI_t'][-1]
        q_LLI_EFC = 1 - states['qLoss_LLI_EFC'][-1]
        q_BreakIn_EFC = 1 - states['qLoss_BreakIn_EFC'][-1]
        q = 1 - states['qLoss_LLI_t'][-1] - states['qLoss_LLI_EFC'][-1] - states['qLoss_BreakIn_EFC'][-1]
        
        # Resistance
        r_LLI_t = 1 + states['rGain_LLI_t'][-1]
        r_LLI_EFC = 1 + states['rGain_LLI_EFC'][-1]
        r = 1 + states['rGain_LLI_t'][-1] + states['rGain_LLI_EFC'][-1]

        # Assemble output
        out = np.array([q, q_LLI_t, q_LLI_EFC, q_BreakIn_EFC, r, r_LLI_t, r_LLI_EFC])
        # Store results
        for k, v in zip(list(self.outputs.keys()), out):
            self.outputs[k] = np.append(self.outputs[k], v)

def sigmoid(x, alpha, beta, gamma):
    return 2*alpha*(1/2 - 1/(1 + np.exp((beta*x)**gamma)))

def skewnormpdf(x, skew, width):
    x_prime = (x-0.5)/width
    return 2 * stats.norm.pdf(x_prime) * stats.norm.cdf(skew * (x_prime))