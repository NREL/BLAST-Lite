# Paul Gasper, NREL
# Model predicting the degradation of Sony-Murata US18650VTC5A 3.5 Ah NCA-GrSi cylindrical cells.
# Relatively high power cells, with 1.4 wt% Si in the graphite-si composite negative electrode.
# (Maximum continuous charge rate of 2C and max continuous discharge rate of 10C, even at 5 degC)
# Data is from Technical University of Munich, reported in studies led by Leo Wildfeuer and Alexander Karger.
# Accelerated aging test data reported in https://doi.org/10.1016/j.jpowsour.2022.232498
# The model here was identifed using AI-Batt at NREL. Note that the TUM authors have put more effort
# into model identification on this data set, see the following papers for more detailed models identified on this data set:
#   Mechanistic cycle aging model (LLI, LAM_PE, LAM_(NE,Gr), LAM_(NE,Si)): https://doi.org/10.1016/j.jpowsour.2023.233947
#   Mechanistic calendar aging model (considers impact of capacity check frequency): https://doi.org/10.1016/j.jpowsour.2023.233208

import numpy as np
from blast.utils.state_functions import update_sigmoid_state
from blast.models.degradation_model import BatteryDegradationModel

# EXPERIMENTAL AGING DATA SUMMARY:
# Aging test matrix varied temperature and state-of-charge for calendar aging, and
# varied depth-of-discharge, average state-of-charge, and C-rates for cycle aging.
# For calendar aging, in addition to temperature and SOC, capacity-check frequency was also varied. 

# MODEL SENSITIVITY
# The model predicts degradation rate versus time as a function of temperature and average
# state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as 
# a function of average state-of-charge during a cycle, depth-of-discharge, and average of the
# charge and discharge C-rates.

# MODEL LIMITATIONS
# I did not model the impact of capacity check frequency, only using the 6 week capacity check data.
# See https://doi.org/10.1016/j.jpowsour.2023.233208 for a detailed consideration of the capacity check frequency impact on aging.
# Only one cycling cell in the data shows 'knee-over' behavior, making empirical model identification of this behavior challenging.
# I simply neglected the knee over; this knee occurs at ~80% capacity fade, so note that predictions of degradation at maximum
# DOD and below 80% capacity should not be believed.

class NCA_GrSi_SonyMurata2p5Ah_Battery(BatteryDegradationModel):

    def __init__(self, degradation_scalar=1, label="NCA-GrSi Sony-Murata"):
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
            'soc_low': np.array([np.nan]),
            'soc_high': np.array([np.nan]),
            'Ua': np.array([np.nan]), 
            'dod': np.array([np.nan]), 
            'Cdis': np.array([np.nan]),
            'Cchg': np.array([np.nan]),
        }

        # Rates: History of stressor-dependent degradation rates
        self.rates = {
            'kcal': np.array([np.nan]),
            'pcal': np.array([np.nan]),
            'kcyc': np.array([np.nan]),
            'pcyc': np.array([np.nan]),
        }

        # Expermental range: details on the range of experimental conditions, i.e.,
        # the range we expect the model to be valid in
        self.experimental_range = {
            'cycling_temperature': [5, 35],
            'dod': [0.2, 1],
            'soc': [0, 1],
            'max_rate_charge': 2,
            'max_rate_discharge': 10,
        }

        # Degradation scalar - scales all state changes by a coefficient
        self._degradation_scalar = degradation_scalar
        # Label for plotting
        self._label = label

    # Nominal capacity
    @property
    def _cap(self):
        return 2.5

    # Define life model parameters
    @property
    def _params_life(self):
        return {
            # Capacity fade parameters
            'q2': 3.41,
            'q1_b0': 0.00748,
            'q1_b1': 8.9,
            'q1_b2': -4.87,
            'p_cal_b0': 1.01,
            'p_cal_b1': -0.133,
            'p_cal_b2': 1.22,
            'p_cal_b3': -0.56,
            'q4': 0.479,
            'q3_b0': 5.49,
            'q3_b1': -64.9,
            'q3_b2': 21.3,
            'q3_b3': 0.000801,
            'p_cyc_b0': 0.422,
            'p_cyc_b1': 2.64,
            'p_cyc_b2': 0.00627,
            'p_cyc_b3': -0.0342,
            'p_cyc_b4': 1.29,
            'p_cyc_b5': -2.14
        }

    # Override super class update_battery_state because this model uses a few unusual stressors
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
        
        stressors = self.extract_stressors(t_secs, soc, T_celsius)
        # Unpack and store some stressors for debugging or plotting
        delta_t_days = stressors["delta_t_days"]
        delta_efc = stressors["delta_efc"]
        TdegK = stressors["TdegK"]
        soc = stressors["soc"]
        Ua = stressors["Ua"]
        dod = stressors["dod"]
        soc_low = stressors["soc_low"]
        soc_high = stressors["soc_high"]
        Cchg = stressors["Cchg"]
        Cdis = stressors["Cdis"]
        t_days = self.stressors['t_days'][-1] + delta_t_days
        efc = self.stressors['efc'][-1] + delta_efc
        stressors_norm = np.array([delta_t_days, t_days, delta_efc, efc, np.mean(TdegK), np.mean(soc), soc_low, soc_high, np.mean(Ua), dod, Cdis, Cchg])
        for k, v in zip(self.stressors.keys(), stressors_norm):
            self.stressors[k] = np.append(self.stressors[k], v)
            
        self.update_rates(stressors)
        # print(f"kcal: {self.rates['kcal'][-1]}, pcal: {self.rates['pcal'][-1]}, kcyc: {self.rates['kcyc'][-1]}, pcyc: {self.rates['pcyc'][-1]}")
        self.update_states(stressors)
        # print(f"qLoss_t: {self.states['qLoss_t'][-1]}, qLoss_EFC: {self.states['qLoss_EFC'][-1]}")
        self.update_outputs(stressors)
        # print(f"q_t: {self.outputs['q_t'][-1]}, q_EFC: {self.outputs['q_EFC'][-1]}, q: {self.outputs['q'][-1]}")

    def update_rates(self, stressors):
        # Calculate and update battery degradation rates based on stressor values
        # Inputs:
        #   stressors (dict): output from extract_stressors

        # Unpack stressors
        t_secs = stressors["t_secs"]
        delta_t_secs = t_secs[-1] - t_secs[0]
        TdegK = stressors["TdegK"]
        soc = stressors["soc"]
        soc_low = stressors["soc_low"]
        soc_high = stressors["soc_high"]
        Ua = stressors["Ua"]
        dod = stressors["dod"]
        Cdis = stressors["Cdis"]
        Cchg = stressors['Cchg']
        # Normalized temperature and Ua
        TdegKN = TdegK / (273.15+45)
        UaN = Ua / 0.123 # Ua at about 50% SOC for Gr

        # Grab parameters
        p = self._params_life

        # Calculate the degradation coefficients
        kcal = (p['q1_b0']
                * np.exp(p['q1_b1']*(TdegKN**2)*(1/(UaN**(1/3))))
                * np.exp(p['q1_b2']*(TdegKN**2)*(1/(UaN**0.5)))
            )
        pcal = (p['p_cal_b0']
            * np.exp(p['p_cal_b1']*(TdegKN**3))
            * np.exp(p['p_cal_b2']*(TdegKN**3)*(soc**3))
            * np.exp(p['p_cal_b3']*(TdegKN**2)*(1/(UaN**3)))
            )
        kcyc = (p['q3_b0']
                + p['q3_b1']*(Cchg**0.5)*(1/(TdegKN**3))*soc
                + p['q3_b2']*np.exp((Cchg**0.5)*(1/(TdegKN**3))*soc)
                + p['q3_b3']*np.exp(Cdis*(1/(TdegKN**3))*(soc**0.5))
            )
        pcyc = np.abs(p['p_cyc_b0']
                    + p['p_cyc_b1']*np.exp(Cdis*(1/(TdegKN**3))*(soc_high**0.5))
                    + p['p_cyc_b2']*np.exp((dod**0.5)*(Cchg**0.5)*(TdegKN**3)*(soc_low**2))
                    + p['p_cyc_b3']*(TdegKN**2)*soc_low
                    + p['p_cyc_b4']*(dod**0.5)*(Cdis**0.5)*(1/TdegKN)*(Cchg**0.5)
                    + p['p_cyc_b5']*np.exp((dod**0.5)*(Cchg**3)*TdegKN*(Cdis**0.5))
            )
        
        # Calculate time based average of each rate
        kcal = np.trapz(kcal, x=t_secs) / delta_t_secs
        pcal = np.trapz(pcal, x=t_secs) / delta_t_secs
        kcyc = np.trapz(kcyc, x=t_secs) / delta_t_secs
        pcyc = np.trapz(pcyc, x=t_secs) / delta_t_secs

        # kcyc linear model went negative in a few extrapolations on model investigation
        # floor at near the minimum value observed in the experimental data set
        kcyc = np.max(np.array((kcyc, 0.15)))
        if kcyc == 1.5:
            print('kcyc hit floor')

        # Store rates
        rates = np.array([kcal, pcal, kcyc, pcyc])
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
        dq_t = self._degradation_scalar * update_sigmoid_state(states['qLoss_t'][-1], delta_t_days/5e3, r['kcal'], p['q2'], r['pcal'])
        dq_EFC = self._degradation_scalar * update_sigmoid_state(states['qLoss_EFC'][-1], delta_efc/1e5, r['kcyc'], p['q4'], r['pcyc'])

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
    
    @staticmethod
    def extract_stressors(t_secs, soc, T_celsius):
        # extract the usual stressors
        stressors = BatteryDegradationModel.extract_stressors(t_secs, soc, T_celsius)
        # model specific stressors: soc_low, soc_high, Cchg, Cdis
        soc_low = np.min(soc)
        soc_high = np.max(soc)
        t_days = t_secs / (24*60*60)
        delta_t_days = t_days[-1] - t_days[0]
        dt = np.diff(t_secs)
        abs_instantaneous_crate = np.abs(np.diff(soc)/np.diff(t_secs/(60*60))) # get instantaneous C-rates
        abs_instantaneous_crate[abs_instantaneous_crate < 1e-2] = 0 # get rid of extremely small values (storage) before calculating mean
        Crate = np.trapz(abs_instantaneous_crate, t_days[1:]) / delta_t_days
        # Check storage condition, which will give nan Crate:
        if np.isnan(Crate):
            Cdis = 0
            Cchg = 0
        else:
            instantaneous_crate = np.diff(soc)/(dt/3600)
            instantaneous_crate[np.abs(instantaneous_crate) < 1e-2] = 0 # threshold tiny values to zero to prevent weirdness
            mask_cchg = instantaneous_crate > 0
            mask_cdis = instantaneous_crate < 0
            instantaneous_cchg = instantaneous_crate[mask_cchg]
            dt_charging = dt[mask_cchg]
            t_secs_charging = np.cumsum(dt_charging)
            if len(instantaneous_cchg) > 1:
                Cchg = np.trapz(instantaneous_cchg, t_secs_charging) / (t_secs_charging[-1] - t_secs_charging[0])
            elif len(instantaneous_cchg) == 1:
                Cchg = instantaneous_cchg[0]
            else: # half cycle with no charge segment
                Cchg = 0
            instantaneous_cdis = instantaneous_crate[mask_cdis]
            dt_discharging = dt[mask_cdis]
            t_secs_discharging = np.cumsum(dt_discharging)
            if len(instantaneous_cdis) > 1:
                Cdis = np.trapz(np.abs(instantaneous_cdis), t_secs_discharging) / (t_secs_discharging[-1] - t_secs_discharging[0])
            elif len(instantaneous_cdis) == 1:
                Cdis = np.abs(instantaneous_cdis)
                Cdis = Cdis[0]
            else: # half cycle with no discharge segment
                Cdis = 0
        
        stressors['soc_low'] = soc_low
        stressors['soc_high'] = soc_high
        stressors['Cchg'] = Cchg
        stressors['Cdis'] = Cdis
        return stressors

