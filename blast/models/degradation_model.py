import numpy as np
from ..functions.rainflow import count_cycles

class BatteryDegradationModel:

    def __init__(self):
        # States: Internal states of the battery model
        self.states = {
        }
        # Outputs: Battery properties derived from state values
        self.outputs = {
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
        }
        # Expermental range: details on the range of experimental conditions, i.e.,
        # the range we expect the model to be valid in
        self.experimental_range = {
        }

    # Functions to calculate voltage, half-cell potentials, or other cell specific
    # features that will be needed for calculating stressors. These equations for
    # calculating anode-to-lithium potential for a graphite electrode are a common
    # example, and are taken from M. Safari and C. Delacourt, Journal of the 
    # Electrochemical Society, 158(5), A562 (2011).
    
    @staticmethod
    def get_Ua(soc):
        # Calculate lithiation fraction from soc
        def get_Xa(soc):
            return 8.5*10**-3 + soc*(0.78 - 8.5*10**-3)
        Xa = get_Xa(soc)
        # Calculate Ua from lithiation fraction
        return (0.6379 + 0.5416*np.exp(-305.5309*Xa) + 
            0.044*np.tanh(-1*(Xa-0.1958)/0.1088) - 0.1978*np.tanh((Xa-1.0571)/0.0854) - 
            0.6875*np.tanh((Xa+0.0117)/0.0529) - 0.0175*np.tanh((Xa-0.5692)/0.0875))

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
        
        stressors = BatteryDegradationModel.extract_stressors(t_secs, soc, T_celsius)
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
            
        self.update_rates(stressors)
        self.update_states(stressors)
        self.update_outputs(stressors)

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

        self.update_states(stressors)
        self.update_outputs(stressors)
    
    def update_rates(self, stressors):
        # Calculate and update battery degradation rates based on stressor values
        # Inputs:
        #   stressors (dict): output from extract_stressors
        pass
    
    def update_states(self, stressors):
        # Update the battery states, based both on the degradation state as well as the battery performance
        # at the ambient temperature, T_celsius
        # Inputs:
            #   stressors (dict): output from extract_stressors
        pass

    def update_outputs(self, stressors):
        # Calculate outputs, based on current battery state (and maybe stressors)
        # Inputs:
            #   stressors (dict): output from extract_stressors
        pass

    @staticmethod
    def extract_stressors(t_secs, soc, T_celsius):
        # Extract stressors
        t_days = t_secs / (24*60*60)
        delta_t_days = t_days[-1] - t_days[0]
        delta_efc = np.sum(np.abs(np.ediff1d(soc, to_begin=0)))/2  # sum the total changes to SOC / 2
        dod = np.max(soc) - np.min(soc)
        abs_instantaneous_crate = np.abs(np.diff(soc)/np.diff(t_secs/(60*60))) # get instantaneous C-rates
        abs_instantaneous_crate[abs_instantaneous_crate < 1e-2] = 0 # get rid of extremely small values (storage) before calculating mean
        Crate = np.trapz(abs_instantaneous_crate, t_days[1:]) / delta_t_days
        # Check storage condition, which will give nan Crate:
        if np.isnan(Crate):
            Crate = 0
        T_kelvin = T_celsius + 273.15
        # Estimate Ua (anode to reference potential) from SOC.
        # Uses the equation from Safari and Delacourt, https://doi.org/10.1149/1.3567007.
        # Anode stoichiometry is assumed to be the same for any chemistry/cell, and is calculated using the equation from Schimpe et al https://doi.org/10.1149/2.1181714jes
        # While this will not be precise, it still helps get a guess as to where the plateaus of the anode-reference potential are.
        Ua = BatteryDegradationModel.get_Ua(soc)

        cycles = count_cycles(soc)
        cycles = sum(i for _, i in cycles)

        stressors = {
            't_secs': t_secs,
            'delta_t_days': delta_t_days,
            'delta_efc': delta_efc,
            'TdegK': T_kelvin,
            'soc': soc,
            'dod': dod,
            'Ua': Ua,
            'Crate': Crate,
            'cycles': cycles
        }

        return stressors