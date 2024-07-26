# Paul Gasper, NREL
import numpy as np
from blast.models.degradation_model import BatteryDegradationModel

# EXPERIMENTAL AGING DATA SUMMARY:
# Aging test matrix varied primarly temperature, with small DOD variation.
# Calendar and cycle aging were performed between 0 and 55 Celsius. C-rates always at 1C,
# except for charging at 0 Celsius, which was conducted at C/3. Depth-of-discharge was 80%
# for nearly all tests (3.4 V - 4.1 V), with one 100% DOD test (3 V - 4.2 V).
# Reported relative capacity was measured at C/5 rate at the aging temperatures. Reported
# relative DC resistance was measured by HPPC using a 10s, 1C DC pulse, averaged between 
# charge and discharge, calculated using a simple ohmic fit of the voltage response.

# MODEL SENSITIVITY
# The model predicts degradation rate versus time as a function of temperature and average
# state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as 
# a function of temperature and depth-of-discharge. Sensitivity to cycling degradation rate
# at low temperature is inferred from physical insight due to limited data.

# MODEL LIMITATIONS
# There is NO C-RATE DEPENDENCE for degradation in this model. THIS IS NOT PHYSICALLY REALISTIC
# AND IS BASED ON LIMITED DATA.


class Nmc111_Gr_Kokam75Ah_Battery(BatteryDegradationModel):
    # Model predicting the degradation of a Kokam 75 Ah NMC-Gr pouch cell.
    # https://ieeexplore.ieee.org/iel7/7951530/7962914/07963578.pdf
    # It is uncertain if the exact NMC composition is 1-1-1, but it this is definitely not a high nickel (>80%) cell.
    # Degradation rate is a function of the aging stressors, i.e., ambient temperature and use.
    # The state of the battery is updated throughout the lifetime of the cell.
    # Performance metrics are capacity and DC resistance. These metrics change as a function of the
    # cell's current degradation state, as well as the ambient temperature. The model predicts time and 
    # cycling dependent degradation, using Loss of Lithium Inventory (LLI) and Loss of Active
    # Material (LAM) degradation modes that interact competitively (cell performance is limited by
    # one or the other.)
    # Parameters to modify to change fade rates:
    #   Calendar capacity loss rate: q1_0
    #   Cycling capacity loss rate (LLI): q3_0
    #   Cycling capacity loss rate (LAM): q5_0, will also effect resistance growth onset due to LAM.
    #   Calendar resistance growth rate (LLI), relative to capacity loss rate: r1
    #   Cycling resistance growth rate (LLI), relative to capacity loss rate: r3

    def __init__(self, degradation_scalar=1, label="NMC111-Gr Kokam"):
        # States: Internal states of the battery model
        self.states = {
            'qLoss_LLI_t': np.array([0]),   # relative Li inventory change, time dependent (SEI)
            'qLoss_LLI_EFC': np.array([0]), # relative Li inventory change, charge-throughput dependent (SEI)
            'qLoss_LAM': np.array([1e-8]),  # relative active material change, charge-throughput dependent (electrode damage)
            'rGain_LLI_t': np.array([0]),   # relative SEI growth, time dependent (SEI)
            'rGain_LLI_EFC': np.array([0]), # relative SEI growth, charge-throughput dependent (SEI)
        }

        # Outputs: Battery properties derived from state values
        self.outputs = {
            'q': np.array([1]),         # relative capacity
            'q_LLI': np.array([1]),     # relative lithium inventory
            'q_LLI_t': np.array([1]),   # relative lithium inventory, time dependent loss
            'q_LLI_EFC': np.array([1]), # relative lithium inventory, charge-throughput dependent loss
            'q_LAM': np.array([1.01]),  # relative active material, charge-throughput dependent loss
            'r': np.array([1]),         # relative resistance
            'r_LLI': np.array([1]),     # relative SEI resistance
            'r_LLI_t': np.array([1]),   # relative SEI resistance, time dependent growth
            'r_LLI_EFC': np.array([1]), # relative SEI resistance, charge-throughput dependent growth
            'r_LAM': np.array([1]),     # relative electrode resistance, q_LAM dependent growth
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
        }

        # Rates: History of stressor-dependent degradation rates
        self.rates = {
            'q1': np.array([np.nan]),
            'q3': np.array([np.nan]),
            'q5': np.array([np.nan]),
        }

        # Expermental range: details on the range of experimental conditions, i.e.,
        # the range we expect the model to be valid in
        self.experimental_range = {
            'cycling_temperature': [0, 45],
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
    def _cap(self):
        return 75

    # Define life model parameters
    @property
    def _params_life(self):
        return {
            'q1_0'  : 2.66e7,   # CHANGE to modify calendar degradation rate (larger = faster degradation)
            'q1_1'  : -17.8,
            'q1_2'  : -5.21,
            'q2'    : 0.357,
            'q3_0'  : 3.80e3,   # CHANGE to modify cycling degradation rate (LLI) (larger = faster degradation)
            'q3_1'  : -18.4,
            'q3_2'  : 1.04,
            'q4'    : 0.778,
            'q5_0'  : 1e4,      # CHANGE to modify cycling degradation rate (LAM) (accelerating fade onset) (larger = faster degradation)
            'q5_1'  : 153,
            'p_LAM' : 10,
            'r1'    : 0.0570,   # CHANGE to modify change of resistance relative to change of capacity (calendar degradation)
            'r2'    : 1.25,
            'r3'    : 4.87,     # CHANGE to modify change of resistance relative to change of capacity (cycling degradation)
            'r4'    : 0.712,
            'r5'    : -0.08,
            'r6'    : 1.09,
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
        TdegC = TdegK - 273.15
        TdegKN = TdegK / (273.15 + 35) # normalized temperature
        UaN = Ua / 0.123 # normalized anode-to-reference potential
        
        # Grab parameters
        p = self._params_life

        # Calculate degradation rates
        q1 = p['q1_0'] * np.exp(p['q1_1'] * (1 / TdegKN)) * np.exp(p['q1_2'] * (UaN / TdegKN))
        q3 = p['q3_0'] * np.exp(p['q3_1'] * (1/TdegKN)) * np.exp(p['q3_2'] * np.exp(dod**2))
        q5 = p['q5_0'] + p['q5_1'] * (TdegC - 55) * dod * (np.mean(soc)/0.6)
        # Calculate time based average of each rate
        q1 = np.trapz(q1, x=t_secs) / delta_t_secs
        q3 = np.trapz(q3, x=t_secs) / delta_t_secs
        q5 = np.trapz(q5, x=t_secs) / delta_t_secs

        # Store rates
        rates = np.array([q1, q3, q5])
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
        dq_LLI_t = self._degradation_scalar * self.update_power_state(states['qLoss_LLI_t'][-1], delta_t_days, r['q1'], p['q2'])
        dq_LLI_EFC = self._degradation_scalar * self.update_power_state(states['qLoss_LLI_EFC'][-1], delta_efc, r['q3'], p['q4'])
        dq_LAM = self._degradation_scalar * self.update_sigmoid_state(states['qLoss_LAM'][-1], delta_efc, 1, 1/r['q5'], p['p_LAM'])
        
        # Resistance
        dr_LLI_t = self._degradation_scalar * self.update_power_state(states['rGain_LLI_t'][-1], delta_t_days, p['r1']*r['q1'], p['r2'])
        dr_LLI_EFC = self._degradation_scalar * self.update_power_state(states['rGain_LLI_EFC'][-1], delta_efc, p['r3']*r['q3'], p['r4'])

        # Accumulate and store states
        dx = np.array([dq_LLI_t, dq_LLI_EFC, dq_LAM, dr_LLI_t, dr_LLI_EFC])
        for k, v in zip(states.keys(), dx):
            x = self.states[k][-1] + v
            self.states[k] = np.append(self.states[k], x)

    def update_outputs(self, stressors):
        # Calculate outputs, based on current battery state
        states = self.states
        p = self._params_life

        # Capacity
        q_LLI = 1 - states['qLoss_LLI_t'][-1] - states['qLoss_LLI_EFC'][-1]
        q_LLI_t = 1 - states['qLoss_LLI_t'][-1]
        q_LLI_EFC = 1 - states['qLoss_LLI_EFC'][-1]
        q_LAM = 1.01 - states['qLoss_LAM'][-1]
        q = np.min(np.array([q_LLI, q_LAM]))
        
        # Resistance
        r_LLI = 1 + states['rGain_LLI_t'][-1] + states['rGain_LLI_EFC'][-1]
        r_LLI_t = 1 + states['rGain_LLI_t'][-1]
        r_LLI_EFC = 1 + states['rGain_LLI_EFC'][-1]
        r_LAM = p['r5'] + p['r6'] * (1 / q_LAM)
        r = np.max(np.array([r_LLI, r_LAM]))

        # Assemble output
        out = np.array([q, q_LLI, q_LLI_t, q_LLI_EFC, q_LAM, r, r_LLI, r_LLI_t, r_LLI_EFC, r_LAM])
        # Store results
        for k, v in zip(list(self.outputs.keys()), out):
            self.outputs[k] = np.append(self.outputs[k], v)