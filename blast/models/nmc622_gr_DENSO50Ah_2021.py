# Paul Gasper and Kandler Smith, NREL
# This model is replicated as reported by Smith et al, J. Echem. Soc. (2021) 168 100530
# http://dx.doi.org/10.1149/1945-7111/ac2ebd

import numpy as np
from blast.models.degradation_model import BatteryDegradationModel


class Nmc622_Gr_DENSO50Ah_Battery(BatteryDegradationModel):
    """
    Model predicting the degradation of NMC622-Gr EV cells from DENSO, reported by Smith et al, J. 
    Echem. Soc. (2021) 168 100530 http://dx.doi.org/10.1149/1945-7111/ac2ebd.

    .. note::
        EXPERIMENTAL AGING DATA SUMMARY:
            Calendar aging across many conditions between 10C and 60C, 10% and 100% SOC.
            Cycle aging between 10C and 60C, varying DOD at 25C and 45C, charging C-rate at 25C to 55C,
            discharge rate constant at 1C for all cycle aging tests.
            Relative discharge capacity is reported from measurements recorded at 25C and C/3 rate.

        MODEL SENSITIVITY
            The model predicts degradation rate versus time as a function of temperature and average
            state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as 
            a function of average voltage and depth-of-discharge.

        MODEL LIMITATIONS
            Cycle aging was always conducted at 50% average SOC. No high C-rate testing (1C charge) was
            conducted at less than 25C, and the experimental data showed no 'knee over' events, so no
            accelerated fade due to lithium plating onset is predicted by this model.
    """
    
    def __init__(self, degradation_scalar: float = 1, label: str = "NMC622-Gr DENSO"):
        # States: Internal states of the battery model
        self.states = {
            'qLoss_t': np.array([0]),
            'qLoss_EFC': np.array([0]),
            'qLoss_BreakIn_t': np.array([0]),     
            'qLoss_BreakIn_EFC': np.array([0]),     
        }

        # Outputs: Battery properties derived from state values
        self.outputs = {
            'q': np.array([1]),
            'q_t': np.array([1]),
            'q_EFC': np.array([1]),
            'q_BreakIn': np.array([1]),
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
            'Cchg': np.array([np.nan]),
        }

        # Rates: History of stressor-dependent degradation rates
        self.rates = {
            'b1_t': np.array([np.nan]),
            'b1_N': np.array([np.nan]),
            'b3_t': np.array([np.nan]),
            'b3_N': np.array([np.nan]),
        }

        # Expermental range: details on the range of experimental conditions, i.e.,
        # the range we expect the model to be valid in
        self.experimental_range = {
            'cycling_temperature': [10, 60],
            'dod': [0.1, 1],
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
    def cap(self):
        return 50

    # Define life model parameters
    @property
    def _params_life(self):
        return {
            # Capacity fade parameters
            'b1t_Ea': 37000,
            'b1t_0': -0.000197,
            'b1t_1': 0.0101,
            'b1t_2': -0.0157,
            'b1t_3': 0.00835,
            'b1t_3': -4.06E-6,
            'b1t_4': 3.32E-5,
            'b1N_Ea_1': -58000,
            'b1N_0': 3.24,
            'b1N_1': 2.1,
            'b1N_2': -0.099,
            'b1N_3': 1.44,
            'b1N_Ea_2': -13000,
            'b1N_4': 0.985,
            'b3t_0': -0.0303,
            'b3t_1': 0.269,
            'b3t_2': -1.360,
            'b3t_3': 0.208,
            'b3t_4': -0.272,
            'b3N_0': 0.0791,
            'b3N_Ea': -8800,
            'b3N_1': 1.143,
            'b3N_2': -0.0386,
            'b3N_3': 0.178,          
        }
    
    def arrhenius(TdegK, Ea):
        return np.exp(-(Ea/8.3144)*((1/TdegK)-(1/328.15)))

    # Battery model
    def update_rates(self, stressors):
        # Calculate and update battery degradation rates based on stressor values
        # Inputs:
        #   stressors (dict): output from extract_stressors

        # Unpack stressors
        t_secs = stressors["t_secs"]
        delta_t_secs = t_secs[-1] - t_secs[0]
        TdegK = stressors["TdegK"]
        soc = np.mean(stressors["soc"])
        dod = stressors["dod"]
        Cchg = stressors["Cchg"]
        stress = (dod*Cchg)**0.5

        # Grab parameters
        p = self._params_life

        # Calculate the degradation coefficients
        b1t = self.arrhenius(TdegK, p['b1t_Ea'])*(p['b1t_0'] + p['b1t_1']*soc + p['b1t_2']*soc**2 + p['b1t_3']*soc**3 + p['b1t_4']*soc*TdegK + p['b1t_5']*np.max([0, TdegK-328.15]))
        b1N = np.max([0, p['b1N_0']*self.arrhenius(TdegK, p['b1N_Ea_1'])*stress*(1+p['b1N_1']*soc) + p['b1N_2']]) + p['b1N_3']*self.arrhenius(TdegK, p['b1N_Ea_2'])*(dod**6)
        b3t = p['b3t_0'] + p['b3t_1']*(1+p['b3t_2']*soc) + p['b3t_3']*np.max([0, soc-0.3]) + p['b3t_4']*np.max([0, 0.9-soc])
        b3N = np.max([0, p['b3N_0']*self.arrhenius(TdegK, p['b3N_Ea'])*stress*(1+p['b3N_1']*soc)+p['b3N_2']]) + p['b3N_3']*np.max([0, dod-0.85])

        # Calculate time based average of each rate
        b1t = np.trapz(b1t, x=t_secs) / delta_t_secs
        b1N = np.trapz(b1N, x=t_secs) / delta_t_secs
        b3t = np.trapz(b3t, x=t_secs) / delta_t_secs
        b3N = np.trapz(b3N, x=t_secs) / delta_t_secs

        # Store rates
        rates = np.array([b1t, b1N, b3t, b3N])
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
        dq_t = self._degradation_scalar * self._update_power_state(states['qLoss_t'][-1], delta_t_days, r['b1t'], 0.5)
        dq_EFC = self._degradation_scalar * self._update_power_state(states['qLoss_EFC'][-1], delta_efc, r['b1t']*r['b1N']*p['b1N_4'], 0.5)
        dq_BreakIn_t = self._degradation_scalar * self._update_exponential_relax_state(states['qLoss_BreakIn_t'][-1], delta_t_days, r['b3t'], 10)
        dq_BreakIn_EFC = self._degradation_scalar * self._update_exponential_relax_state(states['qLoss_BreakIn_EFC'][-1], delta_t_days, r['b3N'], 10)

        # Accumulate and store states
        dx = np.array([dq_t, dq_EFC, dq_BreakIn_t, dq_BreakIn_EFC])
        for k, v in zip(states.keys(), dx):
            x = self.states[k][-1] + v
            self.states[k] = np.append(self.states[k], x)
    
    def update_outputs(self, stressors):
        # Calculate outputs, based on current battery state
        states = self.states

        # Capacity
        q_t = 1 - states['qLoss_t'][-1]
        q_EFC = 1 - states['qLoss_EFC'][-1]
        q_BreakIn = 1 - (states['qLoss_BreakIn_t'][-1] + states['qLoss_BreakIn_EFC'][-1])
        q = 1 - states['qLoss_t'][-1] - states['qLoss_EFC'][-1] - (states['qLoss_BreakIn_t'][-1] + states['qLoss_BreakIn_EFC'][-1])

        # Assemble output
        out = np.array([q, q_t, q_EFC, q_BreakIn])
        # Store results
        for k, v in zip(list(self.outputs.keys()), out):
            self.outputs[k] = np.append(self.outputs[k], v)

    @staticmethod
    def _extract_stressors(t_secs, soc, T_celsius):
        # extract the usual stressors
        stressors = BatteryDegradationModel._extract_stressors(t_secs, soc, T_celsius)
        # model specific stressors: Cchg
        t_days = t_secs / (24*60*60)
        delta_t_days = t_days[-1] - t_days[0]
        dt = np.diff(t_secs)
        abs_instantaneous_crate = np.abs(np.diff(soc)/np.diff(t_secs/(60*60))) # get instantaneous C-rates
        abs_instantaneous_crate[abs_instantaneous_crate < 1e-2] = 0 # get rid of extremely small values (storage) before calculating mean
        Crate = np.trapz(abs_instantaneous_crate, t_days[1:]) / delta_t_days
        # Check storage condition, which will give nan Crate:
        if np.isnan(Crate):
            Cchg = 0
        else:
            instantaneous_crate = np.diff(soc)/(dt/3600)
            instantaneous_crate[np.abs(instantaneous_crate) < 1e-2] = 0 # threshold tiny values to zero to prevent weirdness
            mask_cchg = instantaneous_crate > 0
            instantaneous_cchg = instantaneous_crate[mask_cchg]
            dt_charging = dt[mask_cchg]
            t_secs_charging = np.cumsum(dt_charging)
            if len(instantaneous_cchg) > 1:
                Cchg = np.trapz(instantaneous_cchg, t_secs_charging) / (t_secs_charging[-1] - t_secs_charging[0])
            elif len(instantaneous_cchg) == 1:
                Cchg = instantaneous_cchg[0]
            else: # half cycle with no charge segment
                Cchg = 0
        stressors['Cchg'] = Cchg
        return stressors
