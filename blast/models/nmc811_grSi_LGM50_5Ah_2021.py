# Paul Gasper, NREL
# This model is fit to LG M50 cell aging data reported by Truong Bui et al, https://ieeexplore.ieee.org/document/9617644 
# Cell chemistry details were reported in a different study by O'Regan et al, https://www.sciencedirect.com/science/article/pii/S0013468622008593 
# Cell tests were reported in 2021, so likely 2019 or 2020 LG MJ1 cells.
# These 21700 LG M50 cells are similar to those found in some EVs, such as Lucid and Rivian, though not sure if those have GrSi composite negative electrodes.

import numpy as np
from blast.models.degradation_model import BatteryDegradationModel

class Nmc811_GrSi_LGM50_5Ah_Battery(BatteryDegradationModel):
    """
    This model is fit to LG M50 cell aging data reported by Truong Bui et al, https://ieeexplore.ieee.org/document/9617644 
    Cell chemistry details were reported in a different study by O'Regan et al, https://www.sciencedirect.com/science/article/pii/S0013468622008593 
    Cell tests were reported in 2021, so likely 2019 or 2020 LG MJ1 cells.
    These 21700 LG M50 cells are similar to those found in some EVs, such as Lucid and Rivian, though not sure if those have GrSi composite negative electrodes.

    .. note::
        EXPERIMENTAL AGING DATA SUMMARY:
            Full factorial design of experiments in SOC and Temperature space for calendar aging, very extensive.
            Calendar aging varied SOC (0, 2, 5, 10, 30, 50, 60 70, 80, 85, 90, 95, 100%) and temperature (0, 25, 45, 60), three cell replicates at every combination of SOC and T.
            Cycle aging varied temperature and discharge C-rate; all DOD is 100%. NO CHARGE RATE VARIATION.
            (NOT 100% CERTAIN) Relative discharge capacity is reported from measurements recorded at 25 Celsius and C/3 rate.
            For cycle aging, the charge-throughput in equivalent full cycles (EFC) was reported not using measured charge-throughput, but rather estimated (cycles*DOD).

        MODEL SENSITIVITY
            The model predicts degradation rate versus time as a function of temperature and average
            state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as 
            a function of discharge C-rate and temperature.

        MODEL LIMITATIONS
            Cycle degradation predictions ARE NOT SENSITIVE TO CHARGING C-RATE OR DEPTH OF DISCHARGE. It is expected that these electrode materials
            would be highly sensitive to depth-of-discharge, see other models with NMC or NCA positive electrodes.
            Cycling aging tests were only performed between 0 and 25 degrees C, assuming cells would be thermally controlled in applications.
    """

    def __init__(self, degradation_scalar: float = 1, label: str = "NMC811-GrSi LG M50"):
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
            'Cdis': np.array([np.nan]),
        }

        # Rates: History of stressor-dependent degradation rates
        self.rates = {
            'k_cal': np.array([np.nan]),
            'k_cyc': np.array([np.nan]),
            'p_cyc': np.array([np.nan]),
        }

        # Expermental range: details on the range of experimental conditions, i.e.,
        # the range we expect the model to be valid in
        self.experimental_range = {
            'cycling_temperature': [0, 25],
            'dod': [1],
            'soc': [0, 1],
            'max_rate_charge': 0.3,
            'max_rate_discharge': 2,
        }

        # Degradation scalar - scales all state changes by a coefficient
        self._degradation_scalar = degradation_scalar
        # Label for plotting
        self._label = label

    # Nominal capacity
    @property
    def cap(self):
        return 5

    # Define life model parameters
    @property
    def _params_life(self):
        return {
            # Capacity fade parameters
            'p_cal': 0.637876353435434,
            'q1_b0': -2.25836422099948,
            'q1_b1': 2.30472337794488,
            'q1_b2': -4.19298863694080,
            'q1_b3': 0.418423177144011,
            'q1_b4': -0.173472466016871,
            'q4': 0.755467311373474,
            'q3_b0': 3.22126947867990e+70,
            'q3_b1': -24.4648676310245,
            'q3_b2': -186.132184606326,
            'q3_b3': 48.5227273096471,
            'q3_b4': 107.854915710757,
            'p_cyc_b0':	3.48011676426504e+129,
            'p_cyc_b1':	-223.262735749431,
            'p_cyc_b2':	-83.1491259202394,
            'p_cyc_b3':	89.7053920930768,
            'p_cyc_b4': 7.30252634692416,
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
        Cdis = stressors["Cdis"]
        TdegKN = TdegK / (273.15 + 35)
        UaN = Ua / 0.123
        
        # Grab parameters
        p = self._params_life

        # Calculate the degradation coefficients
        k_cal = (
            p['q1_b0']
            + p['q1_b1']*np.exp((TdegKN**3)*(1/(UaN**(1/3))))
            + p['q1_b2']*(TdegKN**3)*(1/(UaN**0.5))
            + p['q1_b3']*np.exp((1/(TdegKN**2))*(soc**0.5))
            + p['q1_b4']*np.exp((soc**3))
        )
        
        k_cyc  = (
            p['q3_b0']
            * np.exp(p['q3_b1']*TdegKN*(Cdis**0.5))
            * np.exp(p['q3_b2']*(1/(TdegKN^2))*(Cdis**0.5))
            * np.exp(p['q3_b3']*(1/(TdegKN^3))*Cdis)
            * np.exp(p['q3_b4']*np.log((1/(TdegKN**2))*(Cdis**0.5)))
        )

        p_cyc = (
            p['p_cyc_b0']
            * np.exp(p['p_cyc_b1']*(1/TdegKN))
            * np.exp(p['p_cyc_b2']*(1/(TdegKN**0.5))*Cdis)
            * np.exp(p['p_cyc_b3']*np.log((1/(TdegKN**3))*(Cdis**0.5)))
            * np.exp(p['p_cyc_b4']*(1/TdegKN)*(Cdis**3))
        )

        # Calculate time based average of each rate
        k_cal = np.trapz(k_cal, x=t_secs) / delta_t_secs
        k_cyc = np.trapz(k_cyc, x=t_secs) / delta_t_secs
        p_cyc = np.trapz(p_cyc, x=t_secs) / delta_t_secs

        # Store rates
        rates = np.array([k_cal, k_cyc, p_cyc])
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
        dq_t = self._degradation_scalar * self._update_power_state(states['qLoss_t'][-1], delta_t_days/1e4, r['k_cal'], p['p_cal'])
        dq_EFC = self._degradation_scalar * self._update_sigmoid_state(states['qLoss_EFC'][-1], delta_efc, r['k_cyc'], p['q4'], r['p_cyc'])

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
    def _extract_stressors(t_secs, soc, T_celsius):
        # extract the usual stressors
        stressors = BatteryDegradationModel._extract_stressors(t_secs, soc, T_celsius)
        # model specific stressors: Cdis
        t_days = t_secs / (24*60*60)
        delta_t_days = t_days[-1] - t_days[0]
        dt = np.diff(t_secs)
        abs_instantaneous_crate = np.abs(np.diff(soc)/np.diff(t_secs/(60*60))) # get instantaneous C-rates
        abs_instantaneous_crate[abs_instantaneous_crate < 1e-2] = 0 # get rid of extremely small values (storage) before calculating mean
        Crate = np.trapz(abs_instantaneous_crate, t_days[1:]) / delta_t_days
        # Check storage condition, which will give nan Crate:
        if np.isnan(Crate):
            Cdis = 0
        else:
            instantaneous_crate = np.diff(soc)/(dt/3600)
            instantaneous_crate[np.abs(instantaneous_crate) < 1e-2] = 0 # threshold tiny values to zero to prevent weirdness
            mask_cdis = instantaneous_crate < 0
            instantaneous_cdis = instantaneous_crate[mask_cdis]
            dt_discharging = dt[mask_cdis]
            t_secs_discharging = np.cumsum(dt_discharging)
            if len(instantaneous_cdis) > 1:
                Cdis = np.trapz(instantaneous_cdis, t_secs_discharging) / (t_secs_discharging[-1] - t_secs_discharging[0])
            elif len(instantaneous_cdis) == 1:
                Cdis = instantaneous_cdis[0]
            else: # half cycle with no charge segment
                Cdis = 0
        stressors['Crate'] = Crate
        stressors['Cdis'] = Cdis
        return stressors