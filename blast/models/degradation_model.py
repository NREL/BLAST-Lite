import numpy as np
from blast.utils.rainflow import count_cycles, extract_cycles
from blast.utils.functions import rescale_soc
import pandas as pd
from typing import Union
import matplotlib.pyplot as plt


class BatteryDegradationModel:
    def __init__(self):
        self.states = {}
        """dict: Internal states of the battery model"""
        self.outputs = {}
        """dict: Battery properties derived from state values"""
        self.stressors = {
            "delta_t_days": np.array([np.nan]),
            "t_days": np.array([0]),
            "delta_efc": np.array([np.nan]),
            "efc": np.array([0]),
            "TdegK": np.array([np.nan]),
            "soc": np.array([np.nan]),
            "Ua": np.array([np.nan]),
            "dod": np.array([np.nan]),
            "Crate": np.array([np.nan]),
        }
        """dict: History of stressors on the battery"""
        self.rates = {}
        """dict: History of stressor-dependent degradation rates"""
        self.experimental_range = {}
        """dict: details on the range of experimental conditions, i.e.,
        the range we expect the model to be valid in"""

        self._label = ""
        """str: Model name for plotting"""

    # Functions for updating time-varying states
    @staticmethod
    def _update_power_state(y0: float, dx: float, k: float, p: float):
        """
        Update time-varying power state
        Trajectory equation: y = k*x^p
        State equation: dy = k*p*(k/y)^((1-p)/p)

        Args:
            y0 (float): current state value
            dx (float): change of the independent variable over this timestep
            k (float): coefficient of x^p
            p (float): power of x

        Return:
            dy (float): the change of the state over this timestep
        """
        if y0 == 0:
            if dx == 0:
                dydx = 0
            else:
                y0 = k * (dx**p)
                dydx = y0 / dx
        else:
            if dx == 0:
                dydx = 0
            else:
                dydx = k * p * ((y0 / k) ** ((p - 1) / p))
        return dydx * dx

    @staticmethod
    def _update_power_B_state(y0: float, dx: float, k: float, p: float):
        """
        Update time-varying power B state
        Trajectory equation: y = (k*x)^p
        State equation: 
            z = (y^(1/p))/k
            dy = (p*(k*z)^p)/z

        Args:
            y0 (float): current state value
            dx (float): change of the independent variable over this timestep
            k (float): coefficient of x
            p (float): power of k*x

        Return:
            dy (float): the change of the state over this timestep
        """
        if y0 == 0:
            if dx == 0:
                dydx = 0
            else:
                y0 = (k * dx) ** p
                dydx = y0 / dx
        else:
            if dx == 0:
                dydx = 0
            else:
                z = (y0 ** (1 / p)) / k
                dydx = (p * (k * z) ** p) / z
        return dydx * dx

    @staticmethod
    def _update_sigmoid_state(y0: float, dx: float, y_inf: float, k: float, p: float):
        """
        Update time-varying sigmoid state
        Trajectory equation:
            y = 2*y_inf*(1/2-(1/(1+exp((k*x)^p)))
        State equation:
            x_inv = (1/k)*log(-(2*y_inf./(y0-y_inf))-1)^(1/p)
            z = (k*x_inv)^p
            dy = (2*y_inf*p*exp(z)*z)/(x_inv*(exp(z)+1)^2)

        Args:
            y0 (float): current state value
            dx (float): change of the independent variable over this timestep
            y_inf (float): maximum extent of the sigmoid
            k (float): 1/x of the half maximum of the sigmoid
            p (float): curvature of the sigmoid at its half maximum

        Return:
            dy (float): the change of the state over this timestep
        """
        if y0 == 0:
            if dx == 0:
                dydx = 0
            else:
                dy = 2 * y_inf * (1 / 2 - 1 / (1 + np.exp((k * dx) ** p)))
                dydx = dy / dx
        else:
            if dx == 0:
                dydx = 0
            else:
                x_inv = (1 / k) * ((np.log(-(2 * y_inf / (y0 - y_inf)) - 1)) ** (1 / p))
                z = (k * x_inv) ** p
                dydx = (2 * y_inf * p * np.exp(z) * z) / (x_inv * (np.exp(z) + 1) ** 2)
        return dydx * dx
    
    @staticmethod
    def _update_exponential_relax_state(y0: float, dx: float, y_inf: float, tau: float):
        """
        Update time-varying exponential relaxation state
        Trajectory equation: y = yInf*(1-exp(-tau*x))
        State equation: dy = (yInf^2*tau)/(yInf-y)

        Args:
            y0 (float): current state value
            dx (float): change of the independent variable over this timestep
            y_inf (float): maximum extent of the exponential relaxation curve
            tau (float): time constant of the exponential relaxation curve

        Return:
            dy (float): the change of the state over this timestep
        """
        if dx == 0:
            dydx = 0
        else:
            dydx = (tau*y_inf**2) / (y_inf-y0)
            # Avoid negative values (when y > y_inf)
            if dydx < 0:
                dydx = 0
        return dydx * dx

    # Functions to calculate voltage, half-cell potentials, or other cell specific
    # features that will be needed for calculating stressors. These equations for
    # calculating anode-to-lithium potential for a graphite electrode are a common
    # example, and are taken from M. Safari and C. Delacourt, Journal of the
    # Electrochemical Society, 158(5), A562 (2011).
    @staticmethod
    def _get_Ua(soc: np.ndarray):
        """
        Calculate Ua from SOC via lithiation fraction.

        Args:
            soc (np.ndarray)    SOC profile

        Return:
            Ua
        """
        # Calculate lithiation fraction from soc
        def get_Xa(soc):
            return 8.5 * 10**-3 + soc * (0.78 - 8.5 * 10**-3)

        Xa = get_Xa(soc)
        # Calculate Ua from lithiation fraction
        return (
            0.6379
            + 0.5416 * np.exp(-305.5309 * Xa)
            + 0.044 * np.tanh(-1 * (Xa - 0.1958) / 0.1088)
            - 0.1978 * np.tanh((Xa - 1.0571) / 0.0854)
            - 0.6875 * np.tanh((Xa + 0.0117) / 0.0529)
            - 0.0175 * np.tanh((Xa - 0.5692) / 0.0875)
        )

    def simulate_battery_life(
        self,
        input_timeseries: Union[dict, pd.DataFrame],
        threshold_time: float = None,
        threshold_capacity: float = None,
        threshold_efc: float = None,
        is_constant_input: bool = False,
        is_conserve_energy_throughput: bool = True,
        breakpoints_max_time_diff_s: float = 86400,
        breakpoints_max_EFC_diff: float = 1,
    ):
        """
        Run battery life simulation over the input, optionally repeating for any number of years or until
        the battery capacity falls below a certain threshold.

        Updates attributes self.rates, self.stressors, self.outputs, and self.states inplace.

        Args:
            input_timeseries (dict, pd.DataFrame)   Input time, SOC, and temperature series
            threshold_time (float)                  Threshold number of years for stopping simulation
            threshold_capacity (float)              Threshold capacity for stopping simulation
            threshold_efc (float)                   Threshold number of equivalent full cycles for stopping simulation
            is_constant_input (bool)                Whether the input is constant
            is_conserve_energy_throughput (bool)    Whether to conserve energy throughput as simulation proceeds.
                                                    This effects how the SOC profile and stressors evolve over lifetime.
                                                    If true, energy throughput is conserved, that is, the depth-of-discharge
                                                    will increase as the battery degrades. If false, the EFCs (i.e., energy throughput)
                                                    of any given SOC swing will decrease. In real life, the time it takes for a cycle
                                                    to complete would also decrease, but this is beyond scope of this package.
            breakpoints_max_time_diff_s (float)     Max time difference for finding breakpoints
            breakpoints_max_EFC_diff (float)        Max EFC difference for finding breakpoints
        """

        # Check if required keys are in the input_timeseries dictionary
        required_keys = ["Temperature_C", "SOC", "Time_s"]
        for key in required_keys:
            if key not in input_timeseries:
                raise ValueError(
                    f"Required key '{key}' is missing in the input_timeseries dictionary."
                )
        # Check if all values are numpy arrays and of the same length
        lengths = [len(input_timeseries[key]) for key in required_keys]
        if len(set(lengths)) != 1:
            raise ValueError(
                "All numpy arrays in input_timeseries must have the same length."
            )

        # Unpack the inputs
        if isinstance(input_timeseries, dict):
            t_secs = input_timeseries["Time_s"]
            soc = input_timeseries["SOC"]
            temperature = input_timeseries["Temperature_C"]
        elif isinstance(input_timeseries, pd.DataFrame):
            t_secs = input_timeseries["Time_s"].values
            soc = input_timeseries["SOC"].values
            temperature = input_timeseries["Temperature_C"].values
        else:
            raise TypeError("'input_timeseries' was not a dict or a Dataframe.")

        # Check stopping criteria (thresholds)
        # If all are none, infer that we just run the simulation over the input timeseries once, so set threshold_time to t_secs[-1]
        if threshold_time is None and threshold_capacity is None and threshold_efc is None:
            threshold_time = (t_secs[-1] - 1) / (365 *24 * 60 * 60) # 1 s threshold to make sure we cross this limit

        if is_constant_input:
            # Run life sim repeating until simulation is longer than threshold_time
            # Only calculate stressors / degradation rates on the first timestep to dramatically accelerate the simulation
            self.update_battery_state(
                t_secs,
                soc,
                temperature,
            )
            is_simulation_complete = False
            while not is_simulation_complete:
                if threshold_time is not None and self.stressors["t_days"][-1] / 365 > threshold_time:
                    is_simulation_complete = True
                    break
                if threshold_capacity is not None and self.outputs['q'][-1] < threshold_capacity:
                    is_simulation_complete = True
                    break
                if threshold_efc is not None and self.stressors["efc"][-1] > threshold_efc:
                    is_simulation_complete = True
                    break
                self.update_battery_state_repeating(is_conserve_energy_throughput)
    
            return self
        else:
            # Chunk the input timeseries into cycles (SOC turning points) using the rainflow algorithm
            cycles_generator = (
                [span, mean, count, idx_start, idx_end]
                for span, mean, count, idx_start, idx_end in extract_cycles(soc)
            )
            idx_turning_point = []
            for cycle in cycles_generator:
                idx_turning_point.append(cycle[-1])
            idx_turning_point = np.array(idx_turning_point)
            # sort turning points and remove any duplicates (sometimes full cycles and half cycles overlap)
            idx_turning_point = np.unique(idx_turning_point)

            # Low DOD cycles can be extremely short, or the cycle time during long periods of storage excessively long.
            # Find breakpoints for simulating degradation, either simulating degradation when the number of EFCs between
            # turning points is greater than 1, or every day if EFCs accumlate slower than that.
            EFCs = np.cumsum(np.abs(np.ediff1d(soc, to_begin=0))) / 2
            simulation_breakpoints = self._find_breakpoints(
                t_secs,
                EFCs,
                idx_turning_point,
                max_time_diff_s=breakpoints_max_time_diff_s,
                max_EFC_diff=breakpoints_max_EFC_diff,
            )
            # Add the final point manually if it's not there by coinkydink, or if there are no cycles in the inputs
            if simulation_breakpoints == []:
                simulation_breakpoints.append(len(t_secs) - 1)
            elif simulation_breakpoints[-1] != len(t_secs) - 1:
                simulation_breakpoints.append(len(t_secs) - 1)

            # Simulate battery life between breakpoints until we reach the end of the simulation
            is_simulation_complete = False
            soc_original = soc.copy()
            while not is_simulation_complete:
                prior_breakpoint = 0
                soc = soc_original.copy()
                for breakpoint in simulation_breakpoints:
                    if is_conserve_energy_throughput:
                        # Rescale SOC according to SOH to reflect changing battery performance over its life
                        # If this isn't true, the '_extract_stressors' function will automatically decrease dEFC
                        # by the health of the battery to reflect lower charge throughput from a cycle as the battery ages.
                        soc_ = soc[prior_breakpoint:breakpoint + 1]
                        soc_ = rescale_soc(soc_, 1/self.outputs["q"][-1])
                        soc[prior_breakpoint:breakpoint + 1] = soc_
                        if breakpoint < len(soc) - 1:
                            soc[breakpoint + 1:] += (soc_[-1] - soc[breakpoint])
                        # plt.plot(t_secs, soc)
                    # Update battery state between breakpoints
                    self.update_battery_state(
                        t_secs[prior_breakpoint:breakpoint + 1],
                        soc[prior_breakpoint:breakpoint + 1],
                        temperature[prior_breakpoint:breakpoint + 1],
                    )
                    prior_breakpoint = breakpoint + 1
                    if threshold_time is not None and self.stressors["t_days"][-1] / 365 > threshold_time:
                        is_simulation_complete = True
                        break
                    if threshold_capacity is not None and self.outputs['q'][-1] < threshold_capacity:
                        is_simulation_complete = True
                        break
                    if threshold_efc is not None and self.stressors["efc"][-1] > threshold_efc:
                        is_simulation_complete = True
                        break
            
            return self

    @staticmethod
    def _find_breakpoints(
        time_seconds: np.ndarray,
        EFCs: np.ndarray,
        idx_turning_point: np.ndarray,
        max_time_diff_s=86400,
        max_EFC_diff=1,
    ):
        """
        Find breakpoints for simulating battery degradation, where degradation is calculated once
        a certain number of equivalent full cycles has passed, otherwise after a certain time has passed.
        Defaults are to calculate degradation every 1 EFC or every day.

        Args:
            time_seconds (np.ndarray)       Time in seconds
            EFCs (np.ndarray)               Equivalent full cycles
            idx_turning_point (np.ndarray)  List of indices in the SOC profile that represent
                                            the end of cycles
            max_time_diff_s (int)           Maximum time difference in seconds
            max_EFC_diff (int)              Maximum difference in equivalent full cycles

        Return:
            breakpoints (list)
        """

        breakpoints = []
        last_breakpoint = 0
        for i in idx_turning_point:
            # Check if change in time_seconds greater than max_time_diff_s
            if time_seconds[i] - time_seconds[last_breakpoint] > max_time_diff_s:
                # Add a breakpoint at each day
                for i_time in range(last_breakpoint, i):
                    if (
                        time_seconds[i_time] - time_seconds[last_breakpoint]
                        > max_time_diff_s
                    ):
                        breakpoints.append(i_time)
                        last_breakpoint = i_time
                continue

            # Check change in EFCs greater than max_EFC_diff
            if EFCs[i] - EFCs[last_breakpoint] > max_EFC_diff:
                breakpoints.append(i)
                last_breakpoint = i

        return breakpoints

    def update_battery_state(
        self, t_secs: np.ndarray, soc: np.ndarray, T_celsius: np.ndarray
    ):
        """
        Update the battery states, based both on the degradation state as well as the battery performance
        at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.

        Args:
            t_secs (np.ndarray)     Vector of the time in seconds since beginning of life
                                    for the soc_timeseries data points
            soc (np.ndarray):       Vector of the state-of-charge of the battery at each t_sec
            T_celsius (ndarray)     Temperature of the battery during this time period, in Celsius units.

        """

        # Check some input types:
        if not isinstance(t_secs, np.ndarray):
            raise TypeError('Input "t_secs" must be a numpy.ndarray')
        if not isinstance(soc, np.ndarray):
            raise TypeError('Input "soc" must be a numpy.ndarray')
        if not isinstance(T_celsius, np.ndarray):
            raise TypeError('Input "T_celsius" must be a numpy.ndarray')
        if not (len(t_secs) == len(soc) and len(t_secs) == len(T_celsius)):
            raise ValueError("All input timeseries must be the same length")

        stressors = self._extract_stressors(t_secs, soc, T_celsius)
        for key in stressors:
            if key == "delta_t_days":
                # Accumulate value
                self.stressors[key] = np.append(self.stressors[key], stressors[key])
                self.stressors["t_days"] = np.append(self.stressors["t_days"], self.stressors["t_days"][-1] + stressors[key])
            elif key == "delta_efc":
                # Accumulate value
                self.stressors[key] = np.append(self.stressors[key], stressors[key])
                self.stressors["efc"] = np.append(self.stressors["efc"], self.stressors["efc"][-1] + stressors[key])
            elif key == "TdegK" or key == "soc" or key == "Ua":
                # Average value
                if key in self.stressors.keys():
                    self.stressors[key] = np.append(self.stressors[key], np.mean(stressors[key]))
            else:
                if key in self.stressors.keys():
                    self.stressors[key] = np.append(self.stressors[key], stressors[key])

        self.update_rates(stressors)
        self.update_states(stressors)
        self.update_outputs(stressors)

    def update_battery_state_repeating(self, is_conserve_energy_throughput=True):
        """
        Update the battery states, based both on the degradation state as well as the battery performance
        at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
        degradation rates are unchanging for every timestep, and don't need to be calculated again.

        Updates self.states and self.outputs inplace.
        """

        # Just accumulate the cumulative stressors (total time, charge throughput in units of EFCs)
        self.stressors["t_days"] = np.append(
            self.stressors["t_days"],
            self.stressors["t_days"][-1] + self.stressors["delta_t_days"][-1],
        )
        if is_conserve_energy_throughput:
            self.stressors["efc"] = np.append(
                self.stressors["efc"],
                self.stressors["efc"][-1] + self.stressors["delta_efc"][-1],
            )
            self.stressors["dod"] = np.append(
                self.stressors["dod"],
                self.stressors["dod"][1] / self.outputs["q"][-1],
            )
        else:
            self.stressors["efc"] = np.append(
                self.stressors["efc"],
                self.stressors["efc"][-1] + self.stressors["delta_efc"][-1] * self.outputs['q'][-1],
            )

        # copy end values of old stressors into a dict
        stressors = {}
        for k, v in zip(self.stressors.keys(), self.stressors.values()):
            stressors[k] = v[-1]

        self.update_states(stressors)
        self.update_outputs(stressors)

    def update_rates(self, stressors: dict):
        """
        Calculate and update battery degradation rates based on stressor values

        Updates self.rates inplace.

        Args:
            stressors (dict)    Output from extract_stressors()
        """
        pass

    def update_states(self, stressors: dict):
        """
        Update the battery states, based both on the degradation state as well as the battery performance
        at the ambient temperature, T_celsius

        Updates self.states inplace.

        Args:
            stressors (dict)    Output from extract_stressors()
        """
        pass

    def update_outputs(self, stressors: dict):
        """
        Calculate outputs, based on current battery state (and maybe stressors)

        Updates self.outputs inplace.

        Args:
            stressors (dict)    Output from extract_stressors()
        """
        pass

    def _extract_stressors(self, 
        t_secs: np.ndarray, soc: np.ndarray, T_celsius: np.ndarray
    ) -> dict:
        """Extract stressor values including time, temperature, depth-of-discharge,
        temperature, Ua, C-rate, and cycles.

        Args:
            t_secs (np.ndarray)     Array of time in seconds
            soc (np.ndarray)        Array of SOC
            T_celsius (np.ndarray)  Array of temperatures in C

        Return:
            stressors (dict)        Dictionary of stressor values.

        """
        # Extract stressors
        t_days = t_secs / (24 * 60 * 60)
        delta_t_days = t_days[-1] - t_days[0]
        delta_efc = (
            np.sum(np.abs(np.ediff1d(soc, to_begin=0))) / 2
        )  # sum the total changes to SOC / 2
        # EFC is in units of amp*hours per amp*hours nominal, not percent capacity per percent capacity nominal, so rescale by current health
        delta_efc = delta_efc * self.outputs["q"][-1]
        dod = np.max(soc) - np.min(soc)
        abs_instantaneous_crate = np.abs(
            np.diff(soc) / np.diff(t_secs / (60 * 60))
        )  # get instantaneous C-rates
        abs_instantaneous_crate[
            abs_instantaneous_crate < 1e-2
        ] = 0  # get rid of extremely small values (storage) before calculating mean
        Crate = np.trapz(abs_instantaneous_crate, t_days[1:]) / delta_t_days
        # Check storage condition, which will give nan Crate:
        if np.isnan(Crate):
            Crate = 0
        # Similar to EFC, C-rate is in units of Amps/Amp*hours nominal, not percent SOC per hour, so rescale by SOH to correct units
        Crate = Crate * self.outputs["q"][-1]
        # Temperature
        T_kelvin = T_celsius + 273.15
        # Estimate Ua (anode to reference potential) from SOC.
        # Uses the equation from Safari and Delacourt, https://doi.org/10.1149/1.3567007.
        # Anode stoichiometry is assumed to be the same for any chemistry/cell, and is calculated using the equation from Schimpe et al https://doi.org/10.1149/2.1181714jes
        # While this will not be precise, it still helps get a guess as to where the plateaus of the anode-reference potential are.
        Ua = BatteryDegradationModel._get_Ua(soc)

        cycles = count_cycles(soc)
        cycles = sum(i for _, i in cycles)

        stressors = {
            "t_secs": t_secs,
            "delta_t_days": delta_t_days,
            "delta_efc": delta_efc,
            "TdegK": T_kelvin,
            "soc": soc,
            "dod": dod,
            "Ua": Ua,
            "Crate": Crate,
            "cycles": cycles,
        }

        return stressors


"""
# Test of 'find_breakpoints'
from .rainflow import extract_cycles

hours = 6 * 24
t_hours = np.arange(hours + 1)
t_secs = t_hours * 3600
soc = np.tile(
    np.array(
        [1, 0.6, 0.4, 0.4, 0.6, 1,
        1, 0.5, 0, 0.5, 1, 1,
        1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1,
        0.98, 1, 0.98, 1, 0.98, 1,
        0.98, 1, 0.98, 1, 0.98, 1,
        0.98, 1, 0.98, 1, 0.98, 1,
        0.98, 1, 0.98, 1, 0.98, 1,
        0.98, 1, 0.98, 1, 0.98, 1,
        0.98, 1, 0.98, 1, 0.98, 1,
        0.98, 1, 0.98, 1, 0.98, 1,
        0.98, 1, 0.98, 1, 0.98, 1,]
        ),
    2
)
soc = np.append(soc, 1)
EFCs = np.cumsum(np.abs(np.ediff1d(soc, to_begin=0)))/2

# Chunk the input timeseries into cycles (SOC turning points) using the rainflow algorithm
cycles_generator = (
        [span, mean, count, idx_start, idx_end] for span, mean, count, idx_start, idx_end in extract_cycles(soc)
    )
idx_turning_point = []
for cycle in cycles_generator:
    idx_turning_point.append(cycle[-1])
idx_turning_point = np.array(idx_turning_point)

# Find the breakpoints for battery life simulation
idx_simulation_breakpoints = find_breakpoints(t_secs, EFCs, idx_turning_point)

plt.plot(t_hours/24, soc)
plt.plot(t_hours[idx_turning_point]/24, soc[idx_turning_point], '.g')
plt.plot(t_hours[idx_simulation_breakpoints]/24, soc[idx_simulation_breakpoints], 'or')
"""
