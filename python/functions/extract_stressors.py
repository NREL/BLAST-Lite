import numpy as np
import functions.rainflow as rainflow

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
    def get_Xa(soc):
        return 8.5*10**-3 + soc*(0.78 - 8.5*10**-3)
    def get_Ua(Xa):
        return (0.6379 + 0.5416*np.exp(-305.5309*Xa) + 
            0.044*np.tanh(-1*(Xa-0.1958)/0.1088) - 0.1978*np.tanh((Xa-1.0571)/0.0854) - 
            0.6875*np.tanh((Xa+0.0117)/0.0529) - 0.0175*np.tanh((Xa-0.5692)/0.0875))
    Ua = get_Ua(get_Xa(soc))

    cycles = rainflow.count_cycles(soc)
    cycles = sum(i for _, i in cycles)

    return delta_t_days, delta_efc, T_kelvin, soc, Ua, dod, Crate, cycles