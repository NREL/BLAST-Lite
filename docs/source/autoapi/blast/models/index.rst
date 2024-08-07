blast.models
============

.. py:module:: blast.models


Submodules
----------

.. toctree::
   :maxdepth: 1

   /autoapi/blast/models/degradation_model/index
   /autoapi/blast/models/lfp_gr_250AhPrismatic_2019/index
   /autoapi/blast/models/lfp_gr_SonyMurata3Ah_2018/index
   /autoapi/blast/models/lmo_gr_NissanLeaf66Ah_2ndLife_2020/index
   /autoapi/blast/models/nca_gr_Panasonic3Ah_2018/index
   /autoapi/blast/models/nca_grsi_SonyMurata2p5Ah_2023/index
   /autoapi/blast/models/nmc111_gr_Kokam75Ah_2017/index
   /autoapi/blast/models/nmc111_gr_Sanyo2Ah_2014/index
   /autoapi/blast/models/nmc811_grSi_LGMJ1_4Ah_2020/index
   /autoapi/blast/models/nmc_gr_50Ah_B1_2020/index
   /autoapi/blast/models/nmc_gr_50Ah_B2_2020/index
   /autoapi/blast/models/nmc_gr_75Ah_A_2019/index
   /autoapi/blast/models/nmc_lto_10Ah_2020/index


Classes
-------

.. autoapisummary::

   blast.models.BatteryDegradationModel
   blast.models.Lfp_Gr_250AhPrismatic
   blast.models.Lfp_Gr_SonyMurata3Ah_Battery
   blast.models.Lmo_Gr_NissanLeaf66Ah_2ndLife_Battery
   blast.models.NCA_GrSi_SonyMurata2p5Ah_Battery
   blast.models.NMC_Gr_50Ah_B1
   blast.models.NMC_Gr_50Ah_B2
   blast.models.NMC_Gr_75Ah_A
   blast.models.Nca_Gr_Panasonic3Ah_Battery
   blast.models.Nmc111_Gr_Kokam75Ah_Battery
   blast.models.Nmc111_Gr_Sanyo2Ah_Battery
   blast.models.Nmc811_GrSi_LGMJ1_4Ah_Battery
   blast.models.Nmc_Lto_10Ah_Battery


Package Contents
----------------

.. py:class:: BatteryDegradationModel

   .. py:method:: simulate_battery_life(input_timeseries, simulation_years = None, is_constant_input = False, breakpoints_max_time_diff_s = 86400, breakpoints_max_EFC_diff = 1)

      Run battery life simulation over the input, or repeat for the number of years specified.

      Updates attributes self.rates, self.stressors, self.outputs, and self.states inplace.

      :param input_timeseries:
      :type input_timeseries: dict, pd.DataFrame
      :param simulation_years:
      :type simulation_years: float
      :param is_constant_input:
      :type is_constant_input: bool
      :param breakpoints_max_time_diff_s:
      :type breakpoints_max_time_diff_s: float
      :param breakpoints_max_EFC_diff:
      :type breakpoints_max_EFC_diff: float



   .. py:method:: update_battery_state(t_secs, soc, T_celsius)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.

      :param t_secs: for the soc_timeseries data points
      :type t_secs: np.ndarray
      :param soc: Vector of the state-of-charge of the battery at each t_sec
      :type soc: np.ndarray
      :param T_celsius:
      :type T_celsius: ndarray



   .. py:method:: update_battery_state_repeating()

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
      degradation rates are unchanging for every timestep, and don't need to be calculated again.

      Updates self.states and self.outputs inplace.



   .. py:method:: update_outputs(stressors)

      Calculate outputs, based on current battery state (and maybe stressors)

      Updates self.outputs inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_rates(stressors)

      Calculate and update battery degradation rates based on stressor values

      Updates self.rates inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_states(stressors)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius

      Updates self.states inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:attribute:: experimental_range

      details on the range of experimental conditions, i.e.,
      the range we expect the model to be valid in

      :type: dict


   .. py:attribute:: outputs

      Battery properties derived from state values

      :type: dict


   .. py:attribute:: rates

      History of stressor-dependent degradation rates

      :type: dict


   .. py:attribute:: states

      Internal states of the battery model

      :type: dict


   .. py:attribute:: stressors

      History of stressors on the battery

      :type: dict


.. py:class:: Lfp_Gr_250AhPrismatic(degradation_scalar = 1, label = 'LFP-Gr 250Ah')

   Bases: :py:obj:`blast.models.degradation_model.BatteryDegradationModel`


   Model predicting the degradation of large-format prismatic commercial LFP-Gr
   cells from a large manufacturer, with >250 Ah capacity and an
   energy-to-power ratio of 6 h**-1 (high energy density cell, low power).
   Experimental test data is reported in https://doi.org/10.1016/j.est.2023.109042.

   .. note::
       EXPERIMENTAL AGING DATA SUMMARY:
           Experimental test data is reported in https://doi.org/10.1016/j.est.2023.109042.
           Aging test matrix varied temperature and state-of-charge for calendar aging, and
           varied depth-of-discharge, average state-of-charge, and C-rates for cycle aging.
           Charging rate was limited to a maximum of 0.16 C at 10 Celsius, and 0.65 C at all
           higher temperatures, so the model is only fit on low rate charge data.

       MODEL SENSITIVITY
           The model predicts degradation rate versus time as a function of temperature and average
           state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as
           a function of average state-of-charge during a cycle, depth-of-discharge, and average of the
           charge and discharge C-rates.
           Test data showed that the cycling degradation rate was nearly identical for all cells,
           despite varying temperature, average SOC, and dis/charge rates. This means that the cycle aging
           model learned the inverse temperature depedence of the calendar aging model to 'balance' the
           the degradation rate at the tested cycling temperatures (10 Celsius to 45 Celsius). So, the model
           will likely make very poor extrapolations outside of this tested temperature range for cycling fade.

       MODEL LIMITATIONS
           Charging and discharging rates were very conservative, following cycling protocols as suggested by the
           cell manufacturer. Degradation at higher charging or discharging rates will not be simulated accurately.


   .. py:method:: simulate_battery_life(input_timeseries, simulation_years = None, is_constant_input = False, breakpoints_max_time_diff_s = 86400, breakpoints_max_EFC_diff = 1)

      Run battery life simulation over the input, or repeat for the number of years specified.

      Updates attributes self.rates, self.stressors, self.outputs, and self.states inplace.

      :param input_timeseries:
      :type input_timeseries: dict, pd.DataFrame
      :param simulation_years:
      :type simulation_years: float
      :param is_constant_input:
      :type is_constant_input: bool
      :param breakpoints_max_time_diff_s:
      :type breakpoints_max_time_diff_s: float
      :param breakpoints_max_EFC_diff:
      :type breakpoints_max_EFC_diff: float



   .. py:method:: update_battery_state(t_secs, soc, T_celsius)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.

      :param t_secs: for the soc_timeseries data points
      :type t_secs: np.ndarray
      :param soc: Vector of the state-of-charge of the battery at each t_sec
      :type soc: np.ndarray
      :param T_celsius:
      :type T_celsius: ndarray



   .. py:method:: update_battery_state_repeating()

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
      degradation rates are unchanging for every timestep, and don't need to be calculated again.

      Updates self.states and self.outputs inplace.



   .. py:method:: update_outputs(stressors)

      Calculate outputs, based on current battery state (and maybe stressors)

      Updates self.outputs inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_rates(stressors)

      Calculate and update battery degradation rates based on stressor values

      Updates self.rates inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_states(stressors)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius

      Updates self.states inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:property:: cap


   .. py:attribute:: experimental_range


   .. py:attribute:: outputs


   .. py:attribute:: rates


   .. py:attribute:: states


   .. py:attribute:: stressors


.. py:class:: Lfp_Gr_SonyMurata3Ah_Battery(degradation_scalar = 1, label = 'LFP-Gr Sony-Murata')

   Bases: :py:obj:`blast.models.degradation_model.BatteryDegradationModel`


   Model predicting the degradation of Sony-Murata 3 Ah LFP-Gr cylindrical cells.
   Data is from Technical University of Munich, reported in studies led by Maik Naumann.
   Capacity model identification was conducted at NREL. Resistance model is from Naumann et al.
   Naumann et al used an interative fitting procedure, but it was found that lower model error could be
   achieved by simply reoptimizing all resistance growth parameters to the entire data set.
   Calendar aging data source: https://doi.org/10.1016/j.est.2018.01.019
   Cycle aging data source: https://doi.org/10.1016/j.jpowsour.2019.227666
   Model identification source: https://doi.org/10.1149/1945-7111/ac86a8

   Degradation rate is a function of the aging stressors, i.e., ambient temperature and use.
   The state of the battery is updated throughout the lifetime of the cell.
   Performance metrics are capacity and DC resistance. These metrics change as a function of the
   cell's current degradation state, as well as the ambient temperature. The model predicts time and
   cycling dependent degradation. Cycling dependent degradation includes a break-in mechanism as well
   as long term cycling fade; the break-in mechanism strongly influenced results of the accelerated
   aging test, but is not expected to have much influence on real-world applications.
   Parameters to modify to change fade rates:
   - q1_b0: rate of capacity loss due to calendar degradation
   - q5_b0: rate of capacity loss due to cycling degradation
   - k_ref_r_cal: rate of resistance growth due to calendar degradation
   - A_r_cyc: rate of resistance growth due to cycling degradation

   .. note::
       EXPERIMENTAL AGING DATA SUMMARY:
           Aging test matrix varied temperature and state-of-charge for calendar aging, and
           varied depth-of-discharge, average state-of-charge, and C-rates for cycle aging.
           There is NO LOW TEMPERATURE cycling aging data, i.e., no lithium-plating induced by
           kinetic limitations on cell performance; CYCLING WAS ONLY DONE AT 25 CELSIUS AND 45 CELSIUS,
           so any model predictions at low temperature cannot incorporate low temperature degradation modes.
           Discharge capacity

       MODEL SENSITIVITY
           The model predicts degradation rate versus time as a function of temperature and average
           state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as
           a function of average state-of-charge during a cycle, depth-of-discharge, and average of the
           charge and discharge C-rates.

       MODEL LIMITATIONS
           There is no influence of TEMPERATURE on CYCLING DEGRADATION RATE due to limited data. This is
           NOT PHYSICALLY REALISTIC AND IS BASED ON LIMITED DATA.


   .. py:method:: simulate_battery_life(input_timeseries, simulation_years = None, is_constant_input = False, breakpoints_max_time_diff_s = 86400, breakpoints_max_EFC_diff = 1)

      Run battery life simulation over the input, or repeat for the number of years specified.

      Updates attributes self.rates, self.stressors, self.outputs, and self.states inplace.

      :param input_timeseries:
      :type input_timeseries: dict, pd.DataFrame
      :param simulation_years:
      :type simulation_years: float
      :param is_constant_input:
      :type is_constant_input: bool
      :param breakpoints_max_time_diff_s:
      :type breakpoints_max_time_diff_s: float
      :param breakpoints_max_EFC_diff:
      :type breakpoints_max_EFC_diff: float



   .. py:method:: update_battery_state(t_secs, soc, T_celsius)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.

      :param t_secs: for the soc_timeseries data points
      :type t_secs: np.ndarray
      :param soc: Vector of the state-of-charge of the battery at each t_sec
      :type soc: np.ndarray
      :param T_celsius:
      :type T_celsius: ndarray



   .. py:method:: update_battery_state_repeating()

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
      degradation rates are unchanging for every timestep, and don't need to be calculated again.

      Updates self.states and self.outputs inplace.



   .. py:method:: update_outputs(stressors)

      Calculate outputs, based on current battery state (and maybe stressors)

      Updates self.outputs inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_rates(stressors)

      Calculate and update battery degradation rates based on stressor values

      Updates self.rates inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_states(stressors)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius

      Updates self.states inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:property:: cap


   .. py:attribute:: experimental_range


   .. py:attribute:: outputs


   .. py:attribute:: rates


   .. py:attribute:: states


   .. py:attribute:: stressors


.. py:class:: Lmo_Gr_NissanLeaf66Ah_2ndLife_Battery(degradation_scalar = 1, label = 'LMO-Gr Nissan Leaf')

   Bases: :py:obj:`blast.models.degradation_model.BatteryDegradationModel`


   Model fit to SECOND LIFE data on Nissan Leaf half-modules (2p cells) by Braco et al.
   https://doi.org/10.1109/EEEIC/ICPSEUROPE54979.2022.9854784 (calendar aging data)
   https://doi.org/10.1016/j.est.2020.101695 (cycle aging data)
   Note that these cells are already hugely degraded, starting out at an average relative capacity
   of 70%. So the model reports q and qNew, where qNew is relative to initial.

   .. note ::
       EXPERIMENTAL AGING DATA SUMMARY:
           Calendar aging widely varied SOC and temperature.
           Cycle aging is only at a single condition (25 Celsius, 100% DOD, 1C-1C).

       MODEL SENSITIVITY
           The model predicts degradation rate versus time as a function of temperature and average
           state-of-charge and degradation rate is only a function equivalent full cycles.

       MODEL LIMITATIONS
           Cycling degradation IS ONLY A FUNCTION OF CHARGE THROUGHPUT due to limited aging data.
           Cycling degradation predictions ARE ONLY VALID NEAR 25 CELSIUS, 100% DOD, 1 C CHARGE/DISCHARGE RATE.


   .. py:method:: simulate_battery_life(input_timeseries, simulation_years = None, is_constant_input = False, breakpoints_max_time_diff_s = 86400, breakpoints_max_EFC_diff = 1)

      Run battery life simulation over the input, or repeat for the number of years specified.

      Updates attributes self.rates, self.stressors, self.outputs, and self.states inplace.

      :param input_timeseries:
      :type input_timeseries: dict, pd.DataFrame
      :param simulation_years:
      :type simulation_years: float
      :param is_constant_input:
      :type is_constant_input: bool
      :param breakpoints_max_time_diff_s:
      :type breakpoints_max_time_diff_s: float
      :param breakpoints_max_EFC_diff:
      :type breakpoints_max_EFC_diff: float



   .. py:method:: update_battery_state(t_secs, soc, T_celsius)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.

      :param t_secs: for the soc_timeseries data points
      :type t_secs: np.ndarray
      :param soc: Vector of the state-of-charge of the battery at each t_sec
      :type soc: np.ndarray
      :param T_celsius:
      :type T_celsius: ndarray



   .. py:method:: update_battery_state_repeating()

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
      degradation rates are unchanging for every timestep, and don't need to be calculated again.

      Updates self.states and self.outputs inplace.



   .. py:method:: update_outputs(stressors)

      Calculate outputs, based on current battery state (and maybe stressors)

      Updates self.outputs inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_rates(stressors)

      Calculate and update battery degradation rates based on stressor values

      Updates self.rates inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_states(stressors)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius

      Updates self.states inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:property:: cap


   .. py:property:: cap_2ndLife


   .. py:attribute:: experimental_range


   .. py:attribute:: outputs


   .. py:attribute:: rates


   .. py:attribute:: states


   .. py:attribute:: stressors


.. py:class:: NCA_GrSi_SonyMurata2p5Ah_Battery(degradation_scalar = 1, label = 'NCA-GrSi Sony-Murata')

   Bases: :py:obj:`blast.models.degradation_model.BatteryDegradationModel`


   Model predicting the degradation of Sony-Murata US18650VTC5A 3.5 Ah NCA-GrSi cylindrical cells.
   Relatively high power cells, with 1.4 wt% Si in the graphite-si composite negative electrode.
   (Maximum continuous charge rate of 2C and max continuous discharge rate of 10C, even at 5 degC)
   Data is from Technical University of Munich, reported in studies led by Leo Wildfeuer and Alexander Karger.
   Accelerated aging test data reported in https://doi.org/10.1016/j.jpowsour.2022.232498
   The model here was identifed using AI-Batt at NREL. Note that the TUM authors have put more effort
   into model identification on this data set, see the following papers for more detailed models identified on this data set:
   Mechanistic cycle aging model (LLI, LAM_PE, LAM_(NE,Gr), LAM_(NE,Si)): https://doi.org/10.1016/j.jpowsour.2023.233947
   Mechanistic calendar aging model (considers impact of capacity check frequency): https://doi.org/10.1016/j.jpowsour.2023.233208

   .. note::
       EXPERIMENTAL AGING DATA SUMMARY:
           Aging test matrix varied temperature and state-of-charge for calendar aging, and
           varied depth-of-discharge, average state-of-charge, and C-rates for cycle aging.
           For calendar aging, in addition to temperature and SOC, capacity-check frequency was also varied.

       MODEL SENSITIVITY
           The model predicts degradation rate versus time as a function of temperature and average
           state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as
           a function of average state-of-charge during a cycle, depth-of-discharge, and average of the
           charge and discharge C-rates.

       MODEL LIMITATIONS
           I did not model the impact of capacity check frequency, only using the 6 week capacity check data.
           See https://doi.org/10.1016/j.jpowsour.2023.233208 for a detailed consideration of the capacity check frequency impact on aging.
           Only one cycling cell in the data shows 'knee-over' behavior, making empirical model identification of this behavior challenging.
           I simply neglected the knee over; this knee occurs at ~80% capacity fade, so note that predictions of degradation at maximum
           DOD and below 80% capacity should not be believed.


   .. py:method:: simulate_battery_life(input_timeseries, simulation_years = None, is_constant_input = False, breakpoints_max_time_diff_s = 86400, breakpoints_max_EFC_diff = 1)

      Run battery life simulation over the input, or repeat for the number of years specified.

      Updates attributes self.rates, self.stressors, self.outputs, and self.states inplace.

      :param input_timeseries:
      :type input_timeseries: dict, pd.DataFrame
      :param simulation_years:
      :type simulation_years: float
      :param is_constant_input:
      :type is_constant_input: bool
      :param breakpoints_max_time_diff_s:
      :type breakpoints_max_time_diff_s: float
      :param breakpoints_max_EFC_diff:
      :type breakpoints_max_EFC_diff: float



   .. py:method:: update_battery_state(t_secs, soc, T_celsius)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.

      :param t_secs: for the soc_timeseries data points
      :type t_secs: np.ndarray
      :param soc: Vector of the state-of-charge of the battery at each t_sec
      :type soc: np.ndarray
      :param T_celsius:
      :type T_celsius: ndarray



   .. py:method:: update_battery_state_repeating()

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
      degradation rates are unchanging for every timestep, and don't need to be calculated again.

      Updates self.states and self.outputs inplace.



   .. py:method:: update_outputs(stressors)

      Calculate outputs, based on current battery state (and maybe stressors)

      Updates self.outputs inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_rates(stressors)

      Calculate and update battery degradation rates based on stressor values

      Updates self.rates inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_states(stressors)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius

      Updates self.states inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:property:: cap


   .. py:attribute:: experimental_range


   .. py:attribute:: outputs


   .. py:attribute:: rates


   .. py:attribute:: states


   .. py:attribute:: stressors


.. py:class:: NMC_Gr_50Ah_B1(degradation_scalar = 1, label = 'NMC-Gr B1 50Ah')

   Bases: :py:obj:`blast.models.degradation_model.BatteryDegradationModel`


   Model predicting the degradation of large format pouch 'NMC-Gr B1' cells from a
   large manufacturer, with ~50 Ah capacity and an
   energy-to-power ratio of 16 h^-1 (relatively high power, suitable for fast charging).
   Experimental test data is reported in https://doi.org/10.1016/j.est.2023.109042.
   This is for the 'NMC-Gr B1' cell reported in the paper.

   .. note::
       EXPERIMENTAL AGING DATA SUMMARY:
           Experimental test data is reported in https://doi.org/10.1016/j.est.2023.109042.
           Aging test matrix varied temperature and state-of-charge for calendar aging, and
           varied depth-of-discharge, average state-of-charge, and C-rates for cycle aging.

       MODEL SENSITIVITY
           The model predicts degradation rate versus time as a function of temperature and average
           state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as
           a function of average state-of-charge during a cycle, depth-of-discharge, and average of the
           charge and discharge C-rates.

       MODEL LIMITATIONS
           Charging and discharging rates were conservative, following cycling protocols as suggested by the
           cell manufacturer. Degradation at higher charging or discharging rates will not be simulated accurately.
           Charging at 10 degC was limited to a very conservative rate of 0.3C by the manufacturer. Simulations with
           low temperatures and charging rates above C/3 will likely give inaccurate predictions.


   .. py:method:: simulate_battery_life(input_timeseries, simulation_years = None, is_constant_input = False, breakpoints_max_time_diff_s = 86400, breakpoints_max_EFC_diff = 1)

      Run battery life simulation over the input, or repeat for the number of years specified.

      Updates attributes self.rates, self.stressors, self.outputs, and self.states inplace.

      :param input_timeseries:
      :type input_timeseries: dict, pd.DataFrame
      :param simulation_years:
      :type simulation_years: float
      :param is_constant_input:
      :type is_constant_input: bool
      :param breakpoints_max_time_diff_s:
      :type breakpoints_max_time_diff_s: float
      :param breakpoints_max_EFC_diff:
      :type breakpoints_max_EFC_diff: float



   .. py:method:: update_battery_state(t_secs, soc, T_celsius)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.

      :param t_secs: for the soc_timeseries data points
      :type t_secs: np.ndarray
      :param soc: Vector of the state-of-charge of the battery at each t_sec
      :type soc: np.ndarray
      :param T_celsius:
      :type T_celsius: ndarray



   .. py:method:: update_battery_state_repeating()

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
      degradation rates are unchanging for every timestep, and don't need to be calculated again.

      Updates self.states and self.outputs inplace.



   .. py:method:: update_outputs(stressors)

      Calculate outputs, based on current battery state (and maybe stressors)

      Updates self.outputs inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_rates(stressors)

      Calculate and update battery degradation rates based on stressor values

      Updates self.rates inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_states(stressors)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius

      Updates self.states inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:property:: cap


   .. py:attribute:: experimental_range


   .. py:attribute:: outputs


   .. py:attribute:: rates


   .. py:attribute:: states


   .. py:attribute:: stressors


.. py:class:: NMC_Gr_50Ah_B2(degradation_scalar = 1, label = 'NMC-Gr B2 50Ah')

   Bases: :py:obj:`blast.models.degradation_model.BatteryDegradationModel`


   Model predicting the degradation of Large format pouch 'NMC-Gr B2' cells from a large manufacturer,
   with ~50 Ah capacity and an energy-to-power ratio of 14 h^-1 (balance of energy and power,
   questionable durability/performance for fast charging).
   Experimental test data is reported in https://doi.org/10.1016/j.est.2023.109042.
   This is for the 'NMC-Gr B2' cell reported in the paper.

   .. note::
       EXPERIMENTAL AGING DATA SUMMARY:
           Experimental test data is reported in https://doi.org/10.1016/j.est.2023.109042.
           Aging test matrix varied temperature and state-of-charge for calendar aging, and
           varied depth-of-discharge, average state-of-charge, and C-rates for cycle aging.

       MODEL SENSITIVITY
           The model predicts degradation rate versus time as a function of temperature and average
           state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as
           a function of average state-of-charge during a cycle, depth-of-discharge, and average of the
           charge and discharge C-rates.

       MODEL LIMITATIONS
           Charging and discharging rates were conservative, following cycling protocols as suggested by the
           cell manufacturer. Degradation at higher charging or discharging rates will not be simulated accurately.
           Charging at 10 degC was limited to a very conservative rate of 0.3C by the manufacturer. Simulations with
           low temperatures and charging rates above C/3 will likely give inaccurate predictions.


   .. py:method:: simulate_battery_life(input_timeseries, simulation_years = None, is_constant_input = False, breakpoints_max_time_diff_s = 86400, breakpoints_max_EFC_diff = 1)

      Run battery life simulation over the input, or repeat for the number of years specified.

      Updates attributes self.rates, self.stressors, self.outputs, and self.states inplace.

      :param input_timeseries:
      :type input_timeseries: dict, pd.DataFrame
      :param simulation_years:
      :type simulation_years: float
      :param is_constant_input:
      :type is_constant_input: bool
      :param breakpoints_max_time_diff_s:
      :type breakpoints_max_time_diff_s: float
      :param breakpoints_max_EFC_diff:
      :type breakpoints_max_EFC_diff: float



   .. py:method:: update_battery_state(t_secs, soc, T_celsius)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.

      :param t_secs: for the soc_timeseries data points
      :type t_secs: np.ndarray
      :param soc: Vector of the state-of-charge of the battery at each t_sec
      :type soc: np.ndarray
      :param T_celsius:
      :type T_celsius: ndarray



   .. py:method:: update_battery_state_repeating()

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
      degradation rates are unchanging for every timestep, and don't need to be calculated again.

      Updates self.states and self.outputs inplace.



   .. py:method:: update_outputs(stressors)

      Calculate outputs, based on current battery state (and maybe stressors)

      Updates self.outputs inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_rates(stressors)

      Calculate and update battery degradation rates based on stressor values

      Updates self.rates inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_states(stressors)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius

      Updates self.states inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:property:: cap


   .. py:attribute:: experimental_range


   .. py:attribute:: outputs


   .. py:attribute:: rates


   .. py:attribute:: states


   .. py:attribute:: stressors


.. py:class:: NMC_Gr_75Ah_A(degradation_scalar = 1, label = 'NMC-Gr A 75Ah')

   Bases: :py:obj:`blast.models.degradation_model.BatteryDegradationModel`


   Model predicting the degradation of large format pouch 'NMC-Gr A' cells from a
   large manufacturer, with ~75 Ah capacity and an energy-to-power ratio of 16 h^-1
   (relatively high power, suitable for fast charging).
   Experimental test data is reported in https://doi.org/10.1016/j.est.2023.109042.
   This is for the 'NMC-Gr A' cell reported in the paper.

   .. note::
       EXPERIMENTAL AGING DATA SUMMARY:
           Experimental test data is reported in https://doi.org/10.1016/j.est.2023.109042.
           Aging test matrix varied temperature and state-of-charge for calendar aging, and
           varied depth-of-discharge, average state-of-charge, and C-rates for cycle aging.

       MODEL SENSITIVITY
           The model predicts degradation rate versus time as a function of temperature and average
           state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as
           a function of average state-of-charge during a cycle, depth-of-discharge, and average of the
           charge and discharge C-rates.

       MODEL LIMITATIONS
           Charging and discharging rates were conservative, following cycling protocols as suggested by the
           cell manufacturer. Degradation at higher charging or discharging rates will not be simulated accurately.


   .. py:method:: simulate_battery_life(input_timeseries, simulation_years = None, is_constant_input = False, breakpoints_max_time_diff_s = 86400, breakpoints_max_EFC_diff = 1)

      Run battery life simulation over the input, or repeat for the number of years specified.

      Updates attributes self.rates, self.stressors, self.outputs, and self.states inplace.

      :param input_timeseries:
      :type input_timeseries: dict, pd.DataFrame
      :param simulation_years:
      :type simulation_years: float
      :param is_constant_input:
      :type is_constant_input: bool
      :param breakpoints_max_time_diff_s:
      :type breakpoints_max_time_diff_s: float
      :param breakpoints_max_EFC_diff:
      :type breakpoints_max_EFC_diff: float



   .. py:method:: update_battery_state(t_secs, soc, T_celsius)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.

      :param t_secs: for the soc_timeseries data points
      :type t_secs: np.ndarray
      :param soc: Vector of the state-of-charge of the battery at each t_sec
      :type soc: np.ndarray
      :param T_celsius:
      :type T_celsius: ndarray



   .. py:method:: update_battery_state_repeating()

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
      degradation rates are unchanging for every timestep, and don't need to be calculated again.

      Updates self.states and self.outputs inplace.



   .. py:method:: update_outputs(stressors)

      Calculate outputs, based on current battery state (and maybe stressors)

      Updates self.outputs inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_rates(stressors)

      Calculate and update battery degradation rates based on stressor values

      Updates self.rates inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_states(stressors)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius

      Updates self.states inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:property:: cap


   .. py:attribute:: experimental_range


   .. py:attribute:: outputs


   .. py:attribute:: rates


   .. py:attribute:: states


   .. py:attribute:: stressors


.. py:class:: Nca_Gr_Panasonic3Ah_Battery(degradation_scalar = 1, label = 'NCA-Gr Panasonic')

   Bases: :py:obj:`blast.models.degradation_model.BatteryDegradationModel`


   Model fit to Panasonic 18650B NCA-Gr cells. High-ish energy density 18650 cells with
   adequate lifetime.
   Calendar data is reported by Keil et al (https://dx.doi.org/10.1149/2.0411609jes)
   Cycling data is reported by Preger et al (https://doi.org/10.1149/1945-7111/abae37) and
   is available at batteryarchive.com.
   The authors of BLAST-Lite are not aware of any study conducting both calendar aging and
   cycle aging of these cells.

   .. note::
       EXPERIMENTAL AGING DATA SUMMARY:
           Calendar aging widely varied SOC at 25, 40, and 50 Celsius. 300 days max.
           Cycle aging varied temperature and C-rates, and DOD. Some accelerating fade is observed
           at room temperature and high DODs but isn't modeled well here. That's not a huge problem,
           because the modeled lifetime is quite short anyways.

       MODEL SENSITIVITY
           The model predicts degradation rate versus time as a function of temperature and average
           state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as
           a function of C-rate, temperature, and depth-of-discharge (DOD dependence is assumed to be linear, no aging data)

       MODEL LIMITATIONS
           Cycle degradation predictions WILL NOT PREDICT KNEE-POINT due to limited data.
           Cycle aging is only modeled at 25, 35, and 45 Celsius, PREDICTIONS OUTSIDE THIS
           TEMPERATURE RANGE MAY BE OPTIMISTIC.


   .. py:method:: simulate_battery_life(input_timeseries, simulation_years = None, is_constant_input = False, breakpoints_max_time_diff_s = 86400, breakpoints_max_EFC_diff = 1)

      Run battery life simulation over the input, or repeat for the number of years specified.

      Updates attributes self.rates, self.stressors, self.outputs, and self.states inplace.

      :param input_timeseries:
      :type input_timeseries: dict, pd.DataFrame
      :param simulation_years:
      :type simulation_years: float
      :param is_constant_input:
      :type is_constant_input: bool
      :param breakpoints_max_time_diff_s:
      :type breakpoints_max_time_diff_s: float
      :param breakpoints_max_EFC_diff:
      :type breakpoints_max_EFC_diff: float



   .. py:method:: update_battery_state(t_secs, soc, T_celsius)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.

      :param t_secs: for the soc_timeseries data points
      :type t_secs: np.ndarray
      :param soc: Vector of the state-of-charge of the battery at each t_sec
      :type soc: np.ndarray
      :param T_celsius:
      :type T_celsius: ndarray



   .. py:method:: update_battery_state_repeating()

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
      degradation rates are unchanging for every timestep, and don't need to be calculated again.

      Updates self.states and self.outputs inplace.



   .. py:method:: update_outputs(stressors)

      Calculate outputs, based on current battery state (and maybe stressors)

      Updates self.outputs inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_rates(stressors)

      Calculate and update battery degradation rates based on stressor values

      Updates self.rates inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_states(stressors)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius

      Updates self.states inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:property:: cap


   .. py:attribute:: experimental_range


   .. py:attribute:: outputs


   .. py:attribute:: rates


   .. py:attribute:: states


   .. py:attribute:: stressors


.. py:class:: Nmc111_Gr_Kokam75Ah_Battery(degradation_scalar = 1, label = 'NMC111-Gr Kokam')

   Bases: :py:obj:`blast.models.degradation_model.BatteryDegradationModel`


   Model predicting the degradation of a Kokam 75 Ah NMC-Gr pouch cell.
   https://ieeexplore.ieee.org/iel7/7951530/7962914/07963578.pdf
   It is uncertain if the exact NMC composition is 1-1-1, but it this is definitely not a high nickel (>80%) cell.
   Degradation rate is a function of the aging stressors, i.e., ambient temperature and use.
   The state of the battery is updated throughout the lifetime of the cell.
   Performance metrics are capacity and DC resistance. These metrics change as a function of the
   cell's current degradation state, as well as the ambient temperature. The model predicts time and
   cycling dependent degradation, using Loss of Lithium Inventory (LLI) and Loss of Active
   Material (LAM) degradation modes that interact competitively (cell performance is limited by
   one or the other.)
   Parameters to modify to change fade rates:
   - Calendar capacity loss rate: q1_0
   - Cycling capacity loss rate (LLI): q3_0
   - Cycling capacity loss rate (LAM): q5_0, will also effect resistance growth onset due to LAM.
   - Calendar resistance growth rate (LLI), relative to capacity loss rate: r1
   - Cycling resistance growth rate (LLI), relative to capacity loss rate: r3

   .. note::
       EXPERIMENTAL AGING DATA SUMMARY:
           Aging test matrix varied primarly temperature, with small DOD variation.
           Calendar and cycle aging were performed between 0 and 55 Celsius. C-rates always at 1C,
           except for charging at 0 Celsius, which was conducted at C/3. Depth-of-discharge was 80%
           for nearly all tests (3.4 V - 4.1 V), with one 100% DOD test (3 V - 4.2 V).
           Reported relative capacity was measured at C/5 rate at the aging temperatures. Reported
           relative DC resistance was measured by HPPC using a 10s, 1C DC pulse, averaged between
           charge and discharge, calculated using a simple ohmic fit of the voltage response.

       MODEL SENSITIVITY
           The model predicts degradation rate versus time as a function of temperature and average
           state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as
           a function of temperature and depth-of-discharge. Sensitivity to cycling degradation rate
           at low temperature is inferred from physical insight due to limited data.

       MODEL LIMITATIONS
           There is NO C-RATE DEPENDENCE for degradation in this model. THIS IS NOT PHYSICALLY REALISTIC
           AND IS BASED ON LIMITED DATA.


   .. py:method:: simulate_battery_life(input_timeseries, simulation_years = None, is_constant_input = False, breakpoints_max_time_diff_s = 86400, breakpoints_max_EFC_diff = 1)

      Run battery life simulation over the input, or repeat for the number of years specified.

      Updates attributes self.rates, self.stressors, self.outputs, and self.states inplace.

      :param input_timeseries:
      :type input_timeseries: dict, pd.DataFrame
      :param simulation_years:
      :type simulation_years: float
      :param is_constant_input:
      :type is_constant_input: bool
      :param breakpoints_max_time_diff_s:
      :type breakpoints_max_time_diff_s: float
      :param breakpoints_max_EFC_diff:
      :type breakpoints_max_EFC_diff: float



   .. py:method:: update_battery_state(t_secs, soc, T_celsius)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.

      :param t_secs: for the soc_timeseries data points
      :type t_secs: np.ndarray
      :param soc: Vector of the state-of-charge of the battery at each t_sec
      :type soc: np.ndarray
      :param T_celsius:
      :type T_celsius: ndarray



   .. py:method:: update_battery_state_repeating()

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
      degradation rates are unchanging for every timestep, and don't need to be calculated again.

      Updates self.states and self.outputs inplace.



   .. py:method:: update_outputs(stressors)

      Calculate outputs, based on current battery state (and maybe stressors)

      Updates self.outputs inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_rates(stressors)

      Calculate and update battery degradation rates based on stressor values

      Updates self.rates inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_states(stressors)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius

      Updates self.states inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:property:: cap


   .. py:attribute:: experimental_range


   .. py:attribute:: outputs


   .. py:attribute:: rates


   .. py:attribute:: states


   .. py:attribute:: stressors


.. py:class:: Nmc111_Gr_Sanyo2Ah_Battery(degradation_scalar = 1, label = 'NMC111-Gr Sanyo')

   Bases: :py:obj:`blast.models.degradation_model.BatteryDegradationModel`


   Model predicting the degradation of Sanyo UR18650E cells, published by Schmalsteig et al:
   http://dx.doi.org/10.1016/j.jpowsour.2014.02.012.
   More detailed analysis of cell performance and voltage vs. state-of-charge data was copied from
   Ecker et al: http://dx.doi.org/10.1016/j.jpowsour.2013.09.143 (KNEE POINTS OBSERVED IN ECKER ET AL
   AT HIGH DEPTH OF DISCHARGE WERE SIMPLY NOT ADDRESSED DURING MODEL FITTING BY SCHMALSTEIG ET AL).
   Voltage lookup table here use data from Ecker et al for 0/10% SOC, and other values were extracted
   from Figure 1 in Schmalsteig et al using WebPlotDigitizer.

   .. note::
       EXPERIMENTAL AGING DATA SUMMARY:
           Calendar aging varied SOC at 50 Celsius, and temperature at 50% state-of-charge.
           Cycle aging varied depth-of-discharge and average state-of-charge at 35 Celsius at
           charge and discharge rates of 1C.
           Relative discharge capacity is reported from measurements recorded at 35 Celsius and 1C rate.
           Relative DC resistance is reported after fitting of 10s 1C discharge pulses near 50% state-of-charge.

       MODEL SENSITIVITY
           The model predicts degradation rate versus time as a function of temperature and average
           state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as
           a function of average voltage and depth-of-discharge.

       MODEL LIMITATIONS
           Cycle degradation predictions are NOT SENSITIVE TO TEMPERATURE OR C-RATE. Cycling degradation predictions
           are ONLY ACCURATE NEAR 1C RATE AND 35 CELSIUS CELL TEMPERATURE.


   .. py:method:: calc_voltage(soc)


   .. py:method:: simulate_battery_life(input_timeseries, simulation_years = None, is_constant_input = False, breakpoints_max_time_diff_s = 86400, breakpoints_max_EFC_diff = 1)

      Run battery life simulation over the input, or repeat for the number of years specified.

      Updates attributes self.rates, self.stressors, self.outputs, and self.states inplace.

      :param input_timeseries:
      :type input_timeseries: dict, pd.DataFrame
      :param simulation_years:
      :type simulation_years: float
      :param is_constant_input:
      :type is_constant_input: bool
      :param breakpoints_max_time_diff_s:
      :type breakpoints_max_time_diff_s: float
      :param breakpoints_max_EFC_diff:
      :type breakpoints_max_EFC_diff: float



   .. py:method:: update_battery_state(t_secs, soc, T_celsius)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.

      :param t_secs: for the soc_timeseries data points
      :type t_secs: np.ndarray
      :param soc: Vector of the state-of-charge of the battery at each t_sec
      :type soc: np.ndarray
      :param T_celsius:
      :type T_celsius: ndarray



   .. py:method:: update_battery_state_repeating()

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
      degradation rates are unchanging for every timestep, and don't need to be calculated again.

      Updates self.states and self.outputs inplace.



   .. py:method:: update_outputs(stressors)

      Calculate outputs, based on current battery state (and maybe stressors)

      Updates self.outputs inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_rates(stressors)

      Calculate and update battery degradation rates based on stressor values

      Updates self.rates inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_states(stressors)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius

      Updates self.states inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:property:: cap


   .. py:attribute:: experimental_range


   .. py:attribute:: outputs


   .. py:attribute:: rates


   .. py:attribute:: states


   .. py:attribute:: stressors


.. py:class:: Nmc811_GrSi_LGMJ1_4Ah_Battery(degradation_scalar = 1, label = 'NMC811-GrSi LG MJ1')

   Bases: :py:obj:`blast.models.degradation_model.BatteryDegradationModel`


   Model fit to LG MJ1 cell aging data reported as part of the EU EVERLASTING battery project, report D2.3
   https://everlasting-project.eu/wp-content/uploads/2020/03/EVERLASTING_D2.3_final_20200228.pdf
   Cell tests were reported in early 2020, so likely 2018 or 2019 LG MJ1 cells.
   High energy density 18650s but poor cycle life.

   .. note::
       EXPERIMENTAL AGING DATA SUMMARY:
           Calendar aging varied SOC (10%, 70%, 90%) and temperature.
           Cycle aging varied temperature and C-rates; all DOD is 80% (10%-90%). NO ACCELERATED FADE OBSERVED.
           Relative discharge capacity is reported from measurements recorded at 25 Celsius and C/20 rate.

       MODEL SENSITIVITY
           The model predicts degradation rate versus time as a function of temperature and average
           state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as
           a function of C-rate, temperature, and depth-of-discharge (DOD dependence is assumed to be linear, no aging data)

       MODEL LIMITATIONS
           Cycle degradation predictions WILL NOT PREDICT KNEE-POINT due to limited data.
           OPERATION AT HIGH DOD PREDCTIONS ARE LIKELY INACCURATE (it is unclear what voltage window corresponds to SOCs defined in the test data).
           NMC811 is known to degrade quickly at voltages above 4.1 V.


   .. py:method:: simulate_battery_life(input_timeseries, simulation_years = None, is_constant_input = False, breakpoints_max_time_diff_s = 86400, breakpoints_max_EFC_diff = 1)

      Run battery life simulation over the input, or repeat for the number of years specified.

      Updates attributes self.rates, self.stressors, self.outputs, and self.states inplace.

      :param input_timeseries:
      :type input_timeseries: dict, pd.DataFrame
      :param simulation_years:
      :type simulation_years: float
      :param is_constant_input:
      :type is_constant_input: bool
      :param breakpoints_max_time_diff_s:
      :type breakpoints_max_time_diff_s: float
      :param breakpoints_max_EFC_diff:
      :type breakpoints_max_EFC_diff: float



   .. py:method:: update_battery_state(t_secs, soc, T_celsius)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.

      :param t_secs: for the soc_timeseries data points
      :type t_secs: np.ndarray
      :param soc: Vector of the state-of-charge of the battery at each t_sec
      :type soc: np.ndarray
      :param T_celsius:
      :type T_celsius: ndarray



   .. py:method:: update_battery_state_repeating()

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
      degradation rates are unchanging for every timestep, and don't need to be calculated again.

      Updates self.states and self.outputs inplace.



   .. py:method:: update_outputs(stressors)

      Calculate outputs, based on current battery state (and maybe stressors)

      Updates self.outputs inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_rates(stressors)

      Calculate and update battery degradation rates based on stressor values

      Updates self.rates inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_states(stressors)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius

      Updates self.states inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:property:: cap


   .. py:attribute:: experimental_range


   .. py:attribute:: outputs


   .. py:attribute:: rates


   .. py:attribute:: states


   .. py:attribute:: stressors


.. py:class:: Nmc_Lto_10Ah_Battery(degradation_scalar = 1, label = 'NMC-LTO')

   Bases: :py:obj:`blast.models.degradation_model.BatteryDegradationModel`


   Model fit to data reported by Bank et al from commercial NMC-LTO cells.
   https://doi.org/10.1016/j.jpowsour.2020.228566

   .. note::
       EXPERIMENTAL AGING DATA SUMMARY:
           Calendar aging varies temperature and SOC. There is almost no calendar aging impact
           at all until 80 Celsius.
           Cycle aging varies temperature, C-rate, and depth-of-discharge.

       MODEL SENSITIVITY
           The model predicts degradation rate versus time as a function of temperature and average
           state-of-charge and degradation rate versus equivalent full cycles (charge-throughput) as
           a function of C-rate, temperature, and depth-of-discharge (DOD dependence is assumed to be linear, no aging data)

       MODEL LIMITATIONS
           Calendar aging has competition between capacity gain and capacity loss. There is an experimental
           case (80 Celsius, 5% SOC) that has complex behavior not modeled here.
           Astonishingly enough, the cycling degradation model is actually _overestimating_ capacity fade for most cases.
           The exception here is at very high temperature (60+ Celsius), where the fade is high, but not quite as high as observed degradation.


   .. py:method:: simulate_battery_life(input_timeseries, simulation_years = None, is_constant_input = False, breakpoints_max_time_diff_s = 86400, breakpoints_max_EFC_diff = 1)

      Run battery life simulation over the input, or repeat for the number of years specified.

      Updates attributes self.rates, self.stressors, self.outputs, and self.states inplace.

      :param input_timeseries:
      :type input_timeseries: dict, pd.DataFrame
      :param simulation_years:
      :type simulation_years: float
      :param is_constant_input:
      :type is_constant_input: bool
      :param breakpoints_max_time_diff_s:
      :type breakpoints_max_time_diff_s: float
      :param breakpoints_max_EFC_diff:
      :type breakpoints_max_EFC_diff: float



   .. py:method:: update_battery_state(t_secs, soc, T_celsius)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is changing all the time.

      :param t_secs: for the soc_timeseries data points
      :type t_secs: np.ndarray
      :param soc: Vector of the state-of-charge of the battery at each t_sec
      :type soc: np.ndarray
      :param T_celsius:
      :type T_celsius: ndarray



   .. py:method:: update_battery_state_repeating()

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius. This function assumes battery load is repeating, i.e., stressors and
      degradation rates are unchanging for every timestep, and don't need to be calculated again.

      Updates self.states and self.outputs inplace.



   .. py:method:: update_outputs(stressors)

      Calculate outputs, based on current battery state (and maybe stressors)

      Updates self.outputs inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_rates(stressors)

      Calculate and update battery degradation rates based on stressor values

      Updates self.rates inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_states(stressors)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius

      Updates self.states inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:property:: cap


   .. py:attribute:: experimental_range


   .. py:attribute:: outputs


   .. py:attribute:: rates


   .. py:attribute:: states


   .. py:attribute:: stressors


