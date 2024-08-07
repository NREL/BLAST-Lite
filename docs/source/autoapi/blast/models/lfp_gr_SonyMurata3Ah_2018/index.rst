blast.models.lfp_gr_SonyMurata3Ah_2018
======================================

.. py:module:: blast.models.lfp_gr_SonyMurata3Ah_2018


Classes
-------

.. autoapisummary::

   blast.models.lfp_gr_SonyMurata3Ah_2018.Lfp_Gr_SonyMurata3Ah_Battery


Module Contents
---------------

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


