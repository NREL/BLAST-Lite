blast.models.nmc111_gr_Sanyo2Ah_2014
====================================

.. py:module:: blast.models.nmc111_gr_Sanyo2Ah_2014


Classes
-------

.. autoapisummary::

   blast.models.nmc111_gr_Sanyo2Ah_2014.Nmc111_Gr_Sanyo2Ah_Battery


Module Contents
---------------

.. py:class:: Nmc111_Gr_Sanyo2Ah_Battery(degradation_scalar = 1, label = 'NMC111-Gr Sanyo')



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


   .. py:method:: extract_stressors(t_secs, soc, T_celsius)
      :staticmethod:


      Extract stressor values including time, temperature, depth-of-discharge,
      temperature, Ua, C-rate, and cycles.

      :param t_secs:
      :type t_secs: np.ndarray
      :param soc:
      :type soc: np.ndarray
      :param T_celsius:
      :type T_celsius: np.ndarray

      :returns: stressors (dict)        Dictionary of stressor values.



   .. py:method:: find_breakpoints(time_seconds, EFCs, idx_turning_point, max_time_diff_s=86400, max_EFC_diff=1)
      :staticmethod:


      Find breakpoints for simulating battery degradation, where degradation is calculated once
      a certain number of equivalent full cycles has passed, otherwise after a certain time has passed.
      Defaults are to calculate degradation every 1 EFC or every day.

      :param time_seconds:
      :type time_seconds: np.ndarray
      :param EFCs:
      :type EFCs: np.ndarray
      :param idx_turning_point: the end of cycles
      :type idx_turning_point: np.ndarray
      :param max_time_diff_s:
      :type max_time_diff_s: int
      :param max_EFC_diff:
      :type max_EFC_diff: int

      :returns: breakpoints (list)



   .. py:method:: get_Ua(soc)
      :staticmethod:


      Calculate Ua from SOC via lithiation fraction.

      :param soc:
      :type soc: np.ndarray

      :returns: Ua



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



   .. py:method:: update_power_B_state(y0, dx, k, p)
      :staticmethod:


      Update time-varying power B state

      :param TODO Paul:

      :returns: TODO Paul



   .. py:method:: update_power_state(y0, dx, k, p)
      :staticmethod:


      Update time-varying power state

      :param TODO Paul:

      :returns: TODO Paul



   .. py:method:: update_rates(stressors)

      Calculate and update battery degradation rates based on stressor values

      Updates self.rates inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:method:: update_sigmoid_state(y0, dx, y_inf, k, p)
      :staticmethod:


      Update time-varying sigmoid state

      :param TODO Paul:

      :returns: TODO Paul



   .. py:method:: update_states(stressors)

      Update the battery states, based both on the degradation state as well as the battery performance
      at the ambient temperature, T_celsius

      Updates self.states inplace.

      :param stressors:
      :type stressors: dict)    Output from extract_stressors(



   .. py:attribute:: experimental_range


   .. py:attribute:: outputs


   .. py:attribute:: rates


   .. py:attribute:: states


   .. py:attribute:: stressors


