blast.models.degradation_model
==============================

.. py:module:: blast.models.degradation_model


Classes
-------

.. autoapisummary::

   blast.models.degradation_model.BatteryDegradationModel


Module Contents
---------------

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


