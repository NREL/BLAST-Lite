blast.models.lmo_gr_NissanLeaf66Ah_2ndLife_2020
===============================================

.. py:module:: blast.models.lmo_gr_NissanLeaf66Ah_2ndLife_2020


Classes
-------

.. autoapisummary::

   blast.models.lmo_gr_NissanLeaf66Ah_2ndLife_2020.Lmo_Gr_NissanLeaf66Ah_2ndLife_Battery


Module Contents
---------------

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


