blast.models.lfp_gr_250AhPrismatic_2019
=======================================

.. py:module:: blast.models.lfp_gr_250AhPrismatic_2019


Classes
-------

.. autoapisummary::

   blast.models.lfp_gr_250AhPrismatic_2019.Lfp_Gr_250AhPrismatic


Module Contents
---------------

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


