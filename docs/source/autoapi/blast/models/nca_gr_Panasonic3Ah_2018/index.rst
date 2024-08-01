blast.models.nca_gr_Panasonic3Ah_2018
=====================================

.. py:module:: blast.models.nca_gr_Panasonic3Ah_2018


Classes
-------

.. autoapisummary::

   blast.models.nca_gr_Panasonic3Ah_2018.Nca_Gr_Panasonic3Ah_Battery


Module Contents
---------------

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


