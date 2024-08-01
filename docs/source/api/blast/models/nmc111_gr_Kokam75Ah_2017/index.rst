blast.models.nmc111_gr_Kokam75Ah_2017
=====================================

.. py:module:: blast.models.nmc111_gr_Kokam75Ah_2017


Classes
-------

.. autoapisummary::

   blast.models.nmc111_gr_Kokam75Ah_2017.Nmc111_Gr_Kokam75Ah_Battery


Module Contents
---------------

.. py:class:: Nmc111_Gr_Kokam75Ah_Battery(degradation_scalar = 1, label = 'NMC111-Gr Kokam')



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


