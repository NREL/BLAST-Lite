blast.utils.functions
=====================

.. py:module:: blast.utils.functions

.. autoapi-nested-parse::

   Helper functions for fetching and formatting simulation input data.



Functions
---------

.. autoapisummary::

   blast.utils.functions.assemble_one_year_input
   blast.utils.functions.get_nsrdb_temperature_data
   blast.utils.functions.make_inputs_periodic
   blast.utils.functions.tile_to_one_year


Module Contents
---------------

.. py:function:: assemble_one_year_input(soc, climate, time_interval = 3600)

   Construct SOC and climate data into a single dataframe,
   linked by timestamp.

   :param soc: Dataframe with SOC profile and 'Time_s' timestamp
   :type soc: pd.DataFrame
   :param climate: Dataframe with temperature data and 'Time_s' timestamp column
   :type climate: pd.DataFrame
   :param time_interval: Interval at which the final data is resampled to
   :type time_interval: int

   :returns: Combined DataFrame with time, soc, and temperature columns.
   :rtype: pd.DataFrame


.. py:function:: get_nsrdb_temperature_data(location = 'Honolulu, Hawaii')

   Get temperature time-series data from NSRDB at the nearest coordinates to the
   provided location, and format for battery life simulation.

   :param location: Descriptive location string.
   :type location: str

   :returns: Temperature time series at the specified location.
   :rtype: pd.DataFrame


.. py:function:: make_inputs_periodic(input_timeseries, interp_time_window_hours)

   Linearly interpolate the last 'interp_time_window_hours' of the input to ensure periodicity.

   :param input_timeseries: Dictionary with keys ['Temperature_C', 'SOC', 'Time_s']
                            and timeseries values.
   :type input_timeseries: dict
   :param interp_time_window_hours: Number of hours to interpolate, from end of timeseries
   :type interp_time_window_hours: float

   :returns: Dictionary with keys ['Temperature_C', 'SOC', 'Time_s'] after interpolation.
   :rtype: dict


.. py:function:: tile_to_one_year(df)

   Take a dataframe with timestamp column called 'Time_s'
   and resample to 1 year.

   :param df: Dataframe to cut off to 1 year.
   :type df: pd.DataFrame

   :returns: Dataframe after cutting off to 1 year.
   :rtype: pd.DataFrame


