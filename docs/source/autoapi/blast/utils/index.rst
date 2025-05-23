blast.utils
===========

.. py:module:: blast.utils


Submodules
----------

.. toctree::
   :maxdepth: 1

   /autoapi/blast/utils/demo/index
   /autoapi/blast/utils/functions/index
   /autoapi/blast/utils/rainflow/index


Functions
---------

.. autoapisummary::

   blast.utils.assemble_one_year_input
   blast.utils.generate_sample_data
   blast.utils.get_nsrdb_temperature_data


Package Contents
----------------

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


.. py:function:: generate_sample_data(kind = 'synthetic')

   Generate synthetic sample data for demonstration purposes.
   Options:

       1. Completely synthetic data
       2. Data from a small personal EV in Honolulu, Hawaii.
       3. Data from a large personal EV in Honolulu, Hawaii.
       4. Data from a commercial EV in Honolulu, Hawaii.

   :param kind: One of 'synthetic', 'ev_smallbattery', 'ev_largebattery',
   :type kind: str
   :param 'ev_commercial':
   :param 'ev_commercial_lowdod':
   :param 'ev_commercial_lowdod_lowsoc':

   :returns: Dictionary with keys {'Time_s', 'SOC', 'Temperature_C'}
   :rtype: dict


.. py:function:: get_nsrdb_temperature_data(location = 'Honolulu, Hawaii')

   Get temperature time-series data from NSRDB at the nearest coordinates to the
   provided location, and format for battery life simulation.

   :param location: Descriptive location string.
   :type location: str

   :returns: Temperature time series at the specified location.
   :rtype: pd.DataFrame


