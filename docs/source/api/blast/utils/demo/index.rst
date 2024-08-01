blast.utils.demo
================

.. py:module:: blast.utils.demo

.. autoapi-nested-parse::

   Functions helpful for demonstrating the use of BLAST-Lite in Python.



Functions
---------

.. autoapisummary::

   blast.utils.demo.generate_sample_data


Module Contents
---------------

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


