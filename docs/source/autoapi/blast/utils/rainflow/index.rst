blast.utils.rainflow
====================

.. py:module:: blast.utils.rainflow

.. autoapi-nested-parse::

   Implements rainflow cycle counting algorithm for fatigue analysis according to section 5.4.4 in ASTM E1049-85 (2011).



Functions
---------

.. autoapisummary::

   blast.utils.rainflow.count_cycles
   blast.utils.rainflow.extract_cycles
   blast.utils.rainflow.reversals


Module Contents
---------------

.. py:function:: count_cycles(series, ndigits=None, nbins=None, binsize=None)

   Count cycles in the series.

   :param series: iterable sequence of numbers
   :param ndigits: Use a negative value to round to tens, hundreds, etc.
   :type ndigits: int
   :param nbins:
   :type nbins: int
   :param binsize:
   :type binsize: int
   :param Arguments ndigits:
   :param nbins and binsize are mutually exclusive.:

   :returns: A sorted list containing pairs of range and cycle count.
             The counts may not be whole numbers because the rainflow counting
             algorithm may produce half-cycles. If binning is used then ranges
             correspond to the right (high) edge of a bin.


.. py:function:: extract_cycles(series)

   Iterate cycles in the series.
   :param series:
   :type series: iterable sequence of numbers

   :Yields: **cycle** (*tuple*) -- Each tuple contains (range, mean, count, start index, end index).
            Count equals to 1.0 for full cycles and 0.5 for half cycles.


.. py:function:: reversals(series)

   Iterate reversal points in the series.
   A reversal point is a point in the series at which the first derivative
   changes sign. Reversal is undefined at the first (last) point because the
   derivative before (after) this point is undefined. The first and the last
   points are treated as reversals.
   :param series:
   :type series: iterable sequence of numbers

   :Yields: *Reversal points as tuples (index, value).*


