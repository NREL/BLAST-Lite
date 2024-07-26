==========
User Guide
==========


.. code:: python

    from blast import utils, models


Generate input data
--------------------------------

.. code:: python

    # Temperature data
    climate = utils.get_nsrdb_temperature_data('Honolulu, Hawaii')

    # SOC profile - use example or load your own
    ev_smallbattery = utils.generate_sample_profile('ev_smallbattery')

    # Combine the two
    input_ev_smallbatt = utils.assemble_one_year_input(ev_smallbattery, climate)



