==========
User Guide
==========

This is a detailed user guide that demonstrates how to use the BLAST-Lite battery life models,
using the Kokam NMC111|Gr 75Ah battery life model as an example. The battery modeled here is a
high-power cell with long cycle life. Because nominal cell resistance is low, the relative change
of resistance at end-of-life is quite high compared to other cell designs (~300% increase in cell
resistance at 80% capacity if not more). Fade rates can be changed in the code to accomodate other cell models.
Documentation is provided in the life model class. See https://ieeexplore.ieee.org/abstract/document/7963578
for the aging test details and results used to parameterize this model.


.. code:: python
    
    from blast import utils, models


Generating input data
--------------------------------

To run, the life model needs three time series:

#. The time in seconds since beginning-of-life of the battery
#. The state-of-charge profile of the battery (0 to 1)
#. The ambient temperature (or battery temperature, if you have a thermal model).

These three time series should be fed to the life model in the form
of a dictionary with these keys:

* ``'Time_s'``
* ``'SOC'``
* ``'Temperature_C'``

and each time series should have the same length.


You can load your own data or fetch temperature data from the NSRDB
database using the the function ``utils.get_nsrdb_temperature_data()``.

.. code:: python

    climate = utils.get_nsrdb_temperature_data('Honolulu, Hawaii')

For the purpose of demonstration, we have a few sample input datasets:

.. code:: python

    # Completely synthetic data
    input_synthetic = utils.generate_sample_data('synthetic')

    # Small personal EV in Honolulu, Hawaii
    input_ev_smallbattery = utils.generate_sample_data('ev_smallbattery')

    # Large personal EV in Honolulu, Hawaii
    input_ev_largebattery = utils.generate_sample_data('ev_largebattery')

    # Commercial EV in Honolulu, Hawaii
    input_ev_commercial = utils.generate_sample_data('ev_commercial')

Visualize the sample profiles:

.. code:: python

    import matplotlib.pyplot as plt

    plt.plot(input_ev_largebattery['Time_s'] / (24*3600), input_ev_largebattery['SOC'], label='Large battery EV')
    plt.plot(input_ev_smallbattery['Time_s'] / (24*3600), input_ev_smallbattery['SOC'], label='Small battery EV')
    plt.plot(input_ev_commercial['Time_s'] / (24*3600), input_ev_commercial['SOC'], label='Commercial EV')
    plt.xlabel('Time (days)')
    plt.ylabel('State of charge')
    plt.legend()

.. image:: assets/ev_profiles.PNG
    :width: 600


Instantiating a model
--------------------------------

To see a list of available models, run:

.. code-block:: python

    >> models.available_models()
    
    ['Lfp_Gr_250AhPrismatic', 'Lfp_Gr_SonyMurata3Ah_Battery', 'Lmo_Gr_NissanLeaf66Ah_2ndLife_Battery', 'NCA_GrSi_SonyMurata2p5Ah_Battery', 'NMC_Gr_50Ah_B1', 'NMC_Gr_50Ah_B2', 'NMC_Gr_75Ah_A', 'Nca_Gr_Panasonic3Ah_Battery', 'Nmc111_Gr_Kokam75Ah_Battery', 'Nmc111_Gr_Sanyo2Ah_Battery', 'Nmc811_GrSi_LGMJ1_4Ah_Battery', 'Nmc_Lto_10Ah_Battery']

Select a model and instantiate a cell:

.. code:: python

    from blast.models import Nmc111_Gr_Kokam75Ah_Battery
    
    >> cell = Nmc111_Gr_Kokam75Ah_Battery()

All battery models have five builtin properties stored as attributes of the model class:

.. code:: python

    >> cell._cap
    >> cell.states
    >> cell.outputs
    >> cell.stressors
    >> cell.rates

The first is **_cap**, which is the nominal discharge capacity of the cell in Amp hours.
The next four track battery lifetime values, and store the history of the battery as lifetime is simulated  at each timestep/iteration:

- `states`: internal states of the battery model
    - Ex., time-dependent capacity loss
- `outputs`: battery properties calculated from states
    - Ex., relative discharge capacity
- `stressors`: values of stressors used by the model
    - Ex., temperature, depth-of-discharge, charge-throughput
    - Note that degradation rates are calculate from stressor timeseries, and then normalized for the timestep; for example, an Arrhenius expression would be evaluated from the temperature timeseries for the entire timestep, and then normalized by taking the time-based average - this gives a different value than if the Arrhenius expression was evaluated on the averager temperature. Other normalizations can include using the minimum or maximum value over the timestep, or using the root-mean-square.
- `rates`: values of degradation rates
    - Ex., time-dependent degradation rate due to temperature and state-of-charge

Battery models may have other properties, such as the open-circuit voltage as a function of state-of-charge, nominal DC resistance values, or first-life/second-life capacity definitions.

This specific battery model is relatively complex, and has many states and outputs that describe the degradation state of the battery. Properties `states`, `outputs`, `stressors`, and `rates` are all stored as dicts.


Running the simulation
--------------------------------


Evaluating results
--------------------------------
