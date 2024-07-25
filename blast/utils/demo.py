"""Functions helpful for demonstrating the use of BLAST-Lite in Python."""

import numpy as np


def generate_sample_data() -> dict:
    """
    Generate synthetic sample data for demonstration purposes.

    Returns:
        Input       Dictionary with keys {'Time_s', 'SOC', 'Temperature_C'}
    """
    # Set up example time
    years = 10
    hours = years * 365 * 24
    t_hours = np.arange(hours + 1)

    # Set up an example profile
    soc = np.tile(
        np.array(
            [1, 1, 1, 1, 1, 1,
            0.4, 0.4, 0.6, 0.6, 1, 1,
            0.9, 0.8, 0.7, 0.6, 0.5, 0.4,
            0.3, 0.2, 0.2, 0.2, 0.2, 0.6]
            ),
        years * 365
    )
    soc = np.append(soc, 1)

    # Set up example temperature
    TdegC = (
        np.tile(np.concatenate([np.linspace(25, 35, 12), np.linspace(35, 25, 12)]),
        int(len(t_hours)/24))
        + np.tile(np.concatenate([np.linspace(-5, 5, 24*182), np.linspace(5, -5, 24*183)]),
        int(len(t_hours)/(24*365)))
    )
    TdegC = np.append(TdegC, 20)

    t_secs = t_hours * 3600

    input = {
        'Time_s': t_secs,
        'SOC': soc,
        'Temperature_C': TdegC
    }

    return input