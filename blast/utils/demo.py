"""Functions helpful for demonstrating the use of BLAST-Lite in Python."""

import numpy as np
import pandas as pd

from blast.utils.functions import assemble_one_year_input, get_nsrdb_temperature_data


def generate_sample_data(kind: str = "synthetic") -> dict:
    """
    Generate synthetic sample data for demonstration purposes.
    Options:

        1. Completely synthetic data
        2. Data from a small personal EV in Honolulu, Hawaii.
        3. Data from a large personal EV in Honolulu, Hawaii.
        4. Data from a commercial EV in Honolulu, Hawaii.

    Args:
        kind (str): One of 'synthetic', 'ev_smallbattery', 'ev_largebattery',
        'ev_commercial', 'ev_commercial_lowdod', 'ev_commercial_lowdod_lowsoc'

    Returns:
        dict:  Dictionary with keys {'Time_s', 'SOC', 'Temperature_C'}
    """

    _allowed_kinds = {
        "synthetic",
        "ev_smallbattery",
        "ev_largebattery",
        "ev_commercial",
        "ev_commercial_lowdod",
        "ev_commercial_lowdod_lowsoc",
    }

    if kind == "synthetic":
        # Set up example time
        years = 10
        hours = years * 365 * 24
        t_hours = np.arange(hours + 1)

        # Set up an example profile
        soc = np.tile(
            np.array(
                [
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    0.4,
                    0.4,
                    0.6,
                    0.6,
                    1,
                    1,
                    0.9,
                    0.8,
                    0.7,
                    0.6,
                    0.5,
                    0.4,
                    0.3,
                    0.2,
                    0.2,
                    0.2,
                    0.2,
                    0.6,
                ]
            ),
            years * 365,
        )
        soc = np.append(soc, 1)

        # Set up example temperature
        TdegC = np.tile(
            np.concatenate([np.linspace(25, 35, 12), np.linspace(35, 25, 12)]),
            int(len(t_hours) / 24),
        ) + np.tile(
            np.concatenate(
                [np.linspace(-5, 5, 24 * 182), np.linspace(5, -5, 24 * 183)]
            ),
            int(len(t_hours) / (24 * 365)),
        )
        TdegC = np.append(TdegC, 20)

        t_secs = t_hours * 3600

        input = {"Time_s": t_secs, "SOC": soc, "Temperature_C": TdegC}

    elif kind == "ev_smallbattery":
        climate = get_nsrdb_temperature_data("Honolulu, Hawaii")
        ev_smallbattery = pd.read_csv(
            "examples/application profiles/personal_ev_smallbatt.csv"
        )
        ev_smallbattery = ev_smallbattery.iloc[
            np.linspace(0, 24 * 3600 * 7 - 1, 24 * 7)
        ]
        input = assemble_one_year_input(ev_smallbattery, climate)

    elif kind == "ev_largebattery":
        climate = get_nsrdb_temperature_data("Honolulu, Hawaii")
        ev_largebattery = pd.read_csv(
            "examples/application profiles/personal_ev_largebatt.csv"
        )
        ev_largebattery = ev_largebattery.iloc[
            np.linspace(0, 24 * 3600 * 7 - 1, 24 * 7)
        ]
        input = assemble_one_year_input(ev_largebattery, climate)

    elif kind == "ev_commercial":
        climate = get_nsrdb_temperature_data("Honolulu, Hawaii")
        commercial_ev = pd.read_csv("examples/application profiles/commercial_ev.csv")
        commercial_ev = commercial_ev.iloc[np.linspace(0, 24 * 3600 * 7 - 1, 24 * 7)]
        input = assemble_one_year_input(commercial_ev, climate)

    elif kind == "ev_commercial_lowdod":
        climate = get_nsrdb_temperature_data("Honolulu, Hawaii")
        commercial_ev_lowdod = pd.read_csv('blast/application profiles/commercial_ev_lowdod.csv')
        commercial_ev_lowdod = commercial_ev_lowdod.iloc[np.linspace(0, 24*3600*7 - 1, 24*7)]
        input = assemble_one_year_input(commercial_ev_lowdod, climate)

    elif kind == "ev_commercial_lowdod_lowsoc":
        climate = get_nsrdb_temperature_data("Honolulu, Hawaii")
        commercial_ev_lowdod_lowsoc = pd.read_csv('blast/application profiles/commercial_ev_lowdod.csv')
        commercial_ev_lowdod_lowsoc = commercial_ev_lowdod_lowsoc.iloc[np.linspace(0, 24*3600*7 - 1, 24*7)]
        commercial_ev_lowdod_lowsoc['SOC'] += -0.4
        input = assemble_one_year_input(commercial_ev_lowdod_lowsoc, climate)

    else:
        raise ValueError(
            f"Expected `kind` to be one of {_allowed_kinds}. Received {kind}."
        )

    return input
