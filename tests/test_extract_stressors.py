"""Unit tests for BatteryDegradationModel stressor extraction.

These tests verify that the stressor quantities derived from SOC and
temperature timeseries (average SOC, depth-of-discharge, equivalent full
cycles, elapsed time, temperature, and C-rate) are calculated accurately.
"""

import numpy as np
import pandas as pd
import pytest

from blast.models.degradation_model import BatteryDegradationModel


@pytest.fixture
def model():
    """A bare BatteryDegradationModel with a healthy (q=1) state.

    ``_extract_stressors`` depends on ``self.outputs["q"]`` to rescale
    charge-throughput-based stressors by the current state-of-health, so the
    base model (which has empty ``outputs``) must be given a value.
    """
    m = BatteryDegradationModel()
    m.outputs = {"q": np.array([1.0])}
    return m


# ---------------------------------------------------------------------------
# Equivalent full cycles (EFC)
# ---------------------------------------------------------------------------
class TestEquivalentFullCycles:
    def test_single_full_cycle(self, model):
        """One full 0->1->0 SOC swing is exactly one equivalent full cycle."""
        t_secs = np.array([0, 3600, 7200], dtype=float)
        soc = np.array([0.0, 1.0, 0.0])
        T_celsius = np.full_like(soc, 25.0)

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        assert stressors["delta_efc"] == pytest.approx(1.0)

    def test_half_cycle(self, model):
        """A single 0->1 charge (no discharge back) is half an EFC."""
        t_secs = np.array([0, 3600], dtype=float)
        soc = np.array([0.0, 1.0])
        T_celsius = np.full_like(soc, 25.0)

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        assert stressors["delta_efc"] == pytest.approx(0.5)

    def test_multiple_cycles(self, model):
        """Three full 0<->1 swings accumulate to three EFCs."""
        soc = np.array([0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0])
        t_secs = np.arange(len(soc), dtype=float) * 3600
        T_celsius = np.full_like(soc, 25.0)

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        assert stressors["delta_efc"] == pytest.approx(3.0)

    def test_partial_swings(self, model):
        """EFC equals the total SOC path length divided by two."""
        soc = np.array([0.2, 0.8, 0.4, 0.6])
        t_secs = np.arange(len(soc), dtype=float) * 3600
        T_celsius = np.full_like(soc, 25.0)

        expected = (0.6 + 0.4 + 0.2) / 2  # sum(|dSOC|)/2
        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        assert stressors["delta_efc"] == pytest.approx(expected)

    def test_efc_scales_with_state_of_health(self, model):
        """EFC is rescaled by the current state-of-health (q)."""
        model.outputs["q"] = np.array([0.8])
        soc = np.array([0.0, 1.0, 0.0])
        t_secs = np.arange(len(soc), dtype=float) * 3600
        T_celsius = np.full_like(soc, 25.0)

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        assert stressors["delta_efc"] == pytest.approx(1.0 * 0.8)


# ---------------------------------------------------------------------------
# Depth-of-discharge (DOD)
# ---------------------------------------------------------------------------
class TestDepthOfDischarge:
    def test_full_dod(self, model):
        soc = np.array([0.0, 1.0, 0.0])
        t_secs = np.arange(len(soc), dtype=float) * 3600
        T_celsius = np.full_like(soc, 25.0)

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        assert stressors["dod"] == pytest.approx(1.0)

    def test_partial_dod(self, model):
        soc = np.array([0.3, 0.7, 0.3])
        t_secs = np.arange(len(soc), dtype=float) * 3600
        T_celsius = np.full_like(soc, 25.0)

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        assert stressors["dod"] == pytest.approx(0.4)

    def test_zero_dod_storage(self, model):
        soc = np.array([0.5, 0.5, 0.5])
        t_secs = np.arange(len(soc), dtype=float) * 3600
        T_celsius = np.full_like(soc, 25.0)

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        assert stressors["dod"] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Average SOC
# ---------------------------------------------------------------------------
class TestAverageSOC:
    def test_soc_array_preserved(self, model):
        """The soc stressor holds the raw profile; its mean is the average SOC."""
        soc = np.array([0.2, 0.4, 0.6, 0.8])
        t_secs = np.arange(len(soc), dtype=float) * 3600
        T_celsius = np.full_like(soc, 25.0)

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        np.testing.assert_allclose(stressors["soc"], soc)
        assert np.mean(stressors["soc"]) == pytest.approx(0.5)

    def test_average_soc_recorded_by_update(self, model):
        """update_battery_state stores the mean SOC in the stressor history."""
        soc = np.array([0.2, 0.4, 0.6, 0.8])
        t_secs = np.arange(len(soc), dtype=float) * 3600
        T_celsius = np.full_like(soc, 25.0)

        model.update_battery_state(t_secs, soc, T_celsius)
        assert model.stressors["soc"][-1] == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# Elapsed time
# ---------------------------------------------------------------------------
class TestElapsedTime:
    def test_delta_t_days_one_day(self, model):
        t_secs = np.array([0.0, 86400.0])
        soc = np.array([0.5, 0.5])
        T_celsius = np.full_like(soc, 25.0)

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        assert stressors["delta_t_days"] == pytest.approx(1.0)

    def test_delta_t_days_fractional(self, model):
        t_secs = np.array([0.0, 43200.0])  # half a day
        soc = np.array([0.5, 0.5])
        T_celsius = np.full_like(soc, 25.0)

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        assert stressors["delta_t_days"] == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# Temperature
# ---------------------------------------------------------------------------
class TestTemperature:
    def test_celsius_to_kelvin(self, model):
        t_secs = np.array([0.0, 3600.0])
        soc = np.array([0.5, 0.5])
        T_celsius = np.array([25.0, 35.0])

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        np.testing.assert_allclose(stressors["TdegK"], np.array([298.15, 308.15]))


# ---------------------------------------------------------------------------
# C-rate
# ---------------------------------------------------------------------------
class TestCrate:
    def test_one_c(self, model):
        """A full 0->1 SOC change over one hour is a 1C rate."""
        t_secs = np.array([0.0, 3600.0])
        soc = np.array([0.0, 1.0])
        T_celsius = np.full_like(soc, 25.0)

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        assert stressors["Crate"] == pytest.approx(1.0)

    def test_storage_crate_is_zero(self, model):
        """A constant SOC (storage) yields a zero, non-nan C-rate."""
        t_secs = np.array([0.0, 3600.0, 7200.0])
        soc = np.array([0.5, 0.5, 0.5])
        T_celsius = np.full_like(soc, 25.0)

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        assert stressors["Crate"] == pytest.approx(0.0)
        assert not np.isnan(stressors["Crate"])


# ---------------------------------------------------------------------------
# Input type casting on core functions
# ---------------------------------------------------------------------------
class TestInputTypeCasting:
    def test_extract_stressors_accepts_lists(self, model):
        """Non-array (list) inputs are cast to float arrays internally."""
        t_secs = [0, 3600, 7200]
        soc = [0, 1, 0]
        T_celsius = [25, 25, 25]

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        assert stressors["delta_efc"] == pytest.approx(1.0)
        assert stressors["dod"] == pytest.approx(1.0)

    def test_extract_stressors_accepts_integer_arrays(self, model):
        """Integer SOC arrays don't cause integer-division artifacts."""
        t_secs = np.array([0, 3600, 7200])
        soc = np.array([0, 1, 0])
        T_celsius = np.array([25, 25, 25])

        stressors = model._extract_stressors(t_secs, soc, T_celsius)
        assert stressors["soc"].dtype == float
        assert stressors["dod"] == pytest.approx(1.0)

    def test_simulate_battery_life_accepts_lists_in_dict(self):
        """simulate_battery_life casts list inputs and validates keys."""
        model = BatteryDegradationModel()
        model.outputs = {"q": np.array([1.0])}
        input_timeseries = {
            "Time_s": [0, 3600, 7200],
            "SOC": [0.0, 1.0, 0.0],
            "Temperature_C": [25, 25, 25],
        }
        # Should run without raising a type error from list inputs.
        model.simulate_battery_life(input_timeseries)

    def test_simulate_battery_life_accepts_dataframe(self):
        model = BatteryDegradationModel()
        model.outputs = {"q": np.array([1.0])}
        df = pd.DataFrame(
            {
                "Time_s": [0, 3600, 7200],
                "SOC": [0.0, 1.0, 0.0],
                "Temperature_C": [25, 25, 25],
            }
        )
        model.simulate_battery_life(df)

    def test_missing_key_raises(self):
        model = BatteryDegradationModel()
        with pytest.raises(ValueError):
            model.simulate_battery_life({"Time_s": [0, 1], "SOC": [0, 1]})

    def test_bad_input_type_raises(self):
        model = BatteryDegradationModel()
        model.outputs = {"q": np.array([1.0])}
        with pytest.raises(TypeError):
            # A list is neither a dict nor a DataFrame; but the key check runs
            # first, so use a non-mapping object that still has required keys.
            model.simulate_battery_life(
                pd.Series(
                    {
                        "Time_s": [0, 1],
                        "SOC": [0, 1],
                        "Temperature_C": [25, 25],
                    }
                )
            )

    def test_update_battery_state_rejects_non_arrays(self, model):
        with pytest.raises(TypeError):
            model.update_battery_state([0, 1], np.array([0, 1]), np.array([25, 25]))

    def test_update_battery_state_length_mismatch(self, model):
        with pytest.raises(ValueError):
            model.update_battery_state(
                np.array([0, 1, 2]), np.array([0, 1]), np.array([25, 25])
            )
