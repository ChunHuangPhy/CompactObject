import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from TOVsolver.unit import Msun, dyn_cm_2, g_cm_3, km


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(params=["Test_Case", "docs/source"])
def notebook(request):
    path = ROOT / request.param / "test_Inference_polytrope.ipynb"
    cells = json.loads(path.read_text())["cells"]
    namespace = {}
    for index in [1, 3, 5, 7]:
        exec(compile("".join(cells[index]["source"]), str(path), "exec"), namespace)
    yield cells, namespace
    plt.close("all")


def test_polytrope_pressure_uses_all_three_segments_in_consistent_units(notebook):
    _, ns = notebook
    rho1, rho2 = np.array([1.5, 4.0]) * ns["rho_ns"]
    theta = np.array([rho1, rho2, 2.0, 3.0, 2.0, 0.0])
    energy, pressure = ns["eos_from_parameters"](theta)
    rho = ns["core_density"] / g_cm_3
    rho0 = ns["rho_crust_end"] / g_cm_3
    p0 = ns["P_crust_end"] / dyn_cm_2
    p1 = p0 * (rho1 / rho0) ** 2
    p2 = p1 * (rho2 / rho1) ** 3
    expected = np.where(
        rho <= rho1,
        p0 * (rho / rho0) ** 2,
        np.where(rho <= rho2, p1 * (rho / rho1) ** 3, p2 * (rho / rho2) ** 2),
    )
    np.testing.assert_allclose(pressure[-len(rho):] / dyn_cm_2, expected, rtol=1e-12)
    assert np.all(np.diff(energy) > 0)
    assert np.all(np.diff(pressure) > 0)


def test_prior_upper_density_uses_the_grid_actually_solved(notebook, monkeypatch):
    _, ns = notebook
    grid = np.logspace(14.3, 15.6, 6)
    calls = []

    def truncated_sequence(filename, energy, pressure, central_density_range):
        calls.append(central_density_range)
        return np.array([[1.0 * Msun, 12 * km], [2.0 * Msun, 11 * km], [1.9 * Msun, 10 * km]])

    monkeypatch.setattr(ns["main"], "OutputMR", truncated_sequence)
    theta = ns["prior_transform"](np.full(6, 0.5), central_density_grid_cgs=grid)
    np.testing.assert_allclose(calls[0] / g_cm_3, grid)
    assert theta[5] == pytest.approx((np.log10(grid[0]) + np.log10(grid[1])) / 2)


def test_default_unit_check_does_not_integrate_stars(notebook, monkeypatch, capsys):
    cells, ns = notebook

    def unexpected_tov(*args, **kwargs):
        pytest.fail("The default unit check must not call the TOV solver")

    monkeypatch.setattr(ns["main"], "OutputMR", unexpected_tov)
    monkeypatch.setattr(plt, "show", lambda: None)
    exec("".join(cells[9]["source"]), ns)
    assert len(ns["diagnostic_parameters"]) == 3
    assert "Optional M-R survey skipped" in capsys.readouterr().out


def test_empty_stellar_sequence_is_rejected_without_indexing_it(notebook, monkeypatch):
    _, ns = notebook
    monkeypatch.setattr(ns["main"], "OutputMR", lambda *args, **kwargs: np.empty((0, 2)))
    theta = ns["prior_transform"](np.full(6, 0.5))
    assert theta[5] == 0
    assert ns["likelihood_transform"](theta) == ns["REJECTED_LOG_LIKELIHOOD"]
    theta[5] = 14.8
    assert ns["likelihood_transform"](theta) == ns["REJECTED_LOG_LIKELIHOOD"]
