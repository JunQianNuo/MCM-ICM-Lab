# -*- coding: utf-8 -*-
"""Figure-local tests for plot_stair_wear.py.

Run from the figure directory:

    python -m pytest -q test_plot_stair_wear.py

The closed-form demo is deliberately separable, so every projection below has
an exact hand solution.
"""

import numpy as np
import pytest

import plot_stair_wear as f


def _demo_weights():
    i = np.arange(f.DEMO_NX, dtype=float)
    g = np.where(i < f.DEMO_NX // 2, i + 1.0, f.DEMO_NX - i)
    g = g / g.max()
    h = (f.DEMO_NY - np.arange(f.DEMO_NY, dtype=float)) / float(f.DEMO_NY)
    return g, h


def test_demo_field_is_closed_form():
    depth = f.demo_depth_mm()
    assert depth.shape == (f.DEMO_NY, f.DEMO_NX)
    assert depth.max() == pytest.approx(f.DEMO_PEAK_MM, rel=0.0, abs=1e-12)
    # two centre columns at the front edge carry the peak
    assert depth[0, 11] == pytest.approx(3.0, rel=0.0, abs=1e-12)
    assert depth[0, 12] == pytest.approx(3.0, rel=0.0, abs=1e-12)
    # far corner: 3 * (1/16) * (1/12) mm
    assert depth[-1, 0] == pytest.approx(3.0 / 192.0, rel=0.0, abs=1e-15)


def test_demo_projections_match_hand_calculation():
    depth = f.demo_depth_mm()
    g, h = _demo_weights()
    x_proj, y_proj = f.compute_projections(depth, 10.0, 10.0)

    np.testing.assert_allclose(x_proj, f.DEMO_X_PEAK_MM3 * g, rtol=1e-12)
    np.testing.assert_allclose(y_proj, f.DEMO_Y_FRONT_MM3 * h, rtol=1e-12)

    assert x_proj[0] == pytest.approx(f.DEMO_X_EDGE_MM3, rel=1e-12)
    assert x_proj[11] == pytest.approx(f.DEMO_X_PEAK_MM3, rel=1e-12)
    assert y_proj[0] == pytest.approx(f.DEMO_Y_FRONT_MM3, rel=1e-12)
    assert y_proj[-1] == pytest.approx(f.DEMO_Y_BACK_MM3, rel=1e-12)

    # both marginal sums recover the same total volume
    assert x_proj.sum() == pytest.approx(f.DEMO_TOTAL_MM3, rel=1e-12)
    assert y_proj.sum() == pytest.approx(f.DEMO_TOTAL_MM3, rel=1e-12)


def test_uniform_field_hand_solution():
    depth = np.full((3, 4), 2.0)          # 3 rows (y), 4 columns (x)
    x_proj, y_proj = f.compute_projections(depth, 3.0, 5.0)
    np.testing.assert_allclose(x_proj, np.full(4, 90.0), rtol=1e-12)
    np.testing.assert_allclose(y_proj, np.full(3, 120.0), rtol=1e-12)
    assert x_proj.sum() == pytest.approx(360.0, rel=1e-12)
    assert y_proj.sum() == pytest.approx(360.0, rel=1e-12)


def test_projection_sum_invariant():
    rng = np.random.default_rng(20250109)
    depth = rng.random((9, 13)) * 2.5
    x_proj, y_proj = f.compute_projections(depth, 2.0, 4.0)
    expected = depth.sum() * 8.0
    assert x_proj.sum() == pytest.approx(expected, rel=1e-12)
    assert y_proj.sum() == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("bad", [
    np.ones((1, 5)),
    np.ones((5, 1)),
    np.ones((3,)),                       # 1-D
    np.array([[1.0, np.nan], [1.0, 2.0]]),
    np.array([[1.0, np.inf], [1.0, 2.0]]),
    np.array([[1.0, -0.5], [1.0, 2.0]]),
])
def test_validate_depth_rejects_bad_input(bad):
    with pytest.raises(ValueError):
        f.validate_depth(bad)


def test_validate_depth_clips_tiny_negatives():
    depth = f.validate_depth(np.array([[-1e-12, 1.0], [2.0, 3.0]]))
    assert depth.min() >= 0.0


def test_load_depth_csv_roundtrip(tmp_path):
    arr = np.arange(12, dtype=float).reshape(3, 4)
    path = tmp_path / "depth.csv"
    np.savetxt(path, arr, delimiter=",")
    np.testing.assert_allclose(f.load_depth(path), arr)


def test_load_depth_accepts_comma_header(tmp_path):
    path = tmp_path / "depth_header.csv"
    path.write_text("x0,x1,x2\n1,2,3\n4,5,6\n", encoding="utf-8")
    np.testing.assert_allclose(f.load_depth(path), [[1, 2, 3], [4, 5, 6]])


def test_load_depth_rejects_malformed_first_data_row(tmp_path):
    path = tmp_path / "ragged.csv"
    path.write_text("1,2\n3,4,5\n6,7,8\n", encoding="utf-8")
    with pytest.raises(ValueError):
        f.load_depth(path)


def test_projection_cross_check(tmp_path):
    path = tmp_path / "x_wear_volume_mm3.csv"
    expected = np.array([0.0, 212.5, 2550.0, 212.5])
    np.savetxt(path, np.column_stack([np.arange(4.0), expected]), delimiter=",",
               header="x_mm,x_wear_volume_mm3", comments="")
    result = f.compare_projection(expected.copy(), path, "x", 1e-9)
    assert result["ok"] is True
    assert result["max_rel_error"] < 1e-9


def test_self_test_passes(capsys):
    assert f.run_self_test() == 0
    assert "ALL CHECKS PASSED" in capsys.readouterr().out


def test_render_smoke(tmp_path):
    pytest.importorskip("matplotlib")
    depth = f.demo_depth_mm()
    x_proj, y_proj = f.compute_projections(depth, 10.0, 10.0)
    outputs = f.plot_figure(
        depth, x_proj, y_proj, dx_mm=10.0, dy_mm=10.0, x0_mm=0.0, y0_mm=0.0,
        outdir=tmp_path, stem="smoke", source_label="smoke test",
        synthetic=True, editable_svg=False, input_desc="demo",
    )
    assert len(outputs) == 3
    for path in outputs:
        assert path.exists() and path.stat().st_size > 1000
