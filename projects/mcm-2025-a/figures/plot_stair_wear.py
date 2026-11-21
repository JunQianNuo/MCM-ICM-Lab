#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stair tread wear figure: depth heatmap + lateral/longitudinal projections.

Data contract (matches the WDM implementation in MCM-ICM-Lab PR #2):

    depth_mm[y][x]   wear depth in mm, positive = material removed
    y = 0            front edge of the tread (nosing); y increases inward
    x = 0            left edge of the tread;       x increases to the right

    dx_mm, dy_mm     cell size in mm; cell area = dx_mm * dy_mm

Derived projections (the names are deliberately the program's output names):

    x_wear_volume_mm3[x] = sum_y depth_mm[y, x] * dx_mm * dy_mm
    y_wear_volume_mm3[y] = sum_x depth_mm[y, x] * dx_mm * dy_mm

Both are in mm^3 per grid column/row.  The figure shows the field and the two
marginal projections aligned with the axes they are taken along.

The default `--demo` grid is a *closed-form synthetic* field (clearly labelled
SYNTHETIC in every output).  It is not measured data and must not be presented
as the paper's result.

Usage
-----
    # synthetic first look + self check
    python plot_stair_wear.py --demo --self-test --outdir out/2501909

    # real WDM matrix (rows = y, columns = x) + cross-check against program output
    python plot_stair_wear.py \
        --input depth_mm.csv --dx-mm 10 --dy-mm 10 \
        --check-x-proj x_wear_volume_mm3.csv \
        --check-y-proj y_wear_volume_mm3.csv \
        --source-label "WDM depth_mm[y][x]" --outdir out/2501909
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

DEMO_NX = 24
DEMO_NY = 16
DEMO_PEAK_MM = 3.0
DEMO_DX_MM = 10.0
DEMO_DY_MM = 10.0

PROJ_X_NAME = "x_wear_volume_mm3"
PROJ_Y_NAME = "y_wear_volume_mm3"

# Closed-form demo expectations, also used by --self-test and the unit tests.
DEMO_X_EDGE_MM3 = 212.5      # = 3.0 * (1/12) * 8.5 * 100
DEMO_X_PEAK_MM3 = 2550.0     # = 212.5 * 12
DEMO_Y_BACK_MM3 = 243.75     # = 3.0 * 13 * (1/16) * 100
DEMO_Y_FRONT_MM3 = 3900.0    # = 243.75 * 16
DEMO_TOTAL_MM3 = 33150.0     # = 100 * 3.0 * 13 * 8.5


# ---------------------------------------------------------------------------
# Data model and checks
# ---------------------------------------------------------------------------

def demo_depth_mm(nx: int = DEMO_NX, ny: int = DEMO_NY,
                  peak_mm: float = DEMO_PEAK_MM) -> np.ndarray:
    """Closed-form synthetic wear field, shape (ny, nx).

    depth[y, x] = peak_mm * g[x]/max(g) * h[y]/max(h)

    g[x] is a symmetric triangular walking-lane weight (1 .. 12 .. 1 for nx=24)
    and h[y] = (ny - y)/ny is a linear decay from the front edge inward.
    The construction is separable, which makes both projections exact and
    hand-checkable.
    """
    i = np.arange(nx, dtype=float)
    g = np.where(i < nx // 2, i + 1.0, nx - i)
    g = g / g.max()

    j = np.arange(ny, dtype=float)
    h = (ny - j) / float(ny)

    return peak_mm * np.outer(h, g)


def compute_projections(depth_mm: np.ndarray, dx_mm: float,
                        dy_mm: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (x_wear_volume_mm3, y_wear_volume_mm3) from depth_mm[y, x]."""
    cell_area = float(dx_mm) * float(dy_mm)
    x_proj = np.asarray(depth_mm, dtype=float).sum(axis=0) * cell_area
    y_proj = np.asarray(depth_mm, dtype=float).sum(axis=1) * cell_area
    return x_proj, y_proj


def validate_depth(depth_mm: np.ndarray, *, negative_tol: float = 1e-9) -> np.ndarray:
    """Validate a wear-depth matrix and return a float copy with tiny negatives clipped."""
    depth = np.asarray(depth_mm, dtype=float)
    if depth.ndim != 2:
        raise ValueError(f"depth_mm must be 2-D (ny, nx); got shape {depth.shape}")
    ny, nx = depth.shape
    if ny < 2 or nx < 2:
        raise ValueError(f"depth_mm must have at least 2 rows and 2 columns; got {depth.shape}")
    if not np.all(np.isfinite(depth)):
        raise ValueError("depth_mm contains NaN or infinite values")
    if float(depth.min()) < -negative_tol:
        raise ValueError(
            "depth_mm contains negative wear depth (min = "
            f"{depth.min():.6g} mm); wear depth must be >= 0"
        )
    return np.clip(depth, 0.0, None)


def _loadtxt_tolerant(path: str | Path) -> np.ndarray:
    """loadtxt that accepts comma or whitespace and an optional header row."""
    with Path(path).open(encoding="utf-8-sig") as stream:
        for index, line in enumerate(stream):
            first = line.split("#", 1)[0].strip()
            if first:
                break
        else:
            raise ValueError(f"no matrix data in {path}")
    delimiter = "," if "," in first else None
    fields = [field.strip() for field in first.split(delimiter)]
    header = all(field and (field[0].isalpha() or field[0] == "_")
                 and field.lower() not in {"nan", "inf", "infinity"}
                 for field in fields)
    return np.loadtxt(path, delimiter=delimiter, skiprows=index + int(header),
                      dtype=float, ndmin=2, encoding="utf-8-sig")


def load_depth(path: str | Path) -> np.ndarray:
    """Load a depth matrix from .csv/.txt/.dat, .npy or .npz."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".npy":
        return np.asarray(np.load(path), dtype=float)
    if suffix == ".npz":
        with np.load(path) as bundle:
            key = "depth_mm" if "depth_mm" in bundle.files else bundle.files[0]
            return np.asarray(bundle[key], dtype=float)
    if suffix in {".csv", ".txt", ".dat", ""}:
        return _loadtxt_tolerant(path)
    raise ValueError(f"unsupported input type {suffix!r}; use .csv, .txt, .npy or .npz")


def _read_projection_csv(path: str | Path) -> np.ndarray:
    """Read a projection CSV; accept one column or (coordinate, value) columns."""
    arr = _loadtxt_tolerant(path)
    if arr.shape[1] == 1:
        return arr[:, 0]
    if arr.shape[1] == 2:
        return arr[:, 1]
    raise ValueError("projection CSV must have one or two columns")


def compare_projection(actual: np.ndarray, expected_path: str | Path,
                       label: str, rtol: float) -> dict:
    """Compare a computed projection against the program's CSV output."""
    expected = _read_projection_csv(expected_path)
    if expected.shape != actual.shape:
        raise ValueError(
            f"{label}: shape mismatch, computed {actual.shape}, file {expected.shape}"
        )
    denom = np.maximum(np.abs(expected), 1e-12)
    rel = np.abs(actual - expected) / denom
    max_rel = float(rel.max())
    ok = bool(np.allclose(actual, expected, rtol=rtol, atol=0.0))
    return {"label": label, "ok": ok, "max_rel_error": max_rel,
            "rtol": rtol, "n": int(actual.size)}


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def _cjk_font() -> str | None:
    """Return an installed CJK font name, or None for English-only labels."""
    try:
        from matplotlib import font_manager
    except Exception:  # pragma: no cover - matplotlib is a hard dependency
        return None
    names = {f.name for f in font_manager.fontManager.ttflist}
    for candidate in ("Noto Sans CJK SC", "Source Han Sans SC",
                      "WenQuanYi Zen Hei", "SimHei", "Microsoft YaHei"):
        if candidate in names:
            return candidate
    return None


def plot_figure(depth_mm: np.ndarray, x_proj: np.ndarray, y_proj: np.ndarray,
                *, dx_mm: float, dy_mm: float, x0_mm: float, y0_mm: float,
                outdir: Path, stem: str, source_label: str,
                synthetic: bool, editable_svg: bool, input_desc: str) -> list[Path]:
    """Render the heatmap and the two exactly-aligned marginal projections."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    from matplotlib.ticker import MaxNLocator, StrMethodFormatter

    cjk = _cjk_font()
    rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ([cjk, "DejaVu Sans"] if cjk else ["DejaVu Sans"]),
        "axes.unicode_minus": False,
        "font.size": 9,
        "axes.titlesize": 9.5,
        "axes.labelsize": 8.8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.6,
        "ytick.major.size": 2.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none" if editable_svg else "path",
    })

    def L(zh: str, en: str) -> str:
        return zh if cjk else en

    ny, nx = depth_mm.shape
    x_edges = x0_mm + np.arange(nx + 1) * dx_mm
    y_edges = y0_mm + np.arange(ny + 1) * dy_mm
    x_centers = x0_mm + (np.arange(nx) + 0.5) * dx_mm
    y_centers = y0_mm + (np.arange(ny) + 0.5) * dy_mm

    # --- dynamic layout: the heatmap box matches the physical data aspect ----
    data_aspect = (nx * dx_mm) / (ny * dy_mm)  # width / height
    heat_w_in = 5.8
    heat_h_in = heat_w_in / data_aspect
    if heat_h_in > 4.4:
        heat_h_in = 4.4
        heat_w_in = heat_h_in * data_aspect
    elif heat_h_in < 1.6:
        heat_h_in = 1.6
        heat_w_in = heat_h_in * data_aspect
    if heat_w_in > 8.0:
        heat_w_in = 8.0
        heat_h_in = heat_w_in / data_aspect

    left_in, bottom_in, top_in, right_in = 0.95, 0.72, 0.78, 1.05
    gap_xproj, xproj_h = 0.52, 1.05
    gap_yproj, yproj_w = 0.55, 1.85
    fig_w = left_in + heat_w_in + gap_yproj + yproj_w + right_in
    fig_h = top_in + heat_h_in + gap_xproj + xproj_h + bottom_in

    heat_bottom = (bottom_in + xproj_h + gap_xproj) / fig_h
    heat_h = heat_h_in / fig_h
    xleft = left_in / fig_w
    heat_w = heat_w_in / fig_w
    yleft = (left_in + heat_w_in + gap_yproj) / fig_w

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
    ax_h = fig.add_axes([xleft, heat_bottom, heat_w, heat_h])
    ax_x = fig.add_axes([xleft, bottom_in / fig_h, heat_w, xproj_h / fig_h],
                        sharex=ax_h)
    ax_y = fig.add_axes([yleft, heat_bottom, yproj_w / fig_w, heat_h],
                        sharey=ax_h)
    cax = fig.add_axes([yleft, (bottom_in + 0.42) / fig_h,
                        yproj_w / fig_w, 0.20 / fig_h])

    # --- heatmap -------------------------------------------------------------
    vmax = float(depth_mm.max())
    if vmax <= 0:
        vmax = 1.0
    mesh = ax_h.pcolormesh(x_edges, y_edges, depth_mm, cmap="cividis",
                           shading="flat", vmin=0.0, vmax=vmax,
                           rasterized=True)
    ax_h.set_xlim(x_edges[0], x_edges[-1])
    ax_h.set_ylim(y_edges[0], y_edges[-1])
    ax_h.set_ylabel(L("y (mm) — 踏面前后方向（前缘 → 内部）",
                      "y (mm) — front edge -> back"))
    ax_h.set_title(L("(a) 踏面磨损深度场  depth_mm[y][x]",
                     "(a) Tread wear depth field  depth_mm[y][x]"), loc="left", pad=4)
    ax_h.tick_params(axis="x", bottom=False, labelbottom=False)
    for side in ("top", "right"):
        ax_h.spines[side].set_color("#666666")
    ax_h.spines["bottom"].set_linewidth(1.1)  # the front edge / nosing
    ax_h.spines["bottom"].set_color("#222222")
    ax_h.xaxis.set_major_locator(MaxNLocator(6))
    ax_h.yaxis.set_major_locator(MaxNLocator(5))

    # origin and front-edge annotation
    ax_h.plot([x0_mm], [y0_mm], marker="o", ms=3.2, mfc="white", mec="#111111",
              mew=0.8, zorder=5, clip_on=False)
    ax_h.annotate(L(f"O({x0_mm:g},{y0_mm:g}) 前缘·左端",
                     f"O({x0_mm:g},{y0_mm:g}) front-left"),
                  xy=(x0_mm, y0_mm), xytext=(7, 7), textcoords="offset points",
                  fontsize=7.2, color="#111111",
                  bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#bbbbbb",
                            lw=0.4, alpha=0.88))

    total_mm3 = float(depth_mm.sum() * dx_mm * dy_mm)
    ax_h.text(0.985, 0.97,
              L(f"总磨损体积 ΣV = {total_mm3:,.0f} mm³",
                f"Total wear volume = {total_mm3:,.0f} mm³"),
              transform=ax_h.transAxes, ha="right", va="top", fontsize=7.8,
              color="#111111",
              bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#bbbbbb",
                        lw=0.4, alpha=0.88))

    if synthetic:
        from matplotlib import patheffects
        ax_h.text(0.5, 0.5, "SYNTHETIC", transform=ax_h.transAxes,
                  ha="center", va="center", rotation=28, fontsize=25,
                  fontweight="bold", color="white", alpha=0.75, zorder=6,
                  path_effects=[patheffects.withStroke(linewidth=2.4,
                                                       foreground="#333333")])

    # --- lateral projection (x) ---------------------------------------------
    ax_x.bar(x_centers, x_proj, width=0.82 * dx_mm, color="#1f5c8b",
             edgecolor="none", align="center", zorder=2)
    ax_x.set_ylabel("x_wear_volume_mm3\n(mm³)")
    ax_x.set_xlabel(L("x (mm) — 台阶左右方向（左 → 右）",
                      "x (mm) — lateral (left -> right)"))
    ax_x.set_title(L("(b) 横向磨损投影（沿 y 求和）",
                     "(b) Lateral projection (sum over y)"), loc="left", pad=4)
    ax_x.set_ylim(0, max(float(x_proj.max()), 1e-12) * 1.28)
    ax_x.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    ax_x.yaxis.set_major_locator(MaxNLocator(4))
    ax_x.grid(axis="y", color="#e8e8e8", lw=0.5)
    ax_x.set_axisbelow(True)
    ax_x.spines["top"].set_visible(False)
    ax_x.spines["right"].set_visible(False)
    ix = int(np.argmax(x_proj))
    ax_x.annotate(L(f"峰值 {x_proj[ix]:,.0f} mm³", f"peak {x_proj[ix]:,.0f} mm³"),
                  xy=(x_centers[ix], x_proj[ix]), xytext=(5, 4),
                  textcoords="offset points", fontsize=7.2, color="#123a55",
                  annotation_clip=False)

    # --- longitudinal projection (y) ----------------------------------------
    ax_y.barh(y_centers, y_proj, height=0.82 * dy_mm, color="#1f5c8b",
              edgecolor="none", align="center", zorder=2)
    ax_y.set_xlabel("y_wear_volume_mm3 (mm³)")
    ax_y.set_title(L("(c) 纵向磨损投影（沿 x 求和）",
                     "(c) Longitudinal projection (sum over x)"), loc="left", pad=4)
    ax_y.set_xlim(0, max(float(y_proj.max()), 1e-12) * 1.30)
    ax_y.xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    ax_y.xaxis.set_major_locator(MaxNLocator(3))
    ax_y.tick_params(labelleft=False)
    ax_y.grid(axis="x", color="#e8e8e8", lw=0.5)
    ax_y.set_axisbelow(True)
    ax_y.spines["top"].set_visible(False)
    ax_y.spines["right"].set_visible(False)
    iy = int(np.argmax(y_proj))
    ax_y.annotate(L(f"峰值 {y_proj[iy]:,.0f} mm³", f"peak {y_proj[iy]:,.0f} mm³"),
                  xy=(y_proj[iy], y_centers[iy]), xytext=(6, 0),
                  textcoords="offset points", fontsize=7.2, color="#123a55",
                  ha="left", va="center",
                  annotation_clip=False)

    # --- colour bar ----------------------------------------------------------
    cb = fig.colorbar(mesh, cax=cax, orientation="horizontal")
    cb.set_label(L("磨损深度 depth (mm)", "wear depth (mm)"))
    cb.set_ticks(np.linspace(0.0, vmax, 4))
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=7.4, width=0.4, length=2)

    # Bars on shared axes can trigger autoscaling; restore the exact domain.
    ax_x.set_xlim(x_edges[0], x_edges[-1])
    ax_y.set_ylim(y_edges[0], y_edges[-1])

    # --- provenance strip ----------------------------------------------------
    strip = (f"{input_desc}  ·  {nx} × {ny}  ·  "
             f"dx={dx_mm:g} mm, dy={dy_mm:g} mm")
    fig.text(0.008, 0.17 / fig_h, strip, fontsize=6.8, color="#555555",
             ha="left", va="bottom")
    fig.text(0.008, 0.04 / fig_h, source_label, fontsize=6.8, color="#555555",
             ha="left", va="bottom")

    if synthetic:
        fig.text(0.5, 1 - 0.16 / fig_h,
                 L("SYNTHETIC — 合成演示数据，非论文实测结果",
                   "SYNTHETIC — demo data, not a measured result"),
                 ha="center", va="top", fontsize=9.5, fontweight="bold",
                 color="#7a3e00",
                 bbox=dict(boxstyle="round,pad=0.32", fc="#fff3cd",
                           ec="#d39e00", lw=0.6))

    # --- save ----------------------------------------------------------------
    outdir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for ext in ("png", "pdf", "svg"):
        path = outdir / f"{stem}.{ext}"
        fig.savefig(path, dpi=300, facecolor="white")
        outputs.append(path)
    plt.close(fig)
    return outputs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Render a stair tread wear heatmap and its x/y projections.")
    src = p.add_mutually_exclusive_group(required=False)
    src.add_argument("--input", type=Path,
                     help="depth matrix (.csv/.txt: rows=y front->back, cols=x left->right; no header)")
    src.add_argument("--demo", action="store_true",
                     help="use the closed-form synthetic demo grid")
    p.add_argument("--outdir", type=Path, default=Path("figures/out/2501909"))
    p.add_argument("--stem", default=None,
                   help="output file stem (default: auto by synthetic flag)")
    p.add_argument("--dx-mm", type=float, default=DEMO_DX_MM)
    p.add_argument("--dy-mm", type=float, default=DEMO_DY_MM)
    p.add_argument("--x0-mm", type=float, default=0.0)
    p.add_argument("--y0-mm", type=float, default=0.0)
    p.add_argument("--source-label", default=None,
                   help="one-line provenance printed on the figure")
    p.add_argument("--check-x-proj", type=Path, default=None,
                   help="program output x_wear_volume_mm3.csv to cross-check")
    p.add_argument("--check-y-proj", type=Path, default=None,
                   help="program output y_wear_volume_mm3.csv to cross-check")
    p.add_argument("--rtol", type=float, default=1e-6)
    p.add_argument("--editable-svg", action="store_true",
                   help="keep SVG text as text (needs the same fonts to render)")
    p.add_argument("--self-test", action="store_true",
                   help="run closed-form projection checks only")
    args = p.parse_args(argv)

    if not args.self_test and not args.demo and args.input is None:
        p.error("choose --demo or --input <matrix>")
    if args.dx_mm <= 0 or args.dy_mm <= 0:
        p.error("--dx-mm and --dy-mm must be positive")
    return args


def run_self_test() -> int:
    """Closed-form checks for the demo and a uniform field."""
    ok = True

    depth = demo_depth_mm()
    x_proj, y_proj = compute_projections(depth, DEMO_DX_MM, DEMO_DY_MM)
    g = np.where(np.arange(DEMO_NX) < DEMO_NX // 2,
                 np.arange(DEMO_NX) + 1.0, DEMO_NX - np.arange(DEMO_NX))
    g = g / g.max()
    h = (DEMO_NY - np.arange(DEMO_NY)) / float(DEMO_NY)
    checks = [
        ("demo x projection", x_proj, DEMO_X_PEAK_MM3 * g),
        ("demo y projection", y_proj, DEMO_Y_FRONT_MM3 * h),
    ]

    for label, actual, expected in checks:
        good = bool(np.allclose(actual, expected, rtol=1e-12, atol=0.0))
        ok &= good
        print(f"{'PASS' if good else 'FAIL'} {label}: "
              f"max_abs_err={np.max(np.abs(actual - expected)):.3e}")
    for label, value, expected in [
        ("demo total (x sum)", float(x_proj.sum()), DEMO_TOTAL_MM3),
        ("demo total (y sum)", float(y_proj.sum()), DEMO_TOTAL_MM3),
        ("demo peak depth", float(depth.max()), DEMO_PEAK_MM),
        ("demo x edge (x=5 mm)", float(x_proj[0]), DEMO_X_EDGE_MM3),
        ("demo x peak (x=115 mm)", float(x_proj[11]), DEMO_X_PEAK_MM3),
        ("demo y front (y=5 mm)", float(y_proj[0]), DEMO_Y_FRONT_MM3),
        ("demo y back (y=155 mm)", float(y_proj[-1]), DEMO_Y_BACK_MM3),
    ]:
        good = bool(np.isclose(value, expected, rtol=1e-12, atol=0.0))
        ok &= good
        print(f"{'PASS' if good else 'FAIL'} {label}: {value:.6f} "
              f"(expected {expected:.6f})")

    # a uniform 4x3 field has a trivial hand solution
    n_y, n_x, d, dx, dy = 3, 4, 2.0, 3.0, 5.0
    uniform = np.full((n_y, n_x), d)
    ux, uy = compute_projections(uniform, dx, dy)
    for label, value, expected in [
        ("uniform x column", float(ux[0]), d * n_y * dx * dy),
        ("uniform y row", float(uy[0]), d * n_x * dx * dy),
        ("uniform total", float(ux.sum()), d * n_x * n_y * dx * dy),
    ]:
        good = bool(np.isclose(value, expected, rtol=1e-12, atol=0.0))
        ok &= good
        print(f"{'PASS' if good else 'FAIL'} {label}: {value:.6f} "
              f"(expected {expected:.6f})")

    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    if args.self_test:
        return run_self_test()

    if args.demo:
        depth = demo_depth_mm()
        dx_mm, dy_mm = DEMO_DX_MM, DEMO_DY_MM
        synthetic = True
        input_desc = "demo: triangular(x) × linear(y)"
        default_stem = "fig_2501909_tread_wear_synthetic"
        default_source = "Closed-form demo. Peak depth: 3 mm."
        input_info = {
            "kind": "closed_form_demo",
            "generator": "demo_depth_mm",
            "parameters": {"nx": DEMO_NX, "ny": DEMO_NY,
                           "peak_mm": DEMO_PEAK_MM},
        }
    else:
        depth = load_depth(args.input)
        dx_mm, dy_mm = args.dx_mm, args.dy_mm
        synthetic = False
        input_desc = f"input: {Path(args.input).name}"
        default_stem = "fig_2501909_tread_wear"
        default_source = "Input grid. See figure_manifest.json for provenance."
        digest = hashlib.sha256(Path(args.input).read_bytes()).hexdigest()
        input_info = {
            "kind": "matrix",
            "path": str(args.input),
            "sha256": digest,
            "shape_ny_nx": list(np.asarray(depth).shape),
        }

    depth = validate_depth(depth)
    x_proj, y_proj = compute_projections(depth, dx_mm, dy_mm)

    checks = []
    for label, actual, path in (
        (PROJ_X_NAME, x_proj, args.check_x_proj),
        (PROJ_Y_NAME, y_proj, args.check_y_proj),
    ):
        if path is not None:
            result = compare_projection(actual, path, label, args.rtol)
            checks.append(result)
            if not result["ok"]:
                print(f"[FAIL] {label}: max_rel_error="
                      f"{result['max_rel_error']:.3e} > rtol={args.rtol:g}",
                      file=sys.stderr)
                return 2

    stem = args.stem or default_stem
    source_label = args.source_label or default_source
    outdir = args.outdir
    outputs = plot_figure(
        depth, x_proj, y_proj,
        dx_mm=dx_mm, dy_mm=dy_mm, x0_mm=args.x0_mm, y0_mm=args.y0_mm,
        outdir=outdir, stem=stem, source_label=source_label,
        synthetic=synthetic, editable_svg=args.editable_svg,
        input_desc=input_desc,
    )

    x_csv = outdir / f"{PROJ_X_NAME}.csv"
    y_csv = outdir / f"{PROJ_Y_NAME}.csv"
    np.savetxt(x_csv, np.column_stack([args.x0_mm + (np.arange(depth.shape[1]) + 0.5) * dx_mm,
                                       x_proj]),
               delimiter=",", header="x_mm,x_wear_volume_mm3", comments="", fmt="%.10g")
    np.savetxt(y_csv, np.column_stack([args.y0_mm + (np.arange(depth.shape[0]) + 0.5) * dy_mm,
                                       y_proj]),
               delimiter=",", header="y_mm,y_wear_volume_mm3", comments="", fmt="%.10g")

    import matplotlib
    manifest = {
        "figure_stem": stem,
        "synthetic": synthetic,
        "source_label": source_label,
        "input": input_info,
        "grid": {
            "ny": int(depth.shape[0]),
            "nx": int(depth.shape[1]),
            "dx_mm": dx_mm,
            "dy_mm": dy_mm,
            "x0_mm": args.x0_mm,
            "y0_mm": args.y0_mm,
            "axes": {
                "x": "lateral, 0 = left edge, increases rightward",
                "y": "longitudinal, 0 = front edge (nosing), increases inward",
            },
            "units": {"depth": "mm", "cell_area": "mm^2", "projection": "mm^3"},
        },
        "projections": {
            "x_wear_volume_mm3_total": float(x_proj.sum()),
            "y_wear_volume_mm3_total": float(y_proj.sum()),
            "grand_total_mm3": float(depth.sum() * dx_mm * dy_mm),
            "max_depth_mm": float(depth.max()),
            "argmax_depth_xy_mm": [
                float(args.x0_mm + (int(np.unravel_index(np.argmax(depth), depth.shape)[1]) + 0.5) * dx_mm),
                float(args.y0_mm + (int(np.unravel_index(np.argmax(depth), depth.shape)[0]) + 0.5) * dy_mm),
            ],
        },
        "cross_checks": checks,
        "outputs": [str(p) for p in outputs] + [str(x_csv), str(y_csv)],
        "software": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "generated_utc": datetime.now(timezone.utc).isoformat(),
    }
    manifest_path = outdir / "figure_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False),
                             encoding="utf-8")

    print(f"[OK] grid {depth.shape[0]}x{depth.shape[1]}, dx={dx_mm:g} mm, "
          f"dy={dy_mm:g} mm, max depth {depth.max():.6g} mm"
          + ("  [SYNTHETIC]" if synthetic else ""))
    print(f"[OK] total volume {manifest['projections']['grand_total_mm3']:.6f} mm^3 "
          f"(x-sum {x_proj.sum():.6f}, y-sum {y_proj.sum():.6f})")
    for c in checks:
        print(f"[OK] {c['label']} cross-check: max_rel_error={c['max_rel_error']:.3e}")
    for p in outputs:
        print(f"[OK] wrote {p}")
    print(f"[OK] wrote {x_csv}")
    print(f"[OK] wrote {y_csv}")
    print(f"[OK] wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
