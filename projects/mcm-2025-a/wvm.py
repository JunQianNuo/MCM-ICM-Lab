"""Constant-hardness wear volume model from MCM paper 2501909, §3.1."""

import argparse
import json
import math
from pathlib import Path


def positive(value, name):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
    ):
        raise ValueError(f"{name} must be a finite positive number")
    return value


def wear_volume_mm3(depth_mm, cell_area_mm2):
    """Integrate mean depth in each equal-area cell of a rectangular grid."""
    positive(cell_area_mm2, "cell_area_mm2")
    if not isinstance(depth_mm, list) or not depth_mm:
        raise ValueError("depth_mm must be a nonempty array of rows")
    first_row = depth_mm[0]
    if not isinstance(first_row, list) or not first_row:
        raise ValueError("depth_mm rows must be nonempty arrays")
    columns = len(first_row)
    for row in depth_mm:
        if not isinstance(row, list) or len(row) != columns:
            raise ValueError("depth_mm must be a rectangular grid")
        for depth in row:
            if depth != 0 or isinstance(depth, bool):
                positive(depth, "wear depth")
    volume = math.fsum(depth for row in depth_mm for depth in row) * cell_area_mm2
    if not math.isfinite(volume):
        raise ValueError("wear volume exceeds the numeric range")
    return volume


def volume_per_passage_mm3(load_n, wear_coefficient, sliding_mm, hardness_n_per_mm2):
    """Return Archard wear volume for one footfall on this tread."""
    positive(load_n, "load_n")
    positive(wear_coefficient, "wear_coefficient")
    positive(sliding_mm, "sliding_mm")
    positive(hardness_n_per_mm2, "hardness_n_per_mm2")
    volume = wear_coefficient * load_n * sliding_mm / hardness_n_per_mm2
    return positive(volume, "volume per passage")


def infer(depth_mm, cell_area_mm2, archard, *, age_days=None, daily_passages=None):
    """Estimate daily tread passages or age given the other quantity."""
    if (age_days is None) == (daily_passages is None):
        raise ValueError("supply exactly one of age_days or daily_passages")
    volume = wear_volume_mm3(depth_mm, cell_area_mm2)
    per_passage = volume_per_passage_mm3(**archard)
    passages = volume / per_passage
    if age_days is None:
        age_days = passages / positive(daily_passages, "daily_passages")
    else:
        daily_passages = passages / positive(age_days, "age_days")
    return {
        "wear_volume_mm3": volume,
        "volume_per_passage_mm3": per_passage,
        "total_passages": passages,
        "age_days": age_days,
        "daily_passages": daily_passages,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON with a depth grid and model inputs")
    args = parser.parse_args()
    try:
        inputs = json.loads(args.input.read_text(encoding="utf-8"))
        result = infer(**inputs)
        print(json.dumps(result, indent=2, allow_nan=False))
    except (OSError, TypeError, ValueError, ArithmeticError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
