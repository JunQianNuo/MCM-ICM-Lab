"""Compute measured wear marginals from MCM paper 2501909, Eq. (3.17)."""

import argparse
import json
import math
from pathlib import Path

from wvm import wear_volume_mm3


def wear_marginals(depth_mm, cell_area_mm2):
    """Project a grid indexed by [y][x] into lateral and front-to-back wear."""
    total = wear_volume_mm3(depth_mm, cell_area_mm2)
    x_volume = [math.fsum(column) * cell_area_mm2 for column in zip(*depth_mm)]
    y_volume = [math.fsum(row) * cell_area_mm2 for row in depth_mm]
    return {
        "wear_volume_mm3": total,
        "x_wear_volume_mm3": x_volume,
        "y_wear_volume_mm3": y_volume,
        "x_wear_fraction": [volume / total for volume in x_volume] if total else None,
        "y_wear_fraction": [volume / total for volume in y_volume] if total else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON with depth_mm and cell_area_mm2")
    args = parser.parse_args()
    try:
        inputs = json.loads(args.input.read_text(encoding="utf-8"))
        result = wear_marginals(inputs["depth_mm"], inputs["cell_area_mm2"])
        print(json.dumps(result, indent=2, allow_nan=False))
    except (OSError, KeyError, TypeError, ValueError, ArithmeticError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
