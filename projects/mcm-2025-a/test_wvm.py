import json
from pathlib import Path
import unittest

from wvm import infer, volume_per_passage_mm3, wear_volume_mm3


class WearVolumeTests(unittest.TestCase):
    def setUp(self):
        path = Path(__file__).parent / "examples" / "synthetic-grid.json"
        self.inputs = json.loads(path.read_text(encoding="utf-8"))

    def test_grid_to_passage_rate(self):
        result = infer(**self.inputs)
        self.assertAlmostEqual(result["wear_volume_mm3"], 6.0)
        self.assertAlmostEqual(result["volume_per_passage_mm3"], 0.1)
        self.assertAlmostEqual(result["total_passages"], 60.0)
        self.assertAlmostEqual(result["daily_passages"], 6.0)

    def test_grid_to_age_given_daily_rate(self):
        del self.inputs["age_days"]
        self.assertAlmostEqual(infer(**self.inputs, daily_passages=6)["age_days"], 10.0)

    def test_zero_wear_with_known_age(self):
        self.inputs["depth_mm"] = [[0, 0]]
        self.assertEqual(infer(**self.inputs)["daily_passages"], 0.0)

    def test_invalid_grids(self):
        for grid in ([], [[1], [2, 3]], [[-1]], [[float("nan")]], [[False]]):
            with self.subTest(grid=grid), self.assertRaises(ValueError):
                wear_volume_mm3(grid, 1)

    def test_requires_exactly_one_known_quantity(self):
        with self.assertRaises(ValueError):
            infer(**self.inputs, daily_passages=6)
        del self.inputs["age_days"]
        with self.assertRaises(ValueError):
            infer(**self.inputs)

    def test_table_four_units(self):
        # Original paper §3.4 Table 4 uses N, mm and N/mm².
        per_passage = volume_per_passage_mm3(700, 1.8e-7, 0.8, 95)
        self.assertAlmostEqual(per_passage, 1.0610526315789474e-6, places=15)
        self.assertAlmostEqual(per_passage * 260 * 36492, 10.067182484210526)


if __name__ == "__main__":
    unittest.main()
