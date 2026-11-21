import unittest

from wdm import wear_marginals


class WearMarginalTests(unittest.TestCase):
    def test_axes_and_conservation_on_rectangular_grid(self):
        result = wear_marginals([[0, 1, 2], [4, 3, 0]], 2)
        self.assertEqual(result["wear_volume_mm3"], 20)
        self.assertEqual(result["x_wear_volume_mm3"], [8, 8, 4])
        self.assertEqual(result["y_wear_volume_mm3"], [6, 14])
        self.assertEqual(result["x_wear_fraction"], [0.4, 0.4, 0.2])
        self.assertEqual(result["y_wear_fraction"], [0.3, 0.7])

    def test_fractions_are_independent_of_footprint_cell_area(self):
        grid = [[1, 4], [2, 3]]
        reference = wear_marginals(grid, 1)
        scaled = wear_marginals(grid, 5)
        self.assertEqual(scaled["wear_volume_mm3"], reference["wear_volume_mm3"] * 5)
        for axis in ("x", "y"):
            key = f"{axis}_wear_fraction"
            self.assertEqual(scaled[key], reference[key])

    def test_zero_grid_has_no_normalized_distribution(self):
        result = wear_marginals([[0, 0, 0], [0, 0, 0]], 2)
        self.assertEqual(result["x_wear_volume_mm3"], [0, 0, 0])
        self.assertEqual(result["y_wear_volume_mm3"], [0, 0])
        self.assertIsNone(result["x_wear_fraction"])
        self.assertIsNone(result["y_wear_fraction"])

    def test_distinct_spatial_patterns_can_share_both_marginals(self):
        diagonal = wear_marginals([[1, 0], [0, 1]], 1)
        anti_diagonal = wear_marginals([[0, 1], [1, 0]], 1)
        self.assertEqual(diagonal, anti_diagonal)
