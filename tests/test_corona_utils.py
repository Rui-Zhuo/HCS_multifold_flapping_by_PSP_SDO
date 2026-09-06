import unittest

import numpy as np

from Corona_image_process.CoronaImageProcess_utils import (
    build_atrous_coef,
    insert_nan_columns,
    radial_slit,
)


class CoronaUtilityTests(unittest.TestCase):
    def test_level_zero_atrous_kernels_are_normalized(self) -> None:
        for method, size in (("B_spline", 5), ("linear", 3)):
            with self.subTest(method=method):
                kernel = build_atrous_coef(0, method)
                self.assertEqual(kernel.shape, (size, size))
                self.assertTrue(np.isclose(kernel.sum(), 1.0))

    def test_radial_slit_rejects_zero_length(self) -> None:
        with self.assertRaisesRegex(ValueError, "must differ"):
            radial_slit((1, 1), (1, 1), 0, 0)

    def test_insert_nan_columns_does_not_mutate_input(self) -> None:
        source = np.ones((2, 2))
        result = insert_nan_columns(source, [(1, 2)])
        self.assertEqual(result.shape, (2, 4))
        self.assertTrue(np.isnan(result[:, 1:3]).all())
        self.assertEqual(source.shape, (2, 2))


if __name__ == "__main__":
    unittest.main()
