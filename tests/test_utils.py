from datetime import datetime
import unittest

from hcs_flapping.utils import iter_cadence


class UtilityTests(unittest.TestCase):
    def test_iter_cadence_is_inclusive(self) -> None:
        start = datetime(2021, 1, 17, 0, 0)
        end = datetime(2021, 1, 17, 0, 24)
        self.assertEqual(
            list(iter_cadence(start, end, 12)),
            [
                datetime(2021, 1, 17, 0, 0),
                datetime(2021, 1, 17, 0, 12),
                datetime(2021, 1, 17, 0, 24),
            ],
        )

    def test_iter_cadence_rejects_non_positive_step(self) -> None:
        with self.assertRaisesRegex(ValueError, "positive"):
            list(iter_cadence(datetime(2021, 1, 17), datetime(2021, 1, 17), 0))


if __name__ == "__main__":
    unittest.main()
