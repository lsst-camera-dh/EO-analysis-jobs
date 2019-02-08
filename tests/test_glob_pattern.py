"""
Unit tests for get_glob_patterns function.
"""
import unittest
from bot_eo_analyses import GlobPattern

class GlobPatternTestCase(unittest.TestCase):
    """TestCase class for GlobPattern class."""
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_glob_pattern(self):
        """
        Test for GlobPattern class.
        """
        glob_pattern = GlobPattern()
        det_name = 'R22_S11'
        self.assertEqual(glob_pattern('fe55', det_name),
                         'fe55_fe55_*/*_{}.fits'.format(det_name))
        self.assertEqual(glob_pattern('cte_low', det_name),
                         'sflat_*_flat_*L*/*_{}.fits'.format(det_name))

if __name__ == '__main__':
    unittest.main()
