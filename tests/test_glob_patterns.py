"""
Unit tests for get_glob_patterns function.
"""
import unittest
from bot_eo_analyses import get_glob_patterns

class GlobPatternsTestCase(unittest.TestCase):
    """TestCase class for get_glob_patterns function."""
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_glob_pattenrs(self):
        """
        Test for get_glob_patterns function.
        """
        patterns = get_glob_patterns()

        self.assertEqual(patterns['fe55'], 'fe55_fe55_*')
        self.assertEqual(patterns['cte_low'], 'sflat_*_flat_*L*')

if __name__ == '__main__':
    unittest.main()
