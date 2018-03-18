"""
Unit tests for copy_fe55_gain_results.
"""
import os
import unittest
from copy_fe55_gain_results import copy_fe55_gain_results

class CopyFe55GainResultsTestCase(unittest.TestCase):
    "Test case class for copy_fe55_gain_results."
    def setUp(self):
        os.environ['SITENAME'] = 'SLAC'
        os.environ['LCATR_LIMS_URL'] = 'Prod'
        os.environ['LCATR_UNIT_ID'] = 'LCA-11021_RTM-002_ETU1'
        os.environ['LCATR_UNIT_TYPE'] = 'LCA-11021_RTM'

    def tearDown(self):
        pass

    @unittest.skipUnless(os.path.isdir('/gpfs/slac/lsst/fs1/g/data/jobHarness/jh_archive/LCA-11021_RTM/LCA-11021_RTM-002_ETU1/7167/fe55_raft_analysis/v0/45541'),
                         "eotest results files directory not available")
    def test_copy_fe55_gain_results(self):
        run = 7167
        files = copy_fe55_gain_results(run)
        self.assertTrue('ITL-3800C-139_eotest_results.fits' in files)
        for item in files:
            os.remove(item)

        self.assertRaises(RuntimeError, copy_fe55_gain_results, run,
                          num_ccds=10)
