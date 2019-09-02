import os
import unittest
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from bot_eo_analyses import GetAmplifierGains

class GetAmplifierGainsTestCase(unittest.TestCase):
    """TestCase subclass for GetAmplifierGains."""
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_get_amplifier_gains(self):
        """Test the returned gain values against values read directly
        from the eT web interface."""
        det_name = 'R10_S10'
        bot_eo_config_file = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                                          'tests', 'data', 'bot_eo_config_file')
        get_amplifier_gains = GetAmplifierGains(bot_eo_config_file)
        gains = get_amplifier_gains('6549D_R10_S10_eotest_results.fits')
        self.assertAlmostEqual(gains[4], 0.9903897047042847)
        self.assertAlmostEqual(gains[7], 0.9513707160949707)
        self.assertAlmostEqual(gains[16], 0.961654007434845)


if __name__ == '__main__':
    unittest.main()

