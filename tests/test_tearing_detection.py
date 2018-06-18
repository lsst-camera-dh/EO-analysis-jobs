"""
Unit tests for tearing detection code.
"""
import os
import unittest
from tearing_detection import tearing_detection, persist_tearing_png_files

class TearingDetectionTestCase(unittest.TestCase):
    """
    TestCase class for tearing_detection module.
    """
    def setUp(self):
        self.flat_with_tearing = 'E2V-CCD250-229-Dev_flat_with_tearing.fits.gz'
        self.superflat_without_tearing \
            = 'E2V-CCD250-160_sflat_without_tearing.fits.gz'
        self.png_filename \
            = self.flat_with_tearing.split('.')[0] + '_tearing.png'

    def tearDown(self):
        if os.path.isfile(self.png_filename):
            os.remove(self.png_filename)

    def test_tearing_detection(self):
        """
        Test for tearing_detection and persist_tearing_png_files functions.
        """
        fits_files = [os.path.join(os.environ['EOANALYSISJOBSDIR'], 'tests',
                                   'data', x) for x in
                      (self.flat_with_tearing, self.superflat_without_tearing)]
        files_with_tearing, png_files = tearing_detection(fits_files)
        # Check that png_files is empty.
        self.assertFalse(png_files)
        # Check that the expected filename is returned.
        self.assertEqual(files_with_tearing[0], fits_files[0])
        self.assertEqual(len(files_with_tearing), 1)

        # Check output if make_png_files is True.
        files_with_tearing, png_files \
            = tearing_detection(fits_files, make_png_files=True)
        # Check that the expected filenames are returned.
        self.assertListEqual(files_with_tearing, [fits_files[0]])
        self.assertListEqual(png_files, [self.png_filename])

        # Test the persistence function for the png files.
        file_refs = persist_tearing_png_files(png_files)
        self.assertEqual(len(file_refs), 1)
        self.assertEqual(file_refs[0]['path'], png_files[0])
        md = eval(file_refs[0]['metadata'])
        self.assertEqual(md['DATA_PRODUCT'], 'tearing_profiles')
        self.assertEqual(md['LsstId'], 'E2V-CCD250-229-Dev')

if __name__ == '__main__':
    unittest.main()
