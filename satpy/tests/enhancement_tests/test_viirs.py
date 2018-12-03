"""Unit testing for the VIIRS enhancement function
"""

import unittest
import numpy as np
import xarray as xr
from .test_enhancements import TestEnhancementStretch


class TestVIIRSEnhancement(unittest.TestCase):
    """Class for testing the VIIRS enhancement function in satpy.enhancements.viirs"""

    def setUp(self):
        """Setup the test"""
        data = np.arange(15, 301, 15).reshape(2, 10)
        self.da = xr.DataArray(data, dims=('y', 'x'), attrs={'test': 'test'})
        self.palette = {'colors':
                        [[14, [0.0, 0.0, 0.0]],
                         [15, [0.0, 0.0, 0.39215]],
                         [16, [0.76862, 0.63529, 0.44705]],
                         [17, [0.76862, 0.63529, 0.44705]],
                         [18, [0.0, 0.0, 1.0]],
                         [20, [1.0, 1.0, 1.0]],
                         [27, [0.0, 1.0, 1.0]],
                         [30, [0.78431, 0.78431, 0.78431]],
                         [31, [0.39215, 0.39215, 0.39215]],
                         [88, [0.70588, 0.0, 0.90196]],
                         [100, [0.19607, 1.0, 0.39215]],
                         [120, [0.19607, 1.0, 0.39215]],
                         [121, [0.0, 1.0, 0.0]],
                         [130, [0.0, 1.0, 0.0]],
                         [131, [0.78431, 1.0, 0.0]],
                         [140, [0.78431, 1.0, 0.0]],
                         [141, [1.0, 1.0, 0.58823]],
                         [150, [1.0, 1.0, 0.58823]],
                         [151, [1.0, 1.0, 0.0]],
                         [160, [1.0, 1.0, 0.0]],
                         [161, [1.0, 0.78431, 0.0]],
                         [170, [1.0, 0.78431, 0.0]],
                         [171, [1.0, 0.58823, 0.19607]],
                         [180, [1.0, 0.58823, 0.19607]],
                         [181, [1.0, 0.39215, 0.0]],
                         [190, [1.0, 0.39215, 0.0]],
                         [191, [1.0, 0.0, 0.0]],
                         [200, [1.0, 0.0, 0.0]],
                         [201, [0.0, 0.0, 0.0]]],
                        'min_value': 0,
                        'max_value': 201}

    def test_viirs(self):
        from satpy.enhancements.viirs import water_detection
        expected = [[[1, 7, 8, 8, 8, 9, 10, 11, 14, 8],
                     [20, 23, 26, 10, 12, 15, 18, 21, 24, 27]]]
        TestEnhancementStretch._test_enhancement(self, water_detection, self.da,
                                                 expected, palettes=self.palette)

    def tearDown(self):
        """Clean up"""
        pass


def suite():
    """The test suite for test_viirs.
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestVIIRSEnhancement))

    return mysuite
