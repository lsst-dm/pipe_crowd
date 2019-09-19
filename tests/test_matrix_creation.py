
import unittest
import numpy as np
import lsst.utils.tests
from lsst.afw.image import ExposureF
from lsst.pipe.crowd import CrowdedFieldMatrix
from lsst.meas.algorithms.installGaussianPsf import InstallGaussianPsfTask, InstallGaussianPsfConfig
import lsst.pex.exceptions.wrappers
from collections import Counter

class CrowdedFieldMatrixTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.exposure = ExposureF(1000, 1000)
        psfConfig =InstallGaussianPsfConfig()
        psfConfig.fwhm = 4
        psfTask = InstallGaussianPsfTask(config=psfConfig)
        psfTask.run(exposure=self.exposure)

    def test_singleSource(self):

        matrix = CrowdedFieldMatrix(self.exposure,
                                    np.array([200.0]),
                                    np.array([200.0]))
        self.assertGreater(len(matrix.getMatrixEntries()), 0)

    def test_multipleSources(self):

        x_arr = np.linspace(20.0, 60.0, 20)
        y_arr = np.linspace(30.0, 70.0, 20)
        matrix = CrowdedFieldMatrix(self.exposure, x_arr, y_arr)

        self.assertGreater(len(matrix.getMatrixEntries()), 0)

    def test_unequalLengths(self):

        x_arr = np.linspace(20.0, 60.0, 20)
        y_arr = np.linspace(30.0, 70.0, 21)  #  Note 21 != 20
        with self.assertRaises(lsst.pex.exceptions.wrappers.LengthError):
            matrix = CrowdedFieldMatrix(self.exposure, x_arr, y_arr)

    def test_renameMatrixRows(self):
        matrix = CrowdedFieldMatrix(self.exposure,
                                    np.array([200.0, 200.0]),
                                    np.array([204.0, 204.0]))

        pre_rename_entries = matrix.getMatrixEntries()
        matrix.renameMatrixRows()
        post_rename_entries = matrix.getMatrixEntries()


        self.assertEqual(len(pre_rename_entries), len(post_rename_entries))

        # We check that the count of rows that appear a given number of times
        # remains the same, since this is invariant under renaming.
        pre_counter = Counter(x[1] for x in pre_rename_entries)
        post_counter = Counter(x[1] for x in post_rename_entries)

        pre_counts = list(pre_counter.values())
        post_counts = list(post_counter.values())
        pre_counts.sort()
        post_counts.sort()
        self.assertListEqual(pre_counts, post_counts)


    def test_solve(self):
        matrix = CrowdedFieldMatrix(self.exposure,
                                    np.array([200.0, 200.0]),
                                    np.array([204.0, 204.0]))

        matrix.solve()


