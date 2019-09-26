
import unittest
import numpy as np
import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.geom as geom
from lsst.afw.image import ExposureF
from lsst.pipe.crowd import CrowdedFieldMatrix
from lsst.meas.algorithms.installGaussianPsf import InstallGaussianPsfTask, InstallGaussianPsfConfig
import lsst.pex.exceptions.wrappers
from collections import Counter
from lsst.geom import Point2D

def add_psf_image(exposure, x, y, flux):

    psfImg = exposure.getPsf().computeImage(Point2D(x, y))
    psfImg *= flux

    subim = afwImage.ImageF(exposure.getMaskedImage().getImage(),
                            psfImg.getBBox(), afwImage.LOCAL)
    subim += psfImg.convertF()


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

    def test_createFromCatalog(self):

        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField("centroid_x", type=np.float64)
        schema.addField("centroid_y", type=np.float64)
        flux_key = schema.addField("flux_flux", type=np.float32)
        schema.getAliasMap().set("slot_Centroid", "centroid")
        testCatalog = afwTable.SourceCatalog(schema)
        for n in range(20):
            r = testCatalog.addNew()
            r["centroid_x"] = n*3
            r["centroid_y"] = n*3

        matrix = CrowdedFieldMatrix(self.exposure, testCatalog, flux_key)

        self.assertGreater(len(matrix.getMatrixEntries()), 0)

    @unittest.skip
    def test_renameMatrixRows(self):
        #
        # Haven't figured out how to get access to pre-renaming matrix
        # to do the comparision with the post-renaming one.
        #
        pre_rename_entries = matrix.getMatrixEntries()
        matrix.renameMatrixRows()

        matrix = CrowdedFieldMatrix(self.exposure,
                                    np.array([200.0, 200.0]),
                                    np.array([204.0, 204.0]))
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
        exposure = ExposureF(1000, 1000)
        psfConfig = InstallGaussianPsfConfig()
        psfConfig.fwhm = 4
        psfTask = InstallGaussianPsfTask(config=psfConfig)
        psfTask.run(exposure=exposure)

        add_psf_image(exposure, 200.0, 400.0, 600.0)
        add_psf_image(exposure, 210.0, 210.0, 300.0)

        matrix = CrowdedFieldMatrix(exposure,
                                    np.array([200.0, 210.0]),
                                    np.array([400.0, 210.0]))
        result = matrix.solve()

        self.assertFloatsAlmostEqual(result, np.array([600.0, 300.0]), atol=1e-3);

    def test_solve_catalog(self):
        exposure = ExposureF(1000, 1000)
        psfConfig = InstallGaussianPsfConfig()
        psfConfig.fwhm = 4
        psfTask = InstallGaussianPsfTask(config=psfConfig)
        psfTask.run(exposure=exposure)

        add_psf_image(exposure, 200.0, 400.0, 600.0)
        add_psf_image(exposure, 210.0, 210.0, 300.0)

        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField("centroid_x", type=np.float64)
        schema.addField("centroid_y", type=np.float64)
        flux_key = schema.addField("flux_flux", type=np.float32)
        schema.getAliasMap().set("slot_Centroid", "centroid")
        testCatalog = afwTable.SourceCatalog(schema)
        for x, y in zip([200.0, 210.0], [400.0, 210]):
            r = testCatalog.addNew()
            r["centroid_x"] = x
            r["centroid_y"] = y

        matrix = CrowdedFieldMatrix(exposure, testCatalog, flux_key)
        result = matrix.solve()

        self.assertFloatsAlmostEqual(testCatalog["flux_flux"], np.array([600.0, 300.0]), atol=1e-3);

    @unittest.skip
    def test_solve_subimage(self):
        exposure = ExposureF(1000, 1000)
        psfConfig = InstallGaussianPsfConfig()
        psfConfig.fwhm = 4
        psfTask = InstallGaussianPsfTask(config=psfConfig)
        psfTask.run(exposure=exposure)

        self.add_psf_image(exposure, 200.0, 400.0, 600.0)
        self.add_psf_image(exposure, 210.0, 210.0, 300.0)

        matrix = CrowdedFieldMatrix(exposure,
                                    np.array([200.0, 210.0]),
                                    np.array([400.0, 210.0]))
        result = matrix.solve()

        self.assertFloatsAlmostEqual(result, np.array([600.0, 300.0]), atol=1e-3);


