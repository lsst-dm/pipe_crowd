

import numpy as np
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
from lsst.geom import Extent2I
import lsst.afw.detection as afwDetection

from lsst.meas.base import SdssCentroidAlgorithm, SdssCentroidControl
from lsst.meas.base import MeasurementError

class CrowdedCentroidTaskConfig(pexConfig.Config):
    """Config for CrowdedCentroidTask"""

    pass

class CrowdedCentroidTask(pipeBase.Task):
    ConfigClass = CrowdedCentroidTaskConfig
    _DefaultName = "crowdedCentroidTask"

    def __init__(self, schema, **kwargs):
        self.centroid_control = SdssCentroidControl()
        self.sdssCentroid = SdssCentroidAlgorithm(self.centroid_control,
                                                  "centroid",
                                                  schema)
        pipeBase.CmdLineTask.__init__(self, **kwargs)

    @pipeBase.timeMethod
    def run(self, exposure, catalog, flux_key):

        subtracted_exposure = afwImage.ExposureF(exposure, deep=True)
        subtracted_image = subtracted_exposure.getMaskedImage()

        #
        # Totally a repetition from subtraction.py, don't leave this here
        #
        for source in catalog:
            centroid = source.getCentroid()
            psf_image = exposure.getPsf().computeImage(centroid)
            psf_image *= source[flux_key]
            bbox = psf_image.getBBox()
            bbox.clip(subtracted_image.getBBox())
            image_subregion = afwImage.ImageF(subtracted_image.getImage(),
                                    bbox, afwImage.LOCAL)
            image_subregion -= psf_image[bbox].convertF()

        for source in catalog:
            centroid = source.getCentroid()
            psf_image = exposure.getPsf().computeImage(centroid)
            psf_image *= source[flux_key]
            bbox = psf_image.getBBox()
            bbox.clip(subtracted_image.getBBox())

            image_subregion = afwImage.ImageF(subtracted_image.getImage(),
                                         bbox, afwImage.LOCAL)
            image_subregion += psf_image[bbox].convertF()

            spanSet = afwGeom.SpanSet.fromShape(5)
            spanSet = spanSet.shiftedBy(Extent2I(source.getCentroid()))
            footprint = afwDetection.Footprint(spanSet)
            peak = footprint.peaks.addNew()
            peak.setFx(source.getCentroid().getX())
            peak.setFy(source.getCentroid().getY())
            peak.setIx(int(source.getCentroid().getX()))
            peak.setIy(int(source.getCentroid().getY()))
            source.setFootprint(footprint)

            try:
                self.sdssCentroid.measure(source, exposure)
            except (MeasurementError):
                pass


            image_subregion -= psf_image[bbox].convertF()





