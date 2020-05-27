

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

from .modelImage import ModelImageTask

class CrowdedCentroidTaskConfig(pexConfig.Config):
    """Config for CrowdedCentroidTask"""

    modelImage = pexConfig.ConfigurableField(
            target=ModelImageTask,
            doc="Tools for manipulating model images"
    )

class CrowdedCentroidTask(pipeBase.Task):
    ConfigClass = CrowdedCentroidTaskConfig
    _DefaultName = "crowdedCentroidTask"

    def __init__(self, schema, **kwargs):
        pipeBase.CmdLineTask.__init__(self, **kwargs)
        self.centroid_control = SdssCentroidControl()
        self.sdssCentroid = SdssCentroidAlgorithm(self.centroid_control,
                                                  "centroid",
                                                  schema)
        self.makeSubtask("modelImage")

    @pipeBase.timeMethod
    def run(self, exposure, catalog, flux_key):

        subtracted_exposure = afwImage.ExposureF(exposure, deep=True)
        model = self.modelImage.run(subtracted_exposure, catalog, flux_key)

        for source in catalog:
            with self.modelImage.replaced_source(subtracted_exposure,
                                                 source, flux_key):

                # This parameter is the radius of the spanset
                # Changing the radius doesn't seem to affect SDSS Centroid?
                spanSet = afwGeom.SpanSet.fromShape(3)
                spanSet = spanSet.shiftedBy(Extent2I(source.getCentroid()))
                footprint = afwDetection.Footprint(spanSet)
                peak = footprint.peaks.addNew()
                peak.setFx(source.getCentroid().getX())
                peak.setFy(source.getCentroid().getY())

                # Check for NaNs
                if(source.getCentroid().getX() == source.getCentroid().getX()):
                    peak.setIx(int(source.getCentroid().getX()))
                    peak.setIy(int(source.getCentroid().getY()))
                source.setFootprint(footprint)

                try:
                    self.sdssCentroid.measure(source, subtracted_exposure)
                except (MeasurementError):
                    pass







