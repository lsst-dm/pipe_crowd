

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

from .subtraction import CatalogPsfSubtractTask

class CrowdedCentroidTaskConfig(pexConfig.Config):
    """Config for CrowdedCentroidTask"""

    subtraction = pexConfig.ConfigurableField(
            target=CatalogPsfSubtractTask,
            doc="Subtract sources from image"
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
        self.makeSubtask("subtraction")

    @pipeBase.timeMethod
    def run(self, exposure, catalog, flux_key):

        subtracted_exposure = afwImage.ExposureF(exposure, deep=True)
        self.subtraction.run(subtracted_exposure, catalog, flux_key)

        for source in catalog:
            with self.subtraction.replaced_source(subtracted_exposure,
                                                  source, flux_key):
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







