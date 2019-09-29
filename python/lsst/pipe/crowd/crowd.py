

import numpy as np
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
from lsst.meas.algorithms import SourceDetectionTask, SourceDetectionConfig
from lsst.pipe.base import ArgumentParser

from .crowdedFieldMatrix import CrowdedFieldMatrix
from .subtraction import CatalogPsfSubtractTask, CatalogPsfSubtractTaskConfig
from .centroid import CrowdedCentroidTask, CrowdedCentroidTaskConfig

class CrowdedFieldTaskConfig(pexConfig.Config):
    """Config for CrowdedFieldTask"""

    centroid = pexConfig.ConfigurableField(
        target=CrowdedCentroidTask,
        doc="Centroid sources"
    )

    detection = pexConfig.ConfigurableField(
            target=SourceDetectionTask,
            doc="Detect sources"
    )

    subtraction = pexConfig.ConfigurableField(
            target=CatalogPsfSubtractTask,
            doc="Subtract sources from image"
    )

class CrowdedFieldTask(pipeBase.CmdLineTask):
    ConfigClass = CrowdedFieldTaskConfig
    RunnerClass = pipeBase.TaskRunner
    _DefaultName = "crowdedFieldTask"

    def _getConfigName(self):
        return None

    def _getMetadataName(self):
        return None

    @classmethod
    def _makeArgumentParser(cls):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument(name="--id", datasetType="calexp",
                               help="data IDs, e.g. --id visit=12345 ccd=1,2^0,3")
        return parser

    def setDefaults(self):
        super().setDefaults()
        self.detection.thresholdPolarity = "positive"

    def __init__(self, **kwargs):
        pipeBase.CmdLineTask.__init__(self, **kwargs)
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.simultaneousPsfFlux_key = self.schema.addField(
            "crowd_psfFlux_flux", type=np.float32,
            doc="PSF Flux from simultaneous fitting")
        afwTable.Point2DKey.addFields(self.schema,
                                      "coarse_centroid",
                                      "Detection peak", "pixels")
        self.schema.getAliasMap().set("slot_Centroid",
                                      "coarse_centroid")

        self.makeSubtask("detection", schema=self.schema)
        self.makeSubtask("centroid", schema=self.schema)
        self.makeSubtask("subtraction")


    @pipeBase.timeMethod
    def runDataRef(self, sensorRef):

        exposure = sensorRef.get("calexp")

        source_catalog = self.run(exposure)

        sensorRef.put(source_catalog, "crowdedsrc")
        sensorRef.put(exposure, "subtractedimg")


    @pipeBase.timeMethod
    def run(self, exposure):

        source_catalog = afwTable.SourceCatalog(self.schema)
        detRes = self.detection.run(source_catalog, exposure)

        for source in detRes.sources:
            for peak in source.getFootprint().peaks:
                child = source_catalog.addNew()
                child['coarse_centroid_x'] = peak.getFx()
                child['coarse_centroid_y'] = peak.getFy()

        solver_matrix = CrowdedFieldMatrix(exposure, source_catalog,
                                           self.simultaneousPsfFlux_key)
        solver_matrix.solve()

        self.centroid.run(exposure, source_catalog,
                     self.simultaneousPsfFlux_key)


        source_catalog.schema.getAliasMap().set("slot_Centroid",
                                                "centroid")

        # Subtract in-place
        self.subtraction.run(exposure,
                             source_catalog,
                             self.simultaneousPsfFlux_key)

        return source_catalog


