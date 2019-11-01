

import numpy as np
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
from lsst.meas.algorithms import SourceDetectionTask, SourceDetectionConfig
from lsst.pipe.base import ArgumentParser

from .crowdedFieldMatrix import CrowdedFieldMatrix
from .modelImage import ModelImageTask, ModelImageTaskConfig
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

    modelImageTask = pexConfig.ConfigurableField(
            target=ModelImageTask,
            doc="Task for manipulating model images"
    )

    num_iterations = pexConfig.Field(
        dtype=int,
        default=2,
        doc="Number of detect-measure-subtract iterations",
    )

    fitSimultaneousPositions = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Include the source positions when fitting for fluxes?",
    )

    def validate(self):
        super().validate()
        if(self.fitSimultaneousPositions):
           raise ValueError("fitSimultaneousPositions not currently supported.")


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
            "crowd_psfFlux_flux_instFlux", type=np.float64,
            doc="PSF Flux from simultaneous fitting")
        self.schema.getAliasMap().set("slot_PsfFlux",
                                      "crowd_psfFlux_flux")
        self.centroid_key = afwTable.Point2DKey.addFields(self.schema,
                                                          "coarse_centroid",
                                                          "Detection peak", "pixels")
        self.schema.getAliasMap().set("slot_Centroid",
                                      "coarse_centroid")

        self.makeSubtask("detection", schema=self.schema)
        self.makeSubtask("centroid", schema=self.schema)
        self.makeSubtask("modelImageTask")

    @pipeBase.timeMethod
    def runDataRef(self, sensorRef):

        exposure = sensorRef.get("calexp")

        source_catalog = self.run(exposure)

        sensorRef.put(source_catalog, "crowdedsrc")
        sensorRef.put(exposure, "subtractedimg")

    @pipeBase.timeMethod
    def run(self, exposure):

        source_catalog = afwTable.SourceCatalog(self.schema)

        for detection_round in range(1, self.config.num_iterations + 1):

            detection_catalog = afwTable.SourceCatalog(self.schema)
            if(len(source_catalog) > 0):
                residual_exposure = afwImage.ExposureF(exposure, deep=True)

                self.modelImageTask.makeModelSubtractedImage(residual_exposure,
                                                             source_catalog,
                                                             self.simultaneousPsfFlux_key)
            else:
                residual_exposure = exposure


            detRes = self.detection.run(detection_catalog, residual_exposure)

            for source in detRes.sources:
                for peak in source.getFootprint().peaks:
                    child = source_catalog.addNew()
                    child['coarse_centroid_x'] = peak.getFx()
                    child['coarse_centroid_y'] = peak.getFy()

            self.log.info("Source catalog length after detection round %d: %d",
                          detection_round, len(source_catalog))

            solver_matrix = CrowdedFieldMatrix(exposure, source_catalog,
                                               self.simultaneousPsfFlux_key)
            solver_matrix.solve()

            source_catalog.schema.getAliasMap().set("slot_Centroid",
                                                    "coarse_centroid")

            self.centroid.run(exposure, source_catalog,
                         self.simultaneousPsfFlux_key)

            # Move the centroid slot from the coarse peak values to the
            # SdssCentroid values.
            source_catalog.schema.getAliasMap().set("slot_Centroid",
                                                    "centroid")

            # Now that we have more precise centroids, re-fit the fluxes
            solver_matrix = CrowdedFieldMatrix(exposure, source_catalog,
                                               self.simultaneousPsfFlux_key)
            solver_matrix.solve()

        self.log.info("Final source catalog length: %d", len(source_catalog))

        # Subtract in-place
        self.modelImageTask.makeModelSubtractedImage(exposure,
                                                     source_catalog,
                                                     self.simultaneousPsfFlux_key)

        return source_catalog


