

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

    peak_significance_cutoff = pexConfig.Field(
        dtype=float,
        default=0.4,
        doc="Minimum signifiance ratio between new peaks and the existing model",
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

        result = self.run(exposure)
        if(result is None or result.source_catalog is None):
            return

        sensorRef.put(result.source_catalog, "crowdedsrc")
        sensorRef.put(exposure, "subtractedimg")
        sensorRef.put(result.model_image, "modelimg")

    @pipeBase.timeMethod
    def run(self, exposure):

        source_catalog = afwTable.SourceCatalog(self.schema)

        for detection_round in range(1, self.config.num_iterations + 1):

            detection_catalog = afwTable.SourceCatalog(self.schema)
            if(len(source_catalog) > 0):
                residual_exposure = afwImage.ExposureF(exposure, deep=True)

                # This subtracts the model from its input in place.
                # Needs to have a better name than just run.
                model_image = self.modelImageTask.run(residual_exposure,
                                                      source_catalog,
                                                      self.simultaneousPsfFlux_key)

                model_convolution = self.detection.convolveImage(model_image,
                                                                        exposure.getPsf(),
                                                                        doSmooth=True)
                # .middle refers to only the inside area of the convolved image,
                # away from edges, having meaningful data.
                model_significance_image = model_convolution.middle

            else:
                residual_exposure = exposure
                model_image = None
                model_significance_image = None

            detRes = self.detection.run(detection_catalog, residual_exposure)

            for source in detRes.sources:
                for peak in source.getFootprint().peaks:

                    # If this is round >=2, short circuit on insignificant
                    # peaks.
                    if model_image is not None:
                        # Add a 1e-6 epsilon to prevent divide-by-zero
                        model_sig = model_significance_image.getImage()[peak.getF()] + 1e-6
                        value_ratio = peak.getPeakValue()/model_sig
                        if(value_ratio < self.config.peak_significance_cutoff):
                            continue
                    child = source_catalog.addNew()
                    child['coarse_centroid_x'] = peak.getFx()
                    child['coarse_centroid_y'] = peak.getFy()

            self.log.info("Source catalog length after detection round %d: %d",
                          detection_round, len(source_catalog))

            solver_matrix = CrowdedFieldMatrix(exposure, source_catalog,
                                               self.simultaneousPsfFlux_key)
            status = solver_matrix.solve()
            if(status != solver_matrix.SUCCESS):
                self.log.error(f"Matrix solution failed on iteration {detection_round}")
                return None

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
            status = solver_matrix.solve()
            if(status != solver_matrix.SUCCESS):
                self.log.error(f"Matrix solution failed on iteration {detection_round}")
                return None

        self.log.info("Final source catalog length: %d", len(source_catalog))

        # Subtract in-place
        model_image = self.modelImageTask.run(exposure, source_catalog,
                                              self.simultaneousPsfFlux_key)

        return pipeBase.Struct(source_catalog=source_catalog,
                               model_image=model_image)


