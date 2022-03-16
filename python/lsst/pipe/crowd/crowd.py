

import numpy as np
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
from lsst.meas.algorithms import SourceDetectionTask, SourceDetectionConfig
from lsst.pipe.base import ArgumentParser
import lsst.pipe.base.connectionTypes as cT
from lsst.utils.timer import timeMethod

from scipy.spatial import cKDTree


from .crowdedFieldMatrix import CrowdedFieldMatrix
from .modelImage import ModelImageTask, ModelImageTaskConfig
from .centroid import CrowdedCentroidTask, CrowdedCentroidTaskConfig


class CrowdedFieldConnections(pipeBase.PipelineTaskConnections, dimensions=("instrument", "visit",
                                                                            "detector"),
                              defaultTemplates={}):

    calexp = cT.Input(
        doc="Input image to measure.",
        name="calexp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
    )

    crowdedFieldCat = cT.Output(
        doc="Output catalog from simultaneous fitting.",
        name="pipeCrowd_src",
        storageClass="SourceCatalog",
        dimensions=("instrument", "visit", "detector")
    )

    crowdedFieldModel = cT.Output(
        doc="Model image from sources measured in simultaneous fitting.",
        name="pipeCrowd_model",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector")
    )

    crowdedFieldResidual = cT.Output(
        doc="Residual image after subtracting sources.",
        name="pipeCrowd_residual",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector")
    )


class CrowdedFieldTaskConfig(pipeBase.PipelineTaskConfig, pipelineConnections=CrowdedFieldConnections):
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

    minCentroidSeparation = pexConfig.Field(
        dtype=float,
        default=1.0,
        doc="Filter sources before fitting so that none have separation less than minCentroidSeparation",
    )

    def validate(self):
        super().validate()
        if(self.fitSimultaneousPositions):
           raise ValueError("fitSimultaneousPositions not currently supported.")



class CrowdedFieldTask(pipeBase.PipelineTask):
    ConfigClass = CrowdedFieldTaskConfig
    # RunnerClass = pipeBase.TaskRunner
    _DefaultName = "crowdedFieldTask"

    def setDefaults(self):
        super().setDefaults()
        self.detection.thresholdPolarity = "positive"

    def __init__(self, **kwargs):
        pipeBase.PipelineTask.__init__(self, **kwargs)
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

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        outputs = self.run(inputs['calexp'])
        if(outputs is not None):
            butlerQC.put(outputs, outputRefs)

    @timeMethod
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
                self.log.error(f"Matrix solution failed on iteration {detection_round} solve 1")
                return None

            source_catalog.schema.getAliasMap().set("slot_Centroid",
                                                    "coarse_centroid")

            self.centroid.run(exposure, source_catalog,
                         self.simultaneousPsfFlux_key)

            # Move the centroid slot from the coarse peak values to the
            # SdssCentroid values.
            source_catalog.schema.getAliasMap().set("slot_Centroid",
                                                    "centroid")

            # Delete sources with centroids that are too close together.
            # This is pretty ad hoc.
            centroid_tree = cKDTree(np.stack([source_catalog['centroid_x'], source_catalog['centroid_y']], axis=1))

            # TODO: better handle the case of three+ sources inside the matching radius.
            records_to_delete = set(j for (i,j) in centroid_tree.query_pairs(self.config.minCentroidSeparation))
            if(len(records_to_delete) == 0):
                self.log.info("No overlapping sources to delete.")
            else:
                self.log.info("Deleting overlaping sources: " +  ", ".join(f"{x:d}" for x in records_to_delete))

                for n in sorted(list(records_to_delete), reverse=True):
                    del source_catalog[n]

                source_catalog = source_catalog.copy(deep=True)


            double_check_for_pairs = True
            if double_check_for_pairs:
                centroid_tree = cKDTree(np.stack([source_catalog['centroid_x'], source_catalog['centroid_y']], axis=1))
                remaining_pairs = set(j for (i,j) in centroid_tree.query_pairs(self.config.minCentroidSeparation))
                if(len(remaining_pairs) > 0):
                    self.log.warn("Close-pairs remain after filtering: " +  ", ".join(f"{x:d}" for x in remaining_pairs))

            # Now that we have more precise centroids, re-fit the fluxes
            solver_matrix = CrowdedFieldMatrix(exposure, source_catalog,
                                               self.simultaneousPsfFlux_key)

            status = solver_matrix.solve()
            if(status != solver_matrix.SUCCESS):
                self.log.error(f"Matrix solution failed on iteration {detection_round} solve 2")
                raise RuntimeError(f"Matrix solution failed on iteration {detection_round} solve 2")
                return None

        self.log.info("Final source catalog length: %d", len(source_catalog))

        # Subtract in-place
        model_image = self.modelImageTask.run(exposure, source_catalog,
                                              self.simultaneousPsfFlux_key)

        model_exposure = afwImage.ExposureF(model_image, wcs=exposure.getWcs())

        return pipeBase.Struct(crowdedFieldCat=source_catalog,
                               crowdedFieldResidual=exposure,
                               crowdedFieldModel=model_exposure)


