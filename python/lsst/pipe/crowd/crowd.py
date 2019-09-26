

import numpy as np
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable

from .crowdedFieldMatrix import CrowdedFieldMatrix

class CrowdedFieldTaskConfig(pexConfig.Config):
    """Config for CrowdedFieldTask"""

    pass


class CrowdedFieldTask(pipeBase.CmdLineTask):
    ConfigClass = CrowdedFieldTaskConfig
    RunnerClass = pipeBase.TaskRunner
    _DefaultName = "crowdedFieldTask"

    def _getConfigName(self):
        return None

    def _getMetadataName(self):
        return None

    def __init__(self, **kwargs):
        pipeBase.CmdLineTask.__init__(self, **kwargs)

    @pipeBase.timeMethod
    def runDataRef(self, sensorRef):
        exposure = sensorRef.get("calexp")

        sources = sensorRef.get("src")

        mapper = afwTable.SchemaMapper(sources.schema)
        mapper.addMinimalSchema(sources.schema, True)
        simultaneousPsfFlux_key = mapper.editOutputSchema().addField(
            "crowd_psfFlux_flux", type=np.float32,
            doc="PSF Flux from simultaneous fitting")
        mapper.editOutputSchema().setAliasMap(sources.schema.getAliasMap())

        output_catalog = afwTable.SourceCatalog(mapper.getOutputSchema())
        output_catalog.extend(sources, mapper=mapper)

        solver_matrix = CrowdedFieldMatrix(exposure, output_catalog,
                                           simultaneousPsfFlux_key)
        solver_matrix.solve()

        sensorRef.put(output_catalog, "crowdedsrc")

