

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

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


        solver_matrix = CrowdedFieldMatrix(exposure.getMaskedImage(),
                                           exposure.getPsf())
        solver_matrix.addSources(sources['x'], sources['y'])

        self.log.info("end of runDataRef")

