

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class CrowdedFieldTaskConfig(pexConfig.Config):
    """Config for CrowdedFieldTask"""

    pass


class CrowdedFieldTask(pipeBase.CmdLineTask):
    ConfigClass = CrowdedFieldTaskConfig
    RunnerClass = pipeBase.ButlerInitializedTaskRunner
    _DefaultName = "crowdedFieldTask"

    def __init__(self, **kwargs):
        pipeBase.CmdLineTask.__init__(self, **kwargs)

    @pipeBase.timeMethod
    def runDataRef(self, sensorRef):
        exposure = sensorRef.get("calexp")

        sources = sensorRef.get("src")

        self.log.info("end of runDataRef")

