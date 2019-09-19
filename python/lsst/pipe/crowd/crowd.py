

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

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


        solver_matrix = CrowdedFieldMatrix(exposure)

        solver_matrix.addSource(200.0, 200.0)
        #solver_matrix.addSources(sources['slot_Centroid_x'][:10],
        #                         sources['slot_Centroid_y'][:10])
        print(solver_matrix.getMatrixEntries()[:10])
        out = solver_matrix.renameMatrixRows()
        print(solver_matrix.getMatrixEntries()[:10])


        # solver_matrix.solve()

        self.log.info("end of runDataRef")

