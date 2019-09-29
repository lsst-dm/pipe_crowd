
import numpy as np
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage

class CatalogPsfSubtractTaskConfig(pexConfig.Config):
    """Config for CatalogPsfSubtractTask"""

    pass

class CatalogPsfSubtractTask(pipeBase.Task):
    ConfigClass = CatalogPsfSubtractTaskConfig
    _DefaultName = "catalogPsfSubtractTask"


    @pipeBase.timeMethod
    def run(self, exposure, catalog, catalog_key):

        subtracted_image = afwImage.MaskedImageF(exposure.getMaskedImage(),
                                                deep=False)

        for source in catalog:
            centroid = source.getCentroid()
            psf_image = exposure.getPsf().computeImage(centroid)
            psf_image *= source[catalog_key]
            bbox = psf_image.getBBox()
            bbox.clip(subtracted_image.getBBox())
            image_subregion = afwImage.ImageF(subtracted_image.getImage(),
                                    bbox, afwImage.LOCAL)
            image_subregion -= psf_image[bbox].convertF()

        return subtracted_image


