
import numpy as np
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
from lsst.utils.timer import timeMethod

from contextlib import contextmanager

class ModelImageTaskConfig(pexConfig.Config):
    """Config for ModelImageTaskConfig"""

    pass

class ModelImageTask(pipeBase.Task):
    ConfigClass = ModelImageTaskConfig
    _DefaultName = "modelImageTask"


    @timeMethod
    def run(self, exposure, catalog, catalog_key):

        model_image = afwImage.MaskedImageF(exposure.getMaskedImage(),
                                           deep=True)
        model_arr = model_image.image.getArray()
        model_arr[:] = 0.0

        for source in catalog:
            centroid = source.getCentroid()
            psf_image = exposure.getPsf().computeImage(centroid)
            psf_image *= source[catalog_key]
            bbox = psf_image.getBBox()
            bbox.clip(model_image.getBBox())
            image_subregion = afwImage.ImageF(model_image.getImage(),
                                    bbox, afwImage.LOCAL)
            image_subregion += psf_image[bbox].convertF()

        original_image = exposure.getMaskedImage()
        original_image -= model_image
        return model_image

    @timeMethod
    def makeModelSubtractedImage(self, exposure, catalog, catalog_key):

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

    @contextmanager
    def replaced_source(self, exposure, source, flux_key):
        '''Context manager to take a source-subtracted exposure
        and re-insert one source of interest, then re-remove the source
        when finished.
        '''
        subtracted_image = exposure.getMaskedImage()
        centroid = source.getCentroid()
        psf_image = exposure.getPsf().computeImage(centroid)
        psf_image *= source[flux_key]
        bbox = psf_image.getBBox()
        bbox.clip(subtracted_image.getBBox())

        image_subregion = afwImage.ImageF(subtracted_image.getImage(),
                                          bbox, afwImage.LOCAL)
        image_subregion += psf_image[bbox].convertF()

        try:
            yield subtracted_image
        finally:
            image_subregion -= psf_image[bbox].convertF()

