
#include "lsst/base.h"
#include "lsst/afw/table/Exposure.h"
#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/pipe/crowd/CrowdedFieldMatrix.h"
#include "lsst/log/Log.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/Point.h"

using namespace lsst::afw;

namespace lsst {
namespace pipe {
namespace crowd {


LOG_LOGGER _log = LOG_GET("pipe.crowd");

template <typename PixelT>
CrowdedFieldMatrix<PixelT>::CrowdedFieldMatrix(CONST_PTR(afw::image::Exposure<PixelT>) exposure) :
            _exposure(exposure),
            _nStars(0) { };
    
template <typename PixelT>
void CrowdedFieldMatrix<PixelT>::addSource(double x, double y) {
    std::shared_ptr<detection::Psf::Image> psfImage;

    psfImage = _exposure->getPsf()->computeImage(geom::Point2D(x, y));
    LOGL_INFO(_log, "PSF image size %i, %i", psfImage->getWidth(), psfImage->getHeight());

    for (int y = 0; y != psfImage->getHeight(); ++y) {
        for (int x = 0; x != psfImage->getWidth(); ++x) {
            PixelT psfValue = psfImage->get(geom::Point2I(x,y), image::LOCAL);
            int pixelIndex = _exposure->getMaskedImage().getHeight() * psfImage->indexToPosition(x, image::X) + psfImage->indexToPosition(y, image::Y);
            _matrixEntries.push_back(Eigen::Triplet<PixelT>(pixelIndex, _nStars, psfValue));
        }
    }
    _nStars++;
            
}

template <typename PixelT>
void CrowdedFieldMatrix<PixelT>::addSources(ndarray::Array<double const, 1>  &x,
                               ndarray::Array<double const, 1>  &y) {
    LOGL_INFO(_log, "array length %i", x.getSize<0>());
    if(x.getSize<0>() != y.getSize<0>()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthError, "x and y must be the same length.");
    }

    for(int n = 0; n < x.getSize<0>(); n++) {
        addSource(x[n], y[n]);
    }
}

template <typename PixelT>
std::list<std::tuple<int, int, PixelT>> CrowdedFieldMatrix<PixelT>::getMatrixEntries() {
    std::list<std::tuple<int, int, PixelT>> output;
    for (auto ptr = _matrixEntries.begin(); ptr < _matrixEntries.end(); ptr++) {
        output.push_back(std::tuple<int, int, PixelT>(ptr->col(), ptr->row(), ptr->value()));
    }
    return output;
}

template class CrowdedFieldMatrix<float>;


} // namespace crowd
} // namespace pipe
} // namespace lsst
