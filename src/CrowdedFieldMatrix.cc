
#include "lsst/base.h"
#include "lsst/afw/table/Exposure.h"
#include "lsst/afw/table/Key.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/log/Log.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/image/Mask.h"

#include "lsst/pipe/crowd/CrowdedFieldMatrix.h"

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using namespace lsst::afw;

namespace lsst {
namespace pipe {
namespace crowd {


LOG_LOGGER _log = LOG_GET("pipe.crowd");

template <typename PixelT>
CrowdedFieldMatrix<PixelT>::CrowdedFieldMatrix(const afw::image::Exposure<PixelT> &exposure,
                                               ndarray::Array<double const, 1> &x,
                                               ndarray::Array<double const, 1> &y) :
            _exposure(exposure),
            _catalog(NULL),
            _paramTracker(ParameterTracker(1))
{
    _matrixEntries = _makeMatrixEntries(exposure, x, y);
    _dataVector = makeDataVector();
};

template <typename PixelT>
CrowdedFieldMatrix<PixelT>::CrowdedFieldMatrix(const afw::image::Exposure<PixelT> &exposure,
                                               afw::table::SourceCatalog *catalog,
                                               afw::table::Key<float> fluxKey) :
            _exposure(exposure),
            _catalog(catalog),
            _fluxKey(fluxKey),
            _paramTracker(ParameterTracker(1))
{
    _matrixEntries = _makeMatrixEntries(exposure, catalog);
    _dataVector = makeDataVector();
};


template <typename PixelT>
std::vector<Eigen::Triplet<PixelT>> CrowdedFieldMatrix<PixelT>::_makeMatrixEntries(const afw::image::Exposure<PixelT> &exposure,
                                               ndarray::Array<double const, 1> &x,
                                               ndarray::Array<double const, 1> &y) {

    if(x.getSize<0>() != y.getSize<0>()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthError, "x and y must be the same length.");
    }

    std::vector<Eigen::Triplet<PixelT>> matrixEntries;
    for(size_t n = 0; n < x.getSize<0>(); n++) {
        _addSource(exposure, matrixEntries, (int) n, x[n], y[n]);
    }
    return matrixEntries;
}

template <typename PixelT>
std::vector<Eigen::Triplet<PixelT>> CrowdedFieldMatrix<PixelT>::_makeMatrixEntries(const afw::image::Exposure<PixelT> &exposure,
                                            afw::table::SourceCatalog *catalog) {
    if(catalog == NULL) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, "sourceCatalog is NULL");
    }
    std::vector<Eigen::Triplet<PixelT>> matrixEntries;
    size_t n = 0;
    afw::geom::Point2D centroid;
    for(auto rec = catalog->begin(); rec < catalog->end(); rec++, n++) {
        centroid = rec->getCentroid();
        _addSource(exposure, matrixEntries, (int) n, centroid[0], centroid[1]);
    }
    return matrixEntries;
}

template <typename PixelT>
void CrowdedFieldMatrix<PixelT>::_addSource(const afw::image::Exposure<PixelT> &exposure,
                                            std::vector<Eigen::Triplet<PixelT>> &matrixEntries,
                                            int nStar, double x, double y) {
    std::shared_ptr<detection::Psf::Image> psfImage;
    afw::image::MaskPixel maskValue;
    afw::image::Mask<afw::image::MaskPixel> psfShapedMask;
    afw::image::MaskPixel maskFlagsForRejection = afw::image::Mask<afw::image::MaskPixel>::getPlaneBitMask({"SAT", "BAD", "EDGE", "CR"});
    geom::Box2I clippedBBox;

    psfImage = exposure.getPsf()->computeImage(geom::Point2D(x, y));
    clippedBBox = psfImage->getBBox();
    clippedBBox.clip(exposure.getMaskedImage().getBBox());
    psfShapedMask = afw::image::Mask<afw::image::MaskPixel>(*exposure.getMaskedImage().getMask(), clippedBBox);

    _paramTracker.addSource(nStar);

    for (int y = 0; y != psfImage->getHeight(); ++y) {
        for (int x = 0; x != psfImage->getWidth(); ++x) {

            maskValue = psfShapedMask.get(geom::Point2I(x,y), image::LOCAL);
            if((maskValue & maskFlagsForRejection) > 0) {
                continue;
            }

            PixelT psfValue = psfImage->get(geom::Point2I(x,y), image::LOCAL);
            int pixelIndex = _paramTracker.makePixelId(psfImage->indexToPosition(x, image::X),
                                                       psfImage->indexToPosition(y, image::Y));
            int paramIndex = _paramTracker.getSourceParameterId(nStar, 0);
            matrixEntries.push_back(Eigen::Triplet<PixelT>(pixelIndex, paramIndex, psfValue));
        }
    }
}

template <typename PixelT>
const std::list<std::tuple<int, int, PixelT>> CrowdedFieldMatrix<PixelT>::getMatrixEntries() {
    std::list<std::tuple<int, int, PixelT>> output;
    for (auto ptr = _matrixEntries.begin(); ptr < _matrixEntries.end(); ptr++) {
        output.push_back(std::tuple<int, int, PixelT>(ptr->col(), ptr->row(), ptr->value()));
    }
    return output;
}

template <typename PixelT>
const Eigen::Matrix<PixelT, Eigen::Dynamic, 1> CrowdedFieldMatrix<PixelT>::getDataVector() {
    return _dataVector;
}

template <typename PixelT>
const Eigen::Matrix<PixelT, Eigen::Dynamic, 1> CrowdedFieldMatrix<PixelT>::makeDataVector() {

    Eigen::Matrix<PixelT, Eigen::Dynamic, 1> dataMatrix(_paramTracker.nRows(), 1);
    auto img = _exposure.getMaskedImage().getImage();
    int * pixelId;

    for (int y = 0; y != img->getHeight(); ++y) {
        for (auto ptr = img->row_begin(y), end = img->row_end(y), x = 0; ptr != end; ++ptr, ++x) {

            pixelId = _paramTracker.getPixelId(x, y);
            if(pixelId == NULL) {
                continue;
            }
            dataMatrix(*pixelId, 0) = *ptr;
        }
    }
    return dataMatrix;
}

template <typename PixelT>
Eigen::Matrix<PixelT, Eigen::Dynamic, 1> CrowdedFieldMatrix<PixelT>::solve() {


    Eigen::SparseMatrix<PixelT> paramMatrix;
    Eigen::Matrix<PixelT, Eigen::Dynamic, 1> result;

    LOGL_INFO(_log, "parameter matrix size %i rows, %i cols",
              _paramTracker.nRows(), _paramTracker.nColumns());

    paramMatrix = Eigen::SparseMatrix<PixelT>(_paramTracker.nRows(),
                                              _paramTracker.nColumns());
    paramMatrix.setFromTriplets(_matrixEntries.begin(), _matrixEntries.end());

    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<PixelT>> lscg;
    lscg.compute(paramMatrix);
    result = lscg.solve(_dataVector);

    if(_catalog) {
        size_t n = 0;
        for(auto rec = _catalog->begin(); rec < _catalog->end(); rec++, n++) {
            rec->set(_fluxKey, result(n, 0));
        }
    }

    return result;
}

template <typename PixelT>
const std::vector<std::tuple<int, int, PixelT>> CrowdedFieldMatrix<PixelT>::getDebug() {
    return _debugXYValues;
}

template class CrowdedFieldMatrix<float>;


} // namespace crowd
} // namespace pipe
} // namespace lsst
