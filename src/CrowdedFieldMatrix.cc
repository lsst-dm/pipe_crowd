
#include "lsst/base.h"
#include "lsst/afw/table/Exposure.h"
#include "lsst/afw/table/Key.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/log/Log.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/image/Mask.h"

#include "lsst/pipe/crowd/CrowdedFieldMatrix.h"

#include "Eigen/SparseCore"
#include "Eigen/IterativeLinearSolvers"

using namespace lsst;

namespace lsst {
namespace pipe {
namespace crowd {


LOG_LOGGER _log = LOG_GET("lsst.pipe.crowd.CrowdedFieldMatrix");

template <typename PixelT>
CrowdedFieldMatrix<PixelT>::CrowdedFieldMatrix(const afw::image::Exposure<PixelT> &exposure,
                                               ndarray::Array<double const, 1> &x,
                                               ndarray::Array<double const, 1> &y) :
            _exposure(exposure),
            _catalog(NULL),
            _fitCentroids(false),
            _centroidKey(afw::table::PointKey<double>()),
            _paramTracker(ParameterTracker(1)),
            _iterations(0),
            _maxIterations(500)
{
    _matrixEntries = _makeMatrixEntries(exposure, x, y);
    _dataVector = makeDataVector();
};

template <typename PixelT>
CrowdedFieldMatrix<PixelT>::CrowdedFieldMatrix(const afw::image::Exposure<PixelT> &exposure,
                                               afw::table::SourceCatalog *catalog,
                                               afw::table::Key<double> fluxKey,
                                               bool fitCentroids,
                                               afw::table::PointKey<double> centroidKey) :
            _exposure(exposure),
            _catalog(catalog),
            _fluxKey(fluxKey),
            _fitCentroids(fitCentroids),
            _centroidKey(centroidKey),
            _paramTracker(ParameterTracker(fitCentroids ? 3 : 1)),
            _iterations(0),
            _maxIterations(500)
{
    _matrixEntries = _makeMatrixEntries(exposure, catalog);
    _dataVector = makeDataVector();
};


template <typename PixelT>
std::vector<Eigen::Triplet<PixelT>> CrowdedFieldMatrix<PixelT>::_makeMatrixEntries(const afw::image::Exposure<PixelT> &exposure,
                                               ndarray::Array<double const, 1> &x,
                                               ndarray::Array<double const, 1> &y) {

    std::vector<Eigen::Triplet<PixelT>> matrixEntries;

    if(x.getSize<0>() != y.getSize<0>()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthError, "x and y must be the same length.");
    }

    for(size_t n = 0; n < x.getSize<0>(); ++n) {
        _addSource(exposure, matrixEntries, n, x[n], y[n]);
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
    geom::Point2D centroid;
    size_t n = 0;
    for(auto rec = catalog->begin(); rec < catalog->end(); ++rec, ++n) {
        PixelT estFlux = PixelT();
        centroid = rec->getCentroid();
        if(_fitCentroids) {
            estFlux = rec->getPsfInstFlux();
        }
        _addSource(exposure, matrixEntries, n, centroid[0], centroid[1], estFlux);
    }
    return matrixEntries;
}

template <typename PixelT>
void CrowdedFieldMatrix<PixelT>::_addSource(const afw::image::Exposure<PixelT> &exposure,
                                            std::vector<Eigen::Triplet<PixelT>> &matrixEntries,
                                            int nStar, double x, double y, PixelT estFlux) {
    using afw::image::Mask;
    using afw::image::MaskPixel;
    std::shared_ptr<afw::detection::Psf::Image> psfImage, psfImage_dx, psfImage_dy;
    int pixelShift_dx, pixelShift_dy;
    MaskPixel maskValue;
    Mask<MaskPixel> psfShapedMask;
    MaskPixel maskFlagsForRejection = Mask<MaskPixel>::getPlaneBitMask({"SAT", "BAD", "EDGE", "CR"});
    geom::Box2I clippedBBox;

    psfImage = exposure.getPsf()->computeImage(geom::Point2D(x, y));
    clippedBBox = psfImage->getBBox();
    clippedBBox.clip(exposure.getMaskedImage().getBBox());
    psfShapedMask = Mask<MaskPixel>(*exposure.getMaskedImage().getMask(), clippedBBox);

    float pixelNudge = 1.0;

    if(_fitCentroids) {
        psfImage_dx = exposure.getPsf()->computeImage(geom::Point2D(x + pixelNudge, y));
        psfImage_dy = exposure.getPsf()->computeImage(geom::Point2D(x, y + pixelNudge));

        // Assume that the XY0 only changes in the direction of the nudge
        pixelShift_dx = psfImage_dx->getX0() - psfImage->getX0();
        pixelShift_dy = psfImage_dy->getY0() - psfImage->getY0();
    }


    _paramTracker.addSource(nStar);

    for (int y = 0; y != psfImage->getHeight(); ++y) {
        for (int x = 0; x != psfImage->getWidth(); ++x) {

            maskValue = psfShapedMask.get(geom::Point2I(x,y), afw::image::LOCAL);
            if((maskValue & maskFlagsForRejection) > 0) {
                continue;
            }

            PixelT psfValue = psfImage->get(geom::Point2I(x, y), afw::image::LOCAL);
            int pixelIndex = _paramTracker.makePixelId(psfImage->indexToPosition(x, afw::image::X),
                                                       psfImage->indexToPosition(y, afw::image::Y));
            int paramIndex = _paramTracker.getSourceParameterId(nStar, 0);
            matrixEntries.push_back(Eigen::Triplet<PixelT>(pixelIndex, paramIndex, psfValue));

            if(_fitCentroids && (x + pixelShift_dx >= 0) && (x + pixelShift_dx < psfImage->getWidth())) {
                PixelT psfValue_dx = psfImage->get(geom::Point2I(x + pixelShift_dx, y), afw::image::LOCAL);
                PixelT deriv_x = estFlux * (psfValue - psfValue_dx)/pixelNudge;

                int paramIndex = _paramTracker.getSourceParameterId(nStar, 1);
                matrixEntries.push_back(Eigen::Triplet<PixelT>(pixelIndex, paramIndex, deriv_x));
            }

            if(_fitCentroids && (y + pixelShift_dy >= 0) && (y + pixelShift_dy < psfImage->getHeight())) {
                PixelT psfValue_dy = psfImage->get(geom::Point2I(x, y + pixelShift_dy), afw::image::LOCAL);
                PixelT deriv_y = estFlux * (psfValue - psfValue_dy)/pixelNudge;

                int paramIndex = _paramTracker.getSourceParameterId(nStar, 2);
                matrixEntries.push_back(Eigen::Triplet<PixelT>(pixelIndex, paramIndex, deriv_y));
            }

        }
    }
}

template <typename PixelT>
const std::list<std::tuple<int, int, PixelT>> CrowdedFieldMatrix<PixelT>::getMatrixEntries() {
    std::list<std::tuple<int, int, PixelT>> output;
    for (auto ptr = _matrixEntries.begin(); ptr < _matrixEntries.end(); ++ptr) {
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
SolverStatus CrowdedFieldMatrix<PixelT>::solve() {


    Eigen::SparseMatrix<PixelT> paramMatrix;

    LOGL_INFO(_log, "parameter matrix size %i rows, %i cols",
              _paramTracker.nRows(), _paramTracker.nColumns());

    paramMatrix = Eigen::SparseMatrix<PixelT>(_paramTracker.nRows(),
                                              _paramTracker.nColumns());
    paramMatrix.setFromTriplets(_matrixEntries.begin(), _matrixEntries.end());

    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<PixelT>> lscg;
    lscg.setTolerance(1e-6);
    lscg.setMaxIterations(_maxIterations);
    lscg.compute(paramMatrix);
    _result = lscg.solve(_dataVector);


    if(_catalog) {
        size_t n = 0;
        for(auto rec = _catalog->begin(); rec < _catalog->end(); ++rec, ++n) {
            rec->set(_fluxKey, _result(_paramTracker.getSourceParameterId(n, 0), 0));
            if(_fitCentroids && _centroidKey.isValid()) {
                auto deltaCentroid = geom::Extent2D(_result(_paramTracker.getSourceParameterId(n, 1), 0),
                                                    _result(_paramTracker.getSourceParameterId(n, 2), 0));
                rec->set(_centroidKey, rec->getCentroid() + deltaCentroid);
            }
        }
    }

    _iterations = lscg.iterations();
    if(_iterations == _maxIterations) {
        LOGL_WARN(_log, "eigen failed to solve in %i iterations", _iterations);
        return SolverStatus::FAILURE;
    } else {
        LOGL_INFO(_log, "eigen solved in %i iterations", _iterations);
        return SolverStatus::SUCCESS;
    }

}

template <typename PixelT>
Eigen::Matrix<PixelT, Eigen::Dynamic, 1> CrowdedFieldMatrix<PixelT>::result() {
    return _result;
}

template <typename PixelT>
const std::map<std::tuple<int, int>, int> CrowdedFieldMatrix<PixelT>::getParameterMapping() {
    return _paramTracker.getParameterMapping();
}

template <typename PixelT>
const std::map<std::tuple<int, int>, int> CrowdedFieldMatrix<PixelT>::getPixelMapping() {
    return _paramTracker.getPixelMapping();
}

template class CrowdedFieldMatrix<float>;


} // namespace crowd
} // namespace pipe
} // namespace lsst
