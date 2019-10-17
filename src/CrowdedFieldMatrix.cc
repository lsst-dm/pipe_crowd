
#include "lsst/base.h"
#include "lsst/afw/table/Exposure.h"
#include "lsst/afw/table/Key.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/pipe/crowd/CrowdedFieldMatrix.h"
#include "lsst/log/Log.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/image/Mask.h"

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
            _matrixEntries(_makeMatrixEntries(exposure, x, y))
{
    _pixelMapping = renameMatrixRows();
    std::tie(_nRows, _nColumns) = _setNrowsNcols();
    _dataVector = makeDataVector();
};

template <typename PixelT>
CrowdedFieldMatrix<PixelT>::CrowdedFieldMatrix(const afw::image::Exposure<PixelT> &exposure,
                                               afw::table::SourceCatalog *catalog,
                                               afw::table::Key<float> fluxKey) :
            _exposure(exposure),
            _catalog(catalog),
            _fluxKey(fluxKey),
            _matrixEntries(_makeMatrixEntries(exposure, catalog))
{
    _pixelMapping = renameMatrixRows();
    std::tie(_nRows, _nColumns) = _setNrowsNcols();
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


    for (int y = 0; y != psfImage->getHeight(); ++y) {
        for (int x = 0; x != psfImage->getWidth(); ++x) {

            maskValue = psfShapedMask.get(geom::Point2I(x,y), image::LOCAL);
            if((maskValue & maskFlagsForRejection) > 0) {
                continue;
            }

            PixelT psfValue = psfImage->get(geom::Point2I(x,y), image::LOCAL);
            int pixelIndex = exposure.getMaskedImage().getHeight() * psfImage->indexToPosition(x, image::X) + psfImage->indexToPosition(y, image::Y);
            matrixEntries.push_back(Eigen::Triplet<PixelT>(pixelIndex, nStar, psfValue));
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

/*
 *
 * The initial entries in the matrix tuple vector uses pixel numbers, but if we
 * don't touch every pixel then there will be some rows that are empty. We
 * create a remapping dictionary and renumber the entries in the triplets.
 */
template <typename PixelT>
std::map<int, int> CrowdedFieldMatrix<PixelT>::renameMatrixRows() {
    std::map<int, int> outputMap;
    std::map<int, int>::iterator mappingEntry;
    int pixelCounter = 0;

    for(auto ptr = _matrixEntries.begin(); ptr < _matrixEntries.end(); ptr++) {
        mappingEntry = outputMap.find(ptr->row());
        if(mappingEntry != outputMap.end()) {
            *ptr = Eigen::Triplet<PixelT>(mappingEntry->second, ptr->col(), ptr->value());
        } else {
            outputMap.insert({ptr->row(), pixelCounter});
            *ptr = Eigen::Triplet<PixelT>(pixelCounter,  ptr->col(), ptr->value());
            pixelCounter++;
        }
    }

    return outputMap;
}



template <typename PixelT>
const Eigen::Matrix<PixelT, Eigen::Dynamic, 1> CrowdedFieldMatrix<PixelT>::makeDataVector() {

    Eigen::Matrix<PixelT, Eigen::Dynamic, 1> dataMatrix(_nRows, 1);

    for(auto ptr = _pixelMapping.begin(); ptr != _pixelMapping.end(); ptr++) {
        int pixelId = ptr->first;
        int rowId = ptr->second;

        int x = floor(pixelId / _exposure.getMaskedImage().getHeight());
        int y = pixelId % _exposure.getMaskedImage().getHeight();

        if(rowId >= _nRows) {
            LOGL_WARN(_log, "row %i >= _nRows %i", rowId, _nRows);
            throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, "rowId >= _nRows ");
        }
        PixelT pixelValue = _exposure.getMaskedImage().getImage()->get(geom::Point2I(x,y), image::PARENT);
        dataMatrix(rowId, 0) = pixelValue;
        _debugXYValues.push_back(std::tuple<int, int, PixelT>(x, y, pixelValue));
    }
    return dataMatrix;
}

template <typename PixelT>
Eigen::Matrix<PixelT, Eigen::Dynamic, 1> CrowdedFieldMatrix<PixelT>::solve() {


    Eigen::SparseMatrix<PixelT> paramMatrix;
    Eigen::Matrix<PixelT, Eigen::Dynamic, 1> result;

    LOGL_INFO(_log, "parameter matrix size %i rows, %i cols", _nRows, _nColumns);

    paramMatrix = Eigen::SparseMatrix<PixelT>(_nRows, _nColumns);
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
const std::map<int, int> CrowdedFieldMatrix<PixelT>::getPixelMapping() {
    return _pixelMapping;
}


template <typename PixelT>
const std::vector<std::tuple<int, int, PixelT>> CrowdedFieldMatrix<PixelT>::getDebug() {
    return _debugXYValues;
}

template <typename PixelT>
const std::tuple<int, int> CrowdedFieldMatrix<PixelT>::_setNrowsNcols() {
    /*
     * This depends on the result of renameMatrixRows(), not entirely
     * satisfying.
     */
    int max_column = 0;
    int max_row = 0;
    for(auto ptr = _matrixEntries.begin(); ptr < _matrixEntries.end(); ptr++) {
        if(ptr->col() > max_column) { max_column = ptr->col(); };
        if(ptr->row() > max_row) { max_row = ptr->row(); };
    }
    return std::make_tuple(max_row + 1, max_column + 1);
}

template class CrowdedFieldMatrix<float>;


} // namespace crowd
} // namespace pipe
} // namespace lsst
