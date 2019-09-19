
#include "lsst/base.h"
#include "lsst/afw/table/Exposure.h"
#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/pipe/crowd/CrowdedFieldMatrix.h"
#include "lsst/log/Log.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/Point.h"

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using namespace lsst::afw;

namespace lsst {
namespace pipe {
namespace crowd {


LOG_LOGGER _log = LOG_GET("pipe.crowd");

template <typename PixelT>
CrowdedFieldMatrix<PixelT>::CrowdedFieldMatrix(const afw::image::Exposure<PixelT>&  exposure) :
            _exposure(exposure),
            _nStars(0),
            _inputsFrozen(false),
            _nRows(0),
            _nColumns(0) { };
    
template <typename PixelT>
void CrowdedFieldMatrix<PixelT>::addSource(double x, double y) {
    std::shared_ptr<detection::Psf::Image> psfImage;

    if(_inputsFrozen) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, "Cannot add sources after calling solve(); inputs frozen");
    }

    psfImage = _exposure.getPsf()->computeImage(geom::Point2D(x, y));
    LOGL_INFO(_log, "PSF image size %i, %i", psfImage->getWidth(), psfImage->getHeight());

    for (int y = 0; y != psfImage->getHeight(); ++y) {
        for (int x = 0; x != psfImage->getWidth(); ++x) {
            PixelT psfValue = psfImage->get(geom::Point2I(x,y), image::LOCAL);
            int pixelIndex = _exposure.getMaskedImage().getHeight() * psfImage->indexToPosition(x, image::X) + psfImage->indexToPosition(y, image::Y);
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

    for(size_t n = 0; n < x.getSize<0>(); n++) {
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

    _inputsFrozen = true;

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
    /*
     * This part could be broken out into a separate function, or extracted
     * from the previous loop with some thought. 
     */
    int max_column = 0;
    int max_row = 0;
    for(auto ptr = _matrixEntries.begin(); ptr < _matrixEntries.end(); ptr++) {
        if(ptr->col() > max_column) { max_column = ptr->col(); };
        if(ptr->row() > max_row) { max_row = ptr->row(); };
    }
    _nRows = max_row + 1;
    _nColumns = max_column + 1;

    return outputMap;
}



template <typename PixelT>
const Eigen::Matrix<PixelT, Eigen::Dynamic, 1> CrowdedFieldMatrix<PixelT>::makeDataMatrix() {
    
    Eigen::Matrix<PixelT, Eigen::Dynamic, 1> dataMatrix(1, _nColumns);
    for(auto ptr = _pixelMapping.begin(); ptr != _pixelMapping.end(); ptr++) {
        int pixelId = ptr->first;
        int columnId = ptr->second;

        int x = floor(pixelId / _exposure.getMaskedImage().getHeight());
        int y = pixelId % _exposure.getMaskedImage().getHeight();

        dataMatrix(columnId, 1) = _exposure.getMaskedImage().getImage()->get(geom::Point2I(x,y), image::PARENT);
    }
    return dataMatrix;
}

template <typename PixelT>
void CrowdedFieldMatrix<PixelT>::solve() {


    Eigen::SparseMatrix<PixelT> paramMatrix;
    Eigen::SparseMatrix<PixelT> dataMatrix;
    std::map<int, int> pixelMapping;

    _inputsFrozen = true;

    _pixelMapping = renameMatrixRows();

    LOGL_INFO(_log, "parameter matrix size %i rows, %i cols", _nRows, _nColumns);

    paramMatrix = Eigen::SparseMatrix<PixelT>(_nRows, _nColumns);
    paramMatrix.setFromTriplets(_matrixEntries.begin(), _matrixEntries.end());

}

template class CrowdedFieldMatrix<float>;


} // namespace crowd
} // namespace pipe
} // namespace lsst
