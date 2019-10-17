
#include "lsst/base.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/table/Key.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/algorithms/ImagePsf.h"
#include <Eigen/Sparse>
#include <vector>

namespace lsst {
namespace pipe {
namespace crowd {


template <typename PixelT>
class CrowdedFieldMatrix {
public:
    CrowdedFieldMatrix(const afw::image::Exposure<PixelT>& exposure,
                       ndarray::Array<double const, 1>  &x,
                       ndarray::Array<double const, 1>  &y);

    CrowdedFieldMatrix(const afw::image::Exposure<PixelT> &exposure,
                       afw::table::SourceCatalog *catalog,
                       afw::table::Key<float> fluxKey);

    static void _addSource(const afw::image::Exposure<PixelT> &exposure,
                           std::vector<Eigen::Triplet<PixelT>> &matrixEntries,
                           int nStar, double  x, double y);

    static std::vector<Eigen::Triplet<PixelT>> _makeMatrixEntries(
                       const afw::image::Exposure<PixelT> &exposure,
                       ndarray::Array<double const, 1> &x,
                       ndarray::Array<double const, 1> &y);

    static std::vector<Eigen::Triplet<PixelT>> _makeMatrixEntries(
                       const afw::image::Exposure<PixelT> &exposure,
                       afw::table::SourceCatalog *catalog);

    Eigen::Matrix<PixelT, Eigen::Dynamic, 1> solve();

    std::map<int, int> renameMatrixRows();
    const std::list<std::tuple<int, int, PixelT>> getMatrixEntries();
    const Eigen::Matrix<PixelT, Eigen::Dynamic, 1> makeDataVector();
    const Eigen::Matrix<PixelT, Eigen::Dynamic, 1> getDataVector();
    const std::map<int, int> getPixelMapping();
    const std::vector<std::tuple<int, int, PixelT>> getDebug();

private:

    const std::tuple<int, int> _setNrowsNcols();
    const afw::image::Exposure<PixelT> _exposure;
    afw::table::SourceCatalog *_catalog;
    afw::table::Key<float> _fluxKey;
    std::vector<Eigen::Triplet<PixelT>> _matrixEntries;
    Eigen::Matrix<PixelT, Eigen::Dynamic, 1> _dataVector;
    std::map<int, int> _pixelMapping;
    std::vector<std::tuple<int, int, PixelT>> _debugXYValues;
    int _nStars;
    int _nRows;
    int _nColumns;


};

} // namespace crowd
} // namespace pipe
} // namespace lsst
