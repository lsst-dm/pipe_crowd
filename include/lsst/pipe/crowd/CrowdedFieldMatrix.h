
#include "lsst/base.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/table/Key.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/pipe/crowd/ParameterTracker.h"
#include <Eigen/Sparse>
#include <vector>

namespace lsst {
namespace pipe {
namespace crowd {

enum SolverStatus {
    SUCCESS = 0,
    FAILURE
};

template <typename PixelT>
class CrowdedFieldMatrix {
public:
    CrowdedFieldMatrix(const afw::image::Exposure<PixelT>& exposure,
                       ndarray::Array<double const, 1>  &x,
                       ndarray::Array<double const, 1>  &y);

    CrowdedFieldMatrix(const afw::image::Exposure<PixelT> &exposure,
                       afw::table::SourceCatalog *catalog,
                       afw::table::Key<double> fluxKey,
                       bool fitCentroids = true,
                       afw::table::PointKey<double> centroidKey = afw::table::PointKey<double>());

    void _addSource(const afw::image::Exposure<PixelT> &exposure,
                           std::vector<Eigen::Triplet<PixelT>> &matrixEntries,
                           int nStar, double  x, double y,
                           PixelT estFlux=PixelT());

    std::vector<Eigen::Triplet<PixelT>> _makeMatrixEntries(
                       const afw::image::Exposure<PixelT> &exposure,
                       ndarray::Array<double const, 1> &x,
                       ndarray::Array<double const, 1> &y);

    std::vector<Eigen::Triplet<PixelT>> _makeMatrixEntries(
                       const afw::image::Exposure<PixelT> &exposure,
                       afw::table::SourceCatalog *catalog);

    SolverStatus solve();

    const std::list<std::tuple<int, int, PixelT>> getMatrixEntries();
    const Eigen::Matrix<PixelT, Eigen::Dynamic, 1> makeDataVector();
    const Eigen::Matrix<PixelT, Eigen::Dynamic, 1> getDataVector();

    const std::map<std::tuple<int, int>, int> getParameterMapping();
    const std::map<std::tuple<int, int>, int> getPixelMapping();

    int iterations();

    Eigen::Matrix<PixelT, Eigen::Dynamic, 1> result();

private:

    const afw::image::Exposure<PixelT> _exposure;
    afw::table::SourceCatalog *_catalog;
    afw::table::Key<double> _fluxKey;
    const bool _fitCentroids;
    afw::table::PointKey<double> _centroidKey;
    ParameterTracker _paramTracker;
    int _iterations;
    int _maxIterations;
    Eigen::Matrix<PixelT, Eigen::Dynamic, 1> _result;

    std::vector<Eigen::Triplet<PixelT>> _matrixEntries;
    Eigen::Matrix<PixelT, Eigen::Dynamic, 1> _dataVector;

};

} // namespace crowd
} // namespace pipe
} // namespace lsst
