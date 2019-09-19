
#include "lsst/base.h"
#include "lsst/afw/image/Exposure.h"
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

    static void _addSource(const afw::image::Exposure<PixelT> &exposure, 
                           std::vector<Eigen::Triplet<PixelT>> &matrixEntries,
                           int nStar, double  x, double y);

    static std::vector<Eigen::Triplet<PixelT>> _makeMatrixEntries(
                       const afw::image::Exposure<PixelT> &exposure,
                       ndarray::Array<double const, 1> &x,
                       ndarray::Array<double const, 1> &y);

    void solve();

    std::map<int, int> renameMatrixRows();
    std::list<std::tuple<int, int, PixelT>> getMatrixEntries();
    const Eigen::Matrix<PixelT, Eigen::Dynamic, 1> makeDataMatrix();

private:

    const afw::image::Exposure<PixelT> _exposure;
    std::vector<Eigen::Triplet<PixelT>> _matrixEntries;
    std::map<int, int> _pixelMapping;
    int _nStars;
    int _nRows;
    int _nColumns;

};

} // namespace crowd
} // namespace pipe
} // namespace lsst
