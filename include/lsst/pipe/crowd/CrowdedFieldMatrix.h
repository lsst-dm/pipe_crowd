
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
    CrowdedFieldMatrix(CONST_PTR(afw::image::Exposure<PixelT>) exposure);

    void addSource(double  x, double y);

    void addSources(ndarray::Array<double const, 1>  x,
               ndarray::Array<double const, 1>  y);

    std::list<std::tuple<int, int, PixelT>> getMatrixEntries();

private:

    CONST_PTR(afw::image::Exposure<PixelT>) _exposure;
    std::vector<Eigen::Triplet<PixelT>> _matrixEntries;
    int _nStars;

};

} // namespace crowd
} // namespace pipe
} // namespace lsst
