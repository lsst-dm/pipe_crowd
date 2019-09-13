
#include "lsst/base.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/algorithms/ImagePsf.h"

namespace lsst {
namespace pipe {
namespace crowd {

template <typename PixelT>
class CrowdedFieldMatrix {
public:
    CrowdedFieldMatrix(CONST_PTR(afw::image::Exposure<PixelT>) exposure,
                       CONST_PTR(lsst::meas::algorithms::ImagePsf) psf);

    void addSources(ndarray::Array<double const, 1, 1>  x,
               ndarray::Array<double const, 1, 1>  y);

private:

    CONST_PTR(afw::image::Exposure<PixelT>) _exposure;

    CONST_PTR(lsst::meas::algorithms::ImagePsf) _psf;
};

} // namespace crowd
} // namespace pipe
} // namespace lsst
