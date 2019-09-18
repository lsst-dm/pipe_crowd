
#include "lsst/base.h"
#include "lsst/afw/table/Exposure.h"
#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/pipe/crowd/CrowdedFieldMatrix.h"
#include "lsst/log/Log.h"


namespace lsst {
namespace pipe {
namespace crowd {


LOG_LOGGER _log = LOG_GET("pipe.crowd");

template <typename PixelT>
CrowdedFieldMatrix<PixelT>::CrowdedFieldMatrix(CONST_PTR(afw::image::Exposure<PixelT>) exposure) :
            _exposure(exposure) { };
    

template <typename PixelT>
void CrowdedFieldMatrix<PixelT>::addSources(ndarray::Array<double const, 1>  x,
                               ndarray::Array<double const, 1>  y) {
    LOGL_WARN(_log, "Got x,y values %g, %g", x[0], y[0]);
}

template class CrowdedFieldMatrix<float>;


} // namespace crowd
} // namespace pipe
} // namespace lsst
