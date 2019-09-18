
#include "pybind11/pybind11.h"
#include "ndarray/pybind11.h"

#include "lsst/afw/table/io/python.h"
#include "lsst/pipe/crowd/CrowdedFieldMatrix.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace pipe {
namespace crowd {

PYBIND11_MODULE(crowdedFieldMatrix, mod) {

    py::class_<CrowdedFieldMatrix<float>, std::shared_ptr<CrowdedFieldMatrix<float>>>
            clsCrowdedFieldMatrix(mod, "CrowdedFieldMatrix");

    clsCrowdedFieldMatrix.def(py::init<std::shared_ptr<const afw::image::Exposure<float>> &>(),
                              "exposure"_a);
            
    clsCrowdedFieldMatrix.def("addSources", &CrowdedFieldMatrix<float>::addSources);


}
}
}
}
