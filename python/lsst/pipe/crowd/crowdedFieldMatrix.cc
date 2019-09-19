
#include <Eigen/Sparse>
#include <vector>

#include "pybind11/pybind11.h"
#include "ndarray/pybind11.h"
#include "pybind11/stl.h"

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

    clsCrowdedFieldMatrix.def(py::init<const afw::image::Exposure<float> &>(),
                              "exposure"_a);
            
    clsCrowdedFieldMatrix.def("addSource", &CrowdedFieldMatrix<float>::addSource);
    clsCrowdedFieldMatrix.def("addSources", &CrowdedFieldMatrix<float>::addSources);

    clsCrowdedFieldMatrix.def("solve", &CrowdedFieldMatrix<float>::solve);

    clsCrowdedFieldMatrix.def("renameMatrixRows", &CrowdedFieldMatrix<float>::renameMatrixRows);
    clsCrowdedFieldMatrix.def("getMatrixEntries", &CrowdedFieldMatrix<float>::getMatrixEntries);


}
}
}
}
