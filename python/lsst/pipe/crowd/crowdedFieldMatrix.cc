
#include <Eigen/Sparse>
#include <vector>

#include "pybind11/pybind11.h"
#include "ndarray/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/eigen.h"

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

    clsCrowdedFieldMatrix.def(py::init<const afw::image::Exposure<float> &,
                                       ndarray::Array<double const, 1> &,
                                        ndarray::Array<double const, 1> &>(),
                              "exposure"_a, "x"_a, "y"_a);

    clsCrowdedFieldMatrix.def(py::init<const afw::image::Exposure<float> &,
                                       afw::table::SourceCatalog *,
                                       afw::table::Key<float> >(),
                              "exposure"_a, "sourceCatalog"_a, "fluxKey"_a);
            
    clsCrowdedFieldMatrix.def("_addSource", &CrowdedFieldMatrix<float>::_addSource);

    clsCrowdedFieldMatrix.def("solve", &CrowdedFieldMatrix<float>::solve);

    clsCrowdedFieldMatrix.def("renameMatrixRows", &CrowdedFieldMatrix<float>::renameMatrixRows);
    clsCrowdedFieldMatrix.def("getMatrixEntries", &CrowdedFieldMatrix<float>::getMatrixEntries);
    clsCrowdedFieldMatrix.def("getDataVector", &CrowdedFieldMatrix<float>::getDataVector);
    clsCrowdedFieldMatrix.def("getPixelMapping", &CrowdedFieldMatrix<float>::getPixelMapping);
    clsCrowdedFieldMatrix.def("getDebug", &CrowdedFieldMatrix<float>::getDebug);


}
}
}
}
