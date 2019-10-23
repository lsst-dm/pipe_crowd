
#include "lsst/pipe/crowd/ParameterTracker.h"
#include "lsst/pex/exceptions/Runtime.h"

#include <vector>
#include <map>
#include <tuple>

namespace lsst {
namespace pipe {
namespace crowd {

ParameterTracker::ParameterTracker(int nParameters) :
    _nParameters(nParameters),
    _pixelMapping(std::map<std::tuple<int, int>, int>()),
    _sourceParameterMapping(std::map<std::tuple<int, int>, int>())
{
}

void ParameterTracker::addSource(int sourceId) {

    int nextParameterId;

    for(int i = 0; i < _nParameters; i++) {
        nextParameterId = _sourceParameterMapping.size();
        _sourceParameterMapping.insert({std::make_tuple(sourceId, i),
                                       nextParameterId});
    }
}

/*
 * getPixelId can mutate the mapping if a pixel isn't found,
 * because we expect this to get built up before solving and
 * then not need to be accessed again.
 */
int ParameterTracker::makePixelId(int pixelX, int pixelY) {

    int newPixelId;
    auto pixelMapEntry = _pixelMapping.find(std::make_tuple(pixelX, pixelY));

    if(pixelMapEntry != _pixelMapping.end()) {
        return pixelMapEntry->second;
    } else {
        newPixelId = _pixelMapping.size();
        _pixelMapping.insert({std::make_tuple(pixelX, pixelY), newPixelId});
        return newPixelId;
    }

}

int* ParameterTracker::getPixelId(int pixelX, int pixelY) {
    auto pixelMapEntry = _pixelMapping.find(std::make_tuple(pixelX, pixelY));
    if(pixelMapEntry != _pixelMapping.end()) {
        return &pixelMapEntry->second;
    } else {
        return NULL;
    }
}

int ParameterTracker::getSourceParameterId(int sourceId, int param) {
    auto result = _sourceParameterMapping.find(std::make_tuple(sourceId, param));
    if(result == _sourceParameterMapping.end()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, "Request for SourceParameter that does not exist.");
    }
    return result->second;
}


std::map<std::tuple<int, int>, int> ParameterTracker::getParameterMapping() {
    return _sourceParameterMapping;
}

std::map<std::tuple<int, int>, int> ParameterTracker::getPixelMapping() {
    return _pixelMapping;
}

int ParameterTracker::nRows() {
    return _pixelMapping.size();
}

int ParameterTracker::nColumns() {
    return _sourceParameterMapping.size();
}

}
}
}

