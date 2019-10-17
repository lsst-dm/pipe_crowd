
#include <vector>
#include <map>
#include <tuple>

namespace lsst {
namespace pipe {
namespace crowd {


class ParameterTracker {
public:
    ParameterTracker(int nParameters);

    void addSource(int sourceId);

    int makePixelId(int pixelX, int pixelY);
    int* getPixelId(int pixelX, int pixelY);

    int getSourceParameterId(int sourceId, int param);

    int nRows();
    int nColumns();

private:

    int _nParameters;
    int _nSources;
    int _nPixels;

    // image X,Y pixel -> matrix row
    std::map<std::tuple<int, int>, int> _pixelMapping;

    // sourceId, ParamId -> matrix column
    std::map<std::tuple<int, int>, int> _sourceParameterMapping;

};

} // namespace crowd
} // namespace pipe
} // namespace lsst
