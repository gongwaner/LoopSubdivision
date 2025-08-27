#include "AlgorithmHelper.h"

#include <vtkPolyData.h>
#include <vtkPointData.h>


namespace AlgorithmHelper
{
    std::filesystem::path GetDataDir()
    {
        std::filesystem::path cmakeDir{CMAKE_SOURCE_DIR};
        return cmakeDir / "data";
    }

    std::vector<vtkVector3d> GetPoints(vtkPolyData* mesh)
    {
        const auto pointsCnt = mesh->GetNumberOfPoints();
        std::vector<vtkVector3d> points(pointsCnt);

        for(auto i = 0; i < pointsCnt; i++)
        {
            auto point = mesh->GetPoint(i);
            points[i] = vtkVector3d(point[0], point[1], point[2]);
        }

        return points;
    }

    std::unordered_map<std::pair<int, int>, int, PairHash> GetVidsToEdgeMap(const std::vector<std::pair<int, int>>& edgeVidsVector)
    {
        const auto count = edgeVidsVector.size();
        std::unordered_map<std::pair<int, int>, int, PairHash> vidsToEdgeMap;

        //traverse all edges
        for(auto i = 0; i < count; i++)
        {
            const auto startVid = edgeVidsVector[i].first;
            const auto endVid = edgeVidsVector[i].second;
            const auto key = (startVid < endVid) ? std::make_pair(startVid, endVid) : std::make_pair(endVid, startVid);

            vidsToEdgeMap[key] = i;
        }

        return vidsToEdgeMap;
    }
}
