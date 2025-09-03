#include "AlgorithmHelper.h"

#include <vtkPolyData.h>


namespace AlgorithmHelper
{
    std::filesystem::path GetDataDir()
    {
        std::filesystem::path cmakeDir{CMAKE_SOURCE_DIR};
        return cmakeDir / "data";
    }

    std::vector<double> GetPointsAsFlatVector(vtkPolyData* mesh)
    {
        const auto pointsCnt = mesh->GetNumberOfPoints();
        std::vector<double> points(pointsCnt * 3);

        for(auto vid = 0; vid < pointsCnt; vid++)
        {
            auto point = mesh->GetPoint(vid);
            for(int i = 0; i < 3; i++)
                points[vid * 3 + i] = point[i];
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
