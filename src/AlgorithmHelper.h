#pragma once

#include <filesystem>
#include <unordered_map>

#include "PairHash.h"


class vtkPolyData;

namespace AlgorithmHelper
{
    std::filesystem::path GetDataDir();

    std::vector<double> GetPointsAsFlatVector(vtkPolyData* mesh);

    std::unordered_map<std::pair<int, int>, int, PairHash> GetVidsToEdgeMap(const std::vector<std::pair<int, int>>& edgeVidsVector);
}
