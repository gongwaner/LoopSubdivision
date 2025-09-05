#pragma once

#include <filesystem>
#include <unordered_map>

#include <vtkSmartPointer.h>

#include "PairHash.h"


class vtkPolyData;
class vtkCellArray;

namespace AlgorithmHelper
{
    std::filesystem::path GetDataDir();

    std::vector<double> GetPointsAsFlatVector(vtkPolyData* mesh);

    std::vector<int> GetBoundaryEdgeVidsFlatVector(vtkPolyData* mesh);

    std::unordered_map<std::pair<int, int>, int, PairHash> GetVidsToEdgeMap(const std::vector<std::pair<int, int>>& edgeVidsVector);

    vtkSmartPointer<vtkCellArray> GetTriangleTopologyAsCellArray(const std::vector<int>& triangleVids);

    vtkSmartPointer<vtkCellArray> GetTriangleTopologyAsCellArray(const std::vector<int>& triangleV0Ids,
                                                                 const std::vector<int>& triangleV1Ids,
                                                                 const std::vector<int>& triangleV2Ids);
}
