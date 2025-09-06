#pragma once

#include <filesystem>

#include <vtkSmartPointer.h>


class vtkPolyData;
class vtkCellArray;

namespace AlgorithmHelper
{
    std::filesystem::path GetDataDir();

    std::vector<double> GetPointsAsFlatVector(vtkPolyData* mesh);

    vtkSmartPointer<vtkCellArray> GetTriangleTopologyAsCellArray(const std::vector<unsigned int>& triangleVids);

    vtkSmartPointer<vtkCellArray> GetTriangleTopologyAsCellArray(const std::vector<unsigned int>& triangleV0Ids,
                                                                 const std::vector<unsigned int>& triangleV1Ids,
                                                                 const std::vector<unsigned int>& triangleV2Ids);
}
