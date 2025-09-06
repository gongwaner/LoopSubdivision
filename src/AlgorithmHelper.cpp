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

    vtkSmartPointer<vtkCellArray> GetTriangleTopologyAsCellArray(const std::vector<unsigned int>& triangleVids)
    {
        auto cells = vtkSmartPointer<vtkCellArray>::New();
        auto connectivity = vtkSmartPointer<vtkIdTypeArray>::New();

        const auto numTriangles = triangleVids.size() / 3;

        //each triangle needs 3 vertices and a leading '3' for the count.
        connectivity->SetNumberOfComponents(1);
        connectivity->SetNumberOfValues(numTriangles * 4);

        vtkIdType* connData = connectivity->GetPointer(0);
        vtkIdType connIdx = 0;

        //fill the connectivity array directly
        for(int i = 0; i < numTriangles; i++)
        {
            connData[connIdx++] = 3;//points number in the cell
            connData[connIdx++] = triangleVids[i * 3];
            connData[connIdx++] = triangleVids[i * 3 + 1];
            connData[connIdx++] = triangleVids[i * 3 + 2];
        }

        cells->SetCells(numTriangles, connectivity);

        return cells;
    }

    vtkSmartPointer<vtkCellArray> GetTriangleTopologyAsCellArray(const std::vector<unsigned int>& triangleV0Ids,
                                                                 const std::vector<unsigned int>& triangleV1Ids,
                                                                 const std::vector<unsigned int>& triangleV2Ids)
    {
        auto cells = vtkSmartPointer<vtkCellArray>::New();
        auto connectivity = vtkSmartPointer<vtkIdTypeArray>::New();

        const auto numTriangles = triangleV0Ids.size();

        //each triangle needs 3 vertices and a leading '3' for the count.
        connectivity->SetNumberOfComponents(1);
        connectivity->SetNumberOfValues(numTriangles * 4);

        vtkIdType* connData = connectivity->GetPointer(0);
        vtkIdType connIdx = 0;

        //fill the connectivity array directly
        for(int i = 0; i < numTriangles; i++)
        {
            connData[connIdx++] = 3;//points number in the cell
            connData[connIdx++] = triangleV0Ids[i];
            connData[connIdx++] = triangleV1Ids[i];
            connData[connIdx++] = triangleV2Ids[i];
        }

        cells->SetCells(numTriangles, connectivity);

        return cells;
    }
}
