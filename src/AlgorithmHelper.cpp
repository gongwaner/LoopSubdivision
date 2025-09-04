#include "AlgorithmHelper.h"

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkFeatureEdges.h>

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

    std::vector<int> GetBoundaryEdgeVidsFlatVector(vtkPolyData* mesh)
    {
        const auto arrayName = "OriginalIDs";

        auto originalIDs = vtkSmartPointer<vtkIntArray>::New();
        originalIDs->SetName(arrayName);
        originalIDs->SetNumberOfComponents(1);
        originalIDs->SetNumberOfTuples(mesh->GetNumberOfPoints());

        for(vtkIdType i = 0; i < mesh->GetNumberOfPoints(); ++i)
        {
            originalIDs->SetValue(i, i);
        }

        mesh->GetPointData()->AddArray(originalIDs);

        auto featureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
        featureEdges->SetInputData(mesh);
        featureEdges->BoundaryEdgesOn();
        featureEdges->ManifoldEdgesOff();
        featureEdges->NonManifoldEdgesOff();
        featureEdges->FeatureEdgesOff();
        featureEdges->Update();

        auto boundaryEdges = featureEdges->GetOutput();
        vtkCellArray* lines = boundaryEdges->GetLines();

        const auto edgesCnt = lines->GetNumberOfCells();
        vtkDataArray* originalIDsArray = boundaryEdges->GetPointData()->GetArray(arrayName);
        if(!originalIDsArray)
        {
            std::cerr << "Error: original ids array not found in output." << std::endl;
            return {};
        }

        std::vector<int> boundaryEdgeVids(edgesCnt * 2);
        for(auto i = 0; i < edgesCnt; i++)
        {
            auto cell = vtkSmartPointer<vtkGenericCell>::New();
            boundaryEdges->GetCell(i, cell);

            if(cell->GetCellType() == VTK_LINE)
            {
                const auto startVid = originalIDsArray->GetTuple1(cell->GetPointId(0));
                const auto endVid = originalIDsArray->GetTuple1(cell->GetPointId(1));
                boundaryEdgeVids[i * 2] = startVid;
                boundaryEdgeVids[i * 2 + 1] = endVid;
            }
        }

        return boundaryEdgeVids;
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

    vtkSmartPointer<vtkCellArray> GetTriangleTopologyAsCellArray(const std::vector<int>& triangleVids)
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
}
