#include "LoopSubdivision.h"

#include <unordered_set>
#include <numbers>//for PI

#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkTriangle.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include "IOUtil.h"

#include "AlgorithmHelper.h"


namespace Algorithm
{
    struct InputMeshData
    {
        std::vector<double> originalPointsFlatVector;
        std::vector<int> edgeVidsFlatVector;
        std::vector<int> triangleVidsFlatVector;
        std::vector<std::vector<int>> adjacencyMatrix;
        std::unordered_set<int> boundaryEidsSet;
        std::vector<std::vector<int>> edgeNeighborVidsVector;
        std::vector<int> triangleEidsFlatVector;
        int pointsCnt;
        int edgesCnt;
        int cellsCnt;
    };

    void InitializeEdgeTable(vtkPolyData* mesh, std::vector<int>& edgeVidsFlatVector,
                             std::vector<int>& triangleVidsFlatVector,
                             std::unordered_map<std::pair<int, int>, int, PairHash>& vidsToEdgeMap)
    {
        edgeVidsFlatVector.clear();
        vidsToEdgeMap.clear();

        const auto cellsCnt = mesh->GetNumberOfCells();
        triangleVidsFlatVector = std::vector<int>(cellsCnt * 3);

        //iterate through all triangles
        for(auto cellId = 0; cellId < cellsCnt; ++cellId)
        {
            vtkCell* cell = mesh->GetCell(cellId);
            vtkIdList* pointIds = cell->GetPointIds();

            //get the three edges of the current triangle
            for(int i = 0; i < 3; ++i)
            {
                const auto v1 = pointIds->GetId(i);
                const auto v2 = pointIds->GetId((i + 1) % 3);

                //create a canonical edge key
                const auto key = (v1 < v2) ? std::make_pair(v1, v2) : std::make_pair(v2, v1);
                if(!vidsToEdgeMap.count(key))
                {
                    vidsToEdgeMap[key] = edgeVidsFlatVector.size() / 2;

                    edgeVidsFlatVector.push_back(key.first);
                    edgeVidsFlatVector.push_back(key.second);
                }

                triangleVidsFlatVector[cellId * 3 + i] = v1;
            }
        }
    }

    void InitializeAdjacencyMatrix(vtkPolyData* mesh, const std::vector<int>& edgeVidsFlatVector,
                                   const std::vector<int>& boundaryEdgeVidsFlatVector,
                                   std::vector<std::vector<int>>& adjacencyMatrix)
    {
        const auto pointsCnt = mesh->GetNumberOfPoints();
        adjacencyMatrix = std::vector<std::vector<int>>(pointsCnt);

        std::vector<bool> isBoundaryVid(pointsCnt, false);

        //process all boundary edges
        for(auto i = 0; i < boundaryEdgeVidsFlatVector.size() / 2; ++i)
        {
            const auto startVid = boundaryEdgeVidsFlatVector[i * 2];
            const auto endVid = boundaryEdgeVidsFlatVector[i * 2 + 1];

            isBoundaryVid[startVid] = true;
            isBoundaryVid[endVid] = true;

            adjacencyMatrix[startVid].push_back(endVid);
            adjacencyMatrix[endVid].push_back(startVid);
        }

        //traverse all edges
        const auto edgeCnt = edgeVidsFlatVector.size() / 2;
        for(int i = 0; i < edgeCnt; ++i)
        {
            const auto startVid = edgeVidsFlatVector[i * 2];
            const auto endVid = edgeVidsFlatVector[i * 2 + 1];

            //don't count boundary vids' neighbors
            if(!isBoundaryVid[startVid])
                adjacencyMatrix[startVid].push_back(endVid);

            if(!isBoundaryVid[endVid])
                adjacencyMatrix[endVid].push_back(startVid);
        }
    }

    void ProcessTriangles(const std::vector<int>& triangleVidsFlatVector,
                          const std::vector<int>& edgeVidsFlatVector,
                          const std::vector<int>& boundaryEdgeVidsFlatVector,
                          const std::unordered_map<std::pair<int, int>, int, PairHash>& vidsToEdgeMap,
                          std::unordered_set<int>& boundaryEidsSet,
                          std::vector<std::vector<int>>& edgeNeighborVidsVec, std::vector<int>& triangleEidsFlatVec)
    {
        const auto edgeCnt = edgeVidsFlatVector.size() / 2;
        edgeNeighborVidsVec = std::vector<std::vector<int>>(edgeCnt);

        const auto cellsCnt = triangleVidsFlatVector.size() / 3;
        triangleEidsFlatVec = std::vector<int>(cellsCnt * 3);

        //construct boundary eids
        if(!boundaryEdgeVidsFlatVector.empty())
        {
            boundaryEidsSet.clear();

            for(auto i = 0; i < boundaryEdgeVidsFlatVector.size() / 2; ++i)
            {
                const auto startVid = boundaryEdgeVidsFlatVector[i * 2];
                const auto endVid = boundaryEdgeVidsFlatVector[i * 2 + 1];

                const auto key = (startVid < endVid) ? std::make_pair(startVid, endVid) : std::make_pair(endVid, startVid);
                const auto eid = vidsToEdgeMap.at(key);
                boundaryEidsSet.insert(eid);
            }
        }

        //iterate through all triangles
        for(auto cellId = 0; cellId < cellsCnt; ++cellId)
        {
            const int pointIds[3]
            {
                triangleVidsFlatVector[cellId * 3],
                triangleVidsFlatVector[cellId * 3 + 1],
                triangleVidsFlatVector[cellId * 3 + 2]
            };

            //get the three edges of the current triangle
            for(int i = 0; i < 3; ++i)
            {
                const auto v1 = pointIds[i];
                const auto v2 = pointIds[(i + 1) % 3];

                //create a canonical edge key
                const auto key = (v1 < v2) ? std::make_pair(v1, v2) : std::make_pair(v2, v1);
                const auto eid = vidsToEdgeMap.at(key);
                triangleEidsFlatVec[cellId * 3 + i] = eid;

                if(boundaryEidsSet.count(eid))
                    continue;

                const auto v3 = pointIds[(i + 2) % 3];
                edgeNeighborVidsVec[eid].push_back(v3);
            }
        }
    }

    InputMeshData GetInputMeshData(vtkPolyData* mesh)
    {
        const auto originalPointsFlatVector = AlgorithmHelper::GetPointsAsFlatVector(mesh);

        std::vector<int> edgeVidsFlatVector;
        std::vector<int> triangleVidsFlatVector;
        std::unordered_map<std::pair<int, int>, int, PairHash> vidsToEdgeMap;
        InitializeEdgeTable(mesh, edgeVidsFlatVector, triangleVidsFlatVector, vidsToEdgeMap);

        const auto boundaryEdgeVidsFlatVector = AlgorithmHelper::GetBoundaryEdgeVidsFlatVector(mesh);

        std::vector<std::vector<int>> adjacencyMatrix;
        InitializeAdjacencyMatrix(mesh, edgeVidsFlatVector, boundaryEdgeVidsFlatVector, adjacencyMatrix);

        std::unordered_set<int> boundaryEidsSet;
        std::vector<std::vector<int>> edgeNeighborVidsVector;
        std::vector<int> triangleEidsFlatVec;
        ProcessTriangles(triangleVidsFlatVector, edgeVidsFlatVector, boundaryEdgeVidsFlatVector, vidsToEdgeMap,
                         boundaryEidsSet, edgeNeighborVidsVector, triangleEidsFlatVec);

        InputMeshData meshData;
        meshData.originalPointsFlatVector = originalPointsFlatVector;
        meshData.edgeVidsFlatVector = edgeVidsFlatVector;
        meshData.triangleVidsFlatVector = triangleVidsFlatVector;
        meshData.adjacencyMatrix = adjacencyMatrix;
        meshData.boundaryEidsSet = boundaryEidsSet;
        meshData.edgeNeighborVidsVector = edgeNeighborVidsVector;
        meshData.triangleEidsFlatVector = triangleEidsFlatVec;
        meshData.pointsCnt = mesh->GetNumberOfPoints();
        meshData.edgesCnt = edgeVidsFlatVector.size() / 2;
        meshData.cellsCnt = mesh->GetNumberOfCells();

        return meshData;
    }

    std::vector<double> GetUpdatedPointsFlatVector(const InputMeshData& meshData)
    {
        //step1: iterate through edges and add new vertices
        const double endpointsWeight = 0.375;//3/8
        const double neighborWeight = 0.125;//1/8
        const auto twoPI = 2.0 * std::numbers::pi;

        const auto newPointsCnt = meshData.pointsCnt + meshData.edgesCnt;
        std::vector<double> updatedPointsFlatVector(newPointsCnt * 3);

        for(auto edgeId = 0; edgeId < meshData.edgesCnt; ++edgeId)
        {
            const auto vid0 = meshData.edgeVidsFlatVector[edgeId * 2];
            const auto vid1 = meshData.edgeVidsFlatVector[edgeId * 2 + 1];
            const auto index = meshData.pointsCnt + edgeId;

            for(int i = 0; i < 3; ++i)
                updatedPointsFlatVector[index * 3 + i] = meshData.originalPointsFlatVector[vid0 * 3 + i] + meshData.originalPointsFlatVector[vid1 * 3 + i];

            if(meshData.boundaryEidsSet.count(edgeId))//boundary edge
            {
                for(int i = 0; i < 3; ++i)
                    updatedPointsFlatVector[index * 3 + i] *= 0.5;
            }
            else//interior edge
            {
                for(int i = 0; i < 3; ++i)
                    updatedPointsFlatVector[index * 3 + i] *= endpointsWeight;

                for(const auto vid: meshData.edgeNeighborVidsVector[edgeId])
                {
                    for(int i = 0; i < 3; ++i)
                        updatedPointsFlatVector[index * 3 + i] += neighborWeight * meshData.originalPointsFlatVector[vid * 3 + i];
                }
            }
        }

        //step2: update old vertices pos
        for(auto vid = 0; vid < meshData.pointsCnt; vid++)
        {
            const auto& neighborVids = meshData.adjacencyMatrix[vid];
            const auto n = neighborVids.size();

            double beta = 0.0;
            if(n < 3)
            {
                //boundary vertex
                //	p0′=3/4 p0 + 1/8 (p1 + p2)
                for(int i = 0; i < 3; ++i)
                    updatedPointsFlatVector[vid * 3 + i] = 0.75 * meshData.originalPointsFlatVector[vid * 3 + i];

                beta = 0.125;//1/8
            }
            else
            {
                //interior vertex
                //β=1/n {5/8−[3/8 + 1/4 cos(2π/n) ]^2 }
                //v_old' = (1−n⋅β)⋅v_old+Σ_(j=1)^n (β⋅v_j )
                const auto inner_bracket = endpointsWeight + (1.0 / 4.0) * cos(twoPI / n);
                const auto bracket_squared = inner_bracket * inner_bracket;
                beta = (1.0 / n) * ((5.0 / 8.0) - bracket_squared);

                for(int i = 0; i < 3; ++i)
                    updatedPointsFlatVector[vid * 3 + i] = (1.0 - n * beta) * meshData.originalPointsFlatVector[vid * 3 + i];
            }

            for(const auto neighborVid: neighborVids)
            {
                for(int i = 0; i < 3; ++i)
                    updatedPointsFlatVector[vid * 3 + i] += beta * meshData.originalPointsFlatVector[neighborVid * 3 + i];
            }
        }

        return updatedPointsFlatVector;
    }

    std::vector<int> GetUpdatedTriangleVidsFlatVector(const InputMeshData& meshData)
    {
        const int subTrianglesCnt = 4;
        std::vector<int> newTrianglesVidsFlatVector(meshData.cellsCnt * subTrianglesCnt * 3);

        for(auto cellID = 0; cellID < meshData.cellsCnt; ++cellID)
        {
            //newVertices internally maps eid->points. eg, eid 0->newVertices[0]
            //here we need to know the index of point in updatedPoints
            //the index would be eid+oldpoints.count
            int newVids[3];
            for(int i = 0; i < 3; ++i)
            {
                newVids[i] = meshData.triangleEidsFlatVector[cellID * 3 + i] + meshData.pointsCnt;
            }

            auto index = cellID * subTrianglesCnt;

            //interior triangle: {N0, N1, N2}
            for(int i = 0; i < 3; ++i)
                newTrianglesVidsFlatVector[index * 3 + i] = newVids[i];
            index++;

            //corner triangles: {p0, N0, N2}, {p1, N1, N0}, {p2, N2, N1}
            newTrianglesVidsFlatVector[index * 3] = meshData.triangleVidsFlatVector[cellID * 3];
            newTrianglesVidsFlatVector[index * 3 + 1] = newVids[0];
            newTrianglesVidsFlatVector[index * 3 + 2] = newVids[2];
            index++;

            newTrianglesVidsFlatVector[index * 3] = meshData.triangleVidsFlatVector[cellID * 3 + 1];
            newTrianglesVidsFlatVector[index * 3 + 1] = newVids[1];
            newTrianglesVidsFlatVector[index * 3 + 2] = newVids[0];
            index++;

            newTrianglesVidsFlatVector[index * 3] = meshData.triangleVidsFlatVector[cellID * 3 + 2];
            newTrianglesVidsFlatVector[index * 3 + 1] = newVids[2];
            newTrianglesVidsFlatVector[index * 3 + 2] = newVids[1];
        }

        return newTrianglesVidsFlatVector;
    }

    vtkSmartPointer<vtkPolyData> GetLoopSubdivisionMesh(vtkPolyData* originalMesh, const int iteration)
    {
        auto currentMesh = vtkSmartPointer<vtkPolyData>::New();
        currentMesh->DeepCopy(originalMesh);

        for(auto iter = 0; iter < iteration; ++iter)
        {
            //pre-process
            const auto meshData = GetInputMeshData(currentMesh);

            const auto updatedPointsFlatVector = GetUpdatedPointsFlatVector(meshData);
            const auto updatedTriangleVidsFlatVector = GetUpdatedTriangleVidsFlatVector(meshData);

            //set mesh topology
            const auto newPointsCnt = meshData.pointsCnt + meshData.edgesCnt;
            auto vtkPointsVec = vtkSmartPointer<vtkPoints>::New();
            for(auto vid = 0; vid < newPointsCnt; ++vid)
                vtkPointsVec->InsertNextPoint(updatedPointsFlatVector[vid * 3 + 0], updatedPointsFlatVector[vid * 3 + 1], updatedPointsFlatVector[vid * 3 + 2]);

            auto cells = AlgorithmHelper::GetTriangleTopologyAsCellArray(updatedTriangleVidsFlatVector);

            currentMesh->SetPoints(vtkPointsVec);
            currentMesh->SetPolys(cells);
            currentMesh->Modified();
        }

        return currentMesh;
    }
}
