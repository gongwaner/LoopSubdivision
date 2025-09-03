#include "LoopSubdivision.h"

#include <unordered_set>
#include <numbers>//for PI

#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkTriangle.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include "IOUtil.h"
#include "TopologyUtil.h"

#include "AlgorithmHelper.h"


namespace Algorithm
{
    void InitializeEdgeTable(vtkPolyData* mesh, std::vector<std::pair<int, int>>& edgeVidsVector,
                             std::vector<int>& triangleVidsFlatVector)
    {
        edgeVidsVector.clear();

        const auto cellsCnt = mesh->GetNumberOfCells();
        triangleVidsFlatVector = std::vector<int>(cellsCnt * 3);

        std::unordered_set<std::pair<int, int>, PairHash> visited;

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
                if(!visited.count(key))
                {
                    visited.insert(key);
                    edgeVidsVector.push_back(key);
                }

                triangleVidsFlatVector[cellId * 3 + i] = v1;
            }
        }
    }

    void InitializeAdjacencyMatrix(vtkPolyData* mesh, const std::vector<std::pair<int, int>>& edgeVidsVector,
                                   const std::vector<std::pair<int, int>>& boundaryEdgeVids, std::vector<std::vector<int>>& adjacencyMatrix)
    {
        const auto pointsCnt = mesh->GetNumberOfPoints();
        adjacencyMatrix = std::vector<std::vector<int>>(pointsCnt);

        std::vector<bool> isBoundaryVid(pointsCnt, false);

        //process all boundary edges
        for(const auto& [startVid, endVid]: boundaryEdgeVids)
        {
            isBoundaryVid[startVid] = true;
            isBoundaryVid[endVid] = true;

            adjacencyMatrix[startVid].push_back(endVid);
            adjacencyMatrix[endVid].push_back(startVid);
        }

        //traverse all edges
        for(const auto& [startVid, endVid]: edgeVidsVector)
        {
            //don't count boundary vids' neighbors
            if(!isBoundaryVid[startVid])
                adjacencyMatrix[startVid].push_back(endVid);

            if(!isBoundaryVid[endVid])
                adjacencyMatrix[endVid].push_back(startVid);
        }
    }

    void ProcessTriangles(vtkPolyData* mesh, const std::vector<std::pair<int, int>>& edgeVidsVector,
                          const std::vector<std::pair<int, int>>& boundaryEdgeVids,
                          std::unordered_set<int>& boundaryEidsSet,
                          std::vector<std::vector<int>>& edgeNeighborVidsVec, std::vector<int>& triangleEidsFlatVec)
    {
        const auto edgeCnt = edgeVidsVector.size();
        edgeNeighborVidsVec = std::vector<std::vector<int>>(edgeCnt);

        const auto cellsCnt = mesh->GetNumberOfCells();
        triangleEidsFlatVec = std::vector<int>(cellsCnt * 3);

        //construct boundary eids
        const auto vidsToEdgeMap = AlgorithmHelper::GetVidsToEdgeMap(edgeVidsVector);
        boundaryEidsSet.clear();
        for(const auto& [startVid, endVid]: boundaryEdgeVids)
        {
            const auto key = (startVid < endVid) ? std::make_pair(startVid, endVid) : std::make_pair(endVid, startVid);
            const auto eid = vidsToEdgeMap.at(key);
            boundaryEidsSet.insert(eid);
        }

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
                const auto eid = vidsToEdgeMap.at(key);
                triangleEidsFlatVec[cellId * 3 + i] = eid;

                if(boundaryEidsSet.count(eid))
                    continue;

                const auto v3 = pointIds->GetId((i + 2) % 3);
                edgeNeighborVidsVec[eid].push_back(v3);
            }
        }
    }

    vtkSmartPointer<vtkPolyData> GetLoopSubdivisionMesh(vtkPolyData* originalMesh, const int iteration)
    {
        auto currentMesh = vtkSmartPointer<vtkPolyData>::New();
        currentMesh->DeepCopy(originalMesh);

        for(auto iter = 0; iter < iteration; ++iter)
        {
            //pre-process
            const auto originalPointsFlat = AlgorithmHelper::GetPointsAsFlatVector(currentMesh);

            std::vector<std::pair<int, int>> edgeVidsVector;
            std::vector<int> triangleVidsFlatVector;
            InitializeEdgeTable(currentMesh, edgeVidsVector, triangleVidsFlatVector);

            const auto boundaryEdgeVids = TopologyUtil::GetBoundaryEdgeVids(currentMesh);

            std::vector<std::vector<int>> adjacencyMatrix;
            InitializeAdjacencyMatrix(currentMesh, edgeVidsVector, boundaryEdgeVids, adjacencyMatrix);

            std::unordered_set<int> boundaryEidsSet;
            std::vector<std::vector<int>> edgeNeighborVidsVec;
            std::vector<int> triangleEidsFlatVec;
            ProcessTriangles(currentMesh, edgeVidsVector, boundaryEdgeVids, boundaryEidsSet, edgeNeighborVidsVec, triangleEidsFlatVec);

            const auto edgesCnt = edgeVidsVector.size();

            //step1: iterate through edges and add new vertices
            const double endpointsWeight = 0.375;//3/8
            const double neighborWeight = 0.125;//1/8
            const auto twoPI = 2.0 * std::numbers::pi;

            const auto originalPointsCnt = currentMesh->GetNumberOfPoints();
            const auto newPointsCnt = originalPointsCnt + edgesCnt;
            std::vector<double> updatedPointsFlat(newPointsCnt * 3);

            for(auto edgeId = 0; edgeId < edgesCnt; ++edgeId)
            {
                const auto vid0 = edgeVidsVector[edgeId].first;
                const auto vid1 = edgeVidsVector[edgeId].second;
                const auto index = originalPointsCnt + edgeId;

                for(int i = 0; i < 3; ++i)
                    updatedPointsFlat[index * 3 + i] = originalPointsFlat[vid0 * 3 + i] + originalPointsFlat[vid1 * 3 + i];

                if(boundaryEidsSet.count(edgeId))//boundary edge
                {
                    for(int i = 0; i < 3; ++i)
                        updatedPointsFlat[index * 3 + i] *= 0.5;
                }
                else//interior edge
                {
                    for(int i = 0; i < 3; ++i)
                        updatedPointsFlat[index * 3 + i] *= endpointsWeight;

                    for(const auto vid: edgeNeighborVidsVec[edgeId])
                    {
                        for(int i = 0; i < 3; ++i)
                            updatedPointsFlat[index * 3 + i] += neighborWeight * originalPointsFlat[vid * 3 + i];
                    }
                }
            }

            //step2: update old vertices pos
            for(auto vid = 0; vid < originalPointsCnt; vid++)
            {
                const auto& neighborVids = adjacencyMatrix[vid];
                const auto n = neighborVids.size();

                double beta = 0.0;
                if(n < 3)
                {
                    //boundary vertex
                    //	p0′=3/4 p0 + 1/8 (p1 + p2)
                    for(int i = 0; i < 3; ++i)
                        updatedPointsFlat[vid * 3 + i] = 0.75 * originalPointsFlat[vid * 3 + i];

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
                        updatedPointsFlat[vid * 3 + i] = (1.0 - n * beta) * originalPointsFlat[vid * 3 + i];
                }

                for(const auto neighborVid: neighborVids)
                {
                    for(int i = 0; i < 3; ++i)
                        updatedPointsFlat[vid * 3 + i] += beta * originalPointsFlat[neighborVid * 3 + i];
                }
            }

            //step 3: update topology and create new mesh
            //points are already updated

            //triangles
            const auto cellsCnt = currentMesh->GetNumberOfCells();
            const int subTrianglesCnt = 4;
            std::vector<int> newTrianglesVidsFlat(cellsCnt * subTrianglesCnt * 3);

            for(int cellID = 0; cellID < cellsCnt; ++cellID)
            {
                //newVertices internally maps eid->points. eg, eid 0->newVertices[0]
                //here we need to know the index of point in updatedPoints
                //the index would be eid+oldpoints.count
                double newVids[3];
                for(int i = 0; i < 3; ++i)
                {
                    newVids[i] = triangleEidsFlatVec[cellID * 3 + i] + originalPointsCnt;
                }

                auto index = cellID * subTrianglesCnt;

                //interior triangle: {N0, N1, N2}
                for(int i = 0; i < 3; ++i)
                    newTrianglesVidsFlat[index * 3 + i] = newVids[i];
                index++;

                //corner triangles: {p0, N0, N2}, {p1, N1, N0}, {p2, N2, N1}
                newTrianglesVidsFlat[index * 3] = triangleVidsFlatVector[cellID * 3];
                newTrianglesVidsFlat[index * 3 + 1] = newVids[0];
                newTrianglesVidsFlat[index * 3 + 2] = newVids[2];
                index++;

                newTrianglesVidsFlat[index * 3] = triangleVidsFlatVector[cellID * 3 + 1];
                newTrianglesVidsFlat[index * 3 + 1] = newVids[1];
                newTrianglesVidsFlat[index * 3 + 2] = newVids[0];
                index++;

                newTrianglesVidsFlat[index * 3] = triangleVidsFlatVector[cellID * 3 + 2];
                newTrianglesVidsFlat[index * 3 + 1] = newVids[2];
                newTrianglesVidsFlat[index * 3 + 2] = newVids[1];
            }

            //set mesh topology
            auto vtkPointsVec = vtkSmartPointer<vtkPoints>::New();
            for(auto vid = 0; vid < newPointsCnt; ++vid)
                vtkPointsVec->InsertNextPoint(updatedPointsFlat[vid * 3 + 0], updatedPointsFlat[vid * 3 + 1], updatedPointsFlat[vid * 3 + 2]);

            auto cells = AlgorithmHelper::GetTriangleTopologyAsCellArray(newTrianglesVidsFlat);

            auto newMesh = vtkSmartPointer<vtkPolyData>::New();
            newMesh->SetPoints(vtkPointsVec);
            newMesh->SetPolys(cells);
            newMesh->Modified();

            currentMesh->DeepCopy(newMesh);
        }

        return currentMesh;
    }
}
