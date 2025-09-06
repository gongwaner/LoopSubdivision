#include "LoopSubdivisionSoA.h"

#include <unordered_set>
#include <numbers>//for PI

#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkTriangle.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include "IOUtil.h"

#include "AlgorithmHelper.h"
#include "PairHash.h"


namespace AlgorithmSoA
{
    struct InputMeshData
    {
        //points pos
        std::vector<double> positionsX;
        std::vector<double> positionsY;
        std::vector<double> positionsZ;

        //edge end vids
        std::vector<unsigned int> edgeV0Vids;
        std::vector<unsigned int> edgeV1Vids;

        //triangle vertex ids
        std::vector<unsigned int> triangleV0Ids;
        std::vector<unsigned int> triangleV1Ids;
        std::vector<unsigned int> triangleV2Ids;

        //triangle edge ids
        std::vector<unsigned int> triangleE0Ids;
        std::vector<unsigned int> triangleE1Ids;
        std::vector<unsigned int> triangleE2Ids;

        //neighbors
        std::vector<unsigned int> vertexNeighborVidsFlatVector;
        std::vector<unsigned int> vertexNeighborOffsetVector;
        std::vector<unsigned int> edgeNeighborVidsFlatVector;
        std::vector<unsigned int> edgeNeighborOffsetVector;
    };

    void GetMeshPoints(vtkPolyData* mesh,
                       std::vector<double>& positionsX,
                       std::vector<double>& positionsY,
                       std::vector<double>& positionsZ)
    {
        const auto pointsCnt = mesh->GetNumberOfPoints();
        positionsX = std::vector<double>(pointsCnt);
        positionsY = std::vector<double>(pointsCnt);
        positionsZ = std::vector<double>(pointsCnt);

        for(auto vid = 0; vid < pointsCnt; vid++)
        {
            auto point = mesh->GetPoint(vid);
            positionsX[vid] = point[0];
            positionsY[vid] = point[1];
            positionsZ[vid] = point[2];
        }
    }

    void InitializeEdgeTable(vtkPolyData* mesh,
                             std::vector<unsigned int>& edgeV0Vids, std::vector<unsigned int>& edgeV1Vids,
                             std::vector<unsigned int>& triangleV0Ids, std::vector<unsigned int>& triangleV1Ids, std::vector<unsigned int>& triangleV2Ids,
                             std::unordered_map<std::pair<unsigned int, unsigned int>, unsigned int, PairHash>& vidsToEdgeMap,
                             std::vector<unsigned int>& boundaryEdgeVidsFlatVector,
                             std::unordered_set<unsigned int>& boundaryEdgeIdsSet)
    {
        edgeV0Vids.clear();
        edgeV1Vids.clear();
        vidsToEdgeMap.clear();

        const auto cellsCnt = mesh->GetNumberOfCells();
        triangleV0Ids = std::vector<unsigned int>(cellsCnt);
        triangleV1Ids = std::vector<unsigned int>(cellsCnt);
        triangleV2Ids = std::vector<unsigned int>(cellsCnt);

        std::unordered_map<std::pair<unsigned int, unsigned int>, unsigned int, PairHash> edgeVidsCountMap;

        //iterate through all triangles
        vtkCellArray* polys = mesh->GetPolys();
        polys->InitTraversal();

        auto pointIds = vtkSmartPointer<vtkIdList>::New();

        size_t cellId = 0;
        while(polys->GetNextCell(pointIds))
        {
            if(pointIds->GetNumberOfIds() != 3)//a triangle has 3 points
                continue;

            triangleV0Ids[cellId] = pointIds->GetId(0);
            triangleV1Ids[cellId] = pointIds->GetId(1);
            triangleV2Ids[cellId] = pointIds->GetId(2);

            //get the three edges of the current triangle
            for(int i = 0; i < 3; ++i)
            {
                const auto v1 = pointIds->GetId(i);
                const auto v2 = pointIds->GetId((i + 1) % 3);

                //create a canonical edge key
                const auto key = (v1 < v2) ? std::make_pair(v1, v2) : std::make_pair(v2, v1);
                edgeVidsCountMap[key]++;

                if(!vidsToEdgeMap.count(key))
                {
                    vidsToEdgeMap[key] = edgeV0Vids.size();

                    edgeV0Vids.push_back(key.first);
                    edgeV1Vids.push_back(key.second);
                }
            }

            ++cellId;
        }

        //construct boundary edge ids
        boundaryEdgeVidsFlatVector.clear();
        boundaryEdgeIdsSet.clear();
        for(const auto& [edgeVids, count]: edgeVidsCountMap)
        {
            if(count == 1)
            {
                const auto eid = vidsToEdgeMap.at(edgeVids);
                boundaryEdgeVidsFlatVector.push_back(edgeVids.first);
                boundaryEdgeVidsFlatVector.push_back(edgeVids.second);
                boundaryEdgeIdsSet.insert(eid);
            }
        }
    }

    void InitializeAdjacencyList(const size_t pointsCnt,
                                 const std::vector<unsigned int>& edgeStartVids, const std::vector<unsigned int>& edgeEndVids,
                                 const std::vector<unsigned int>& boundaryEdgeVidsFlatVector,
                                 std::vector<unsigned int>& vertexNeighborVidsFlatVector,
                                 std::vector<unsigned int>& vertexNeighborOffsetVector)
    {
        std::vector<std::vector<int>> adjacencyMatrix(pointsCnt);
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
        const auto edgeCnt = edgeStartVids.size();
        for(int i = 0; i < edgeCnt; ++i)
        {
            const auto startVid = edgeStartVids[i];
            const auto endVid = edgeEndVids[i];

            //don't count boundary vids' neighbors
            if(!isBoundaryVid[startVid])
                adjacencyMatrix[startVid].push_back(endVid);

            if(!isBoundaryVid[endVid])
                adjacencyMatrix[endVid].push_back(startVid);
        }

        //flatten adjacency matrix
        vertexNeighborVidsFlatVector.clear();
        vertexNeighborOffsetVector = std::vector<unsigned int>(pointsCnt + 1);
        vertexNeighborOffsetVector[0] = 0;//always starts with 0

        for(auto vid = 0; vid < pointsCnt; ++vid)
        {
            const auto& neighborVids = adjacencyMatrix[vid];
            for(const auto neighborVid: neighborVids)
            {
                vertexNeighborVidsFlatVector.push_back(neighborVid);
            }

            vertexNeighborOffsetVector[vid + 1] = vertexNeighborOffsetVector[vid] + neighborVids.size();
        }
    }

    void ProcessTriangles(const size_t edgesCnt,
                          const std::vector<unsigned int>& triangleV0Ids,
                          const std::vector<unsigned int>& triangleV1Ids,
                          const std::vector<unsigned int>& triangleV2Ids,
                          const std::unordered_set<unsigned int>& boundaryEdgeIdsSet,
                          const std::unordered_map<std::pair<unsigned int, unsigned int>, unsigned int, PairHash>& vidsToEdgeMap,
                          std::vector<unsigned int>& triangleE0Ids,
                          std::vector<unsigned int>& triangleE1Ids,
                          std::vector<unsigned int>& triangleE2Ids,
                          std::vector<unsigned int>& edgeNeighborVidsFlatVector,
                          std::vector<unsigned int>& edgeNeighborOffsetVector)
    {
        std::vector<std::vector<unsigned int>> edgeNeighborVidsVec(edgesCnt);

        const auto cellsCnt = triangleV0Ids.size();
        triangleE0Ids = std::vector<unsigned int>(cellsCnt);
        triangleE1Ids = std::vector<unsigned int>(cellsCnt);
        triangleE2Ids = std::vector<unsigned int>(cellsCnt);

        //iterate through all triangles
        for(auto cellId = 0; cellId < cellsCnt; ++cellId)
        {
            const unsigned int pointIds[3]
            {
                triangleV0Ids[cellId],
                triangleV1Ids[cellId],
                triangleV2Ids[cellId],
            };

            //get the three edges of the current triangle
            unsigned int eids[3];
            for(int i = 0; i < 3; ++i)
            {
                const auto v1 = pointIds[i];
                const auto v2 = pointIds[(i + 1) % 3];

                //create a canonical edge key
                const auto key = (v1 < v2) ? std::make_pair(v1, v2) : std::make_pair(v2, v1);
                const auto eid = vidsToEdgeMap.at(key);
                eids[i] = eid;

                if(boundaryEdgeIdsSet.contains(eid))
                    continue;

                const auto v3 = pointIds[(i + 2) % 3];
                edgeNeighborVidsVec[eid].push_back(v3);
            }

            triangleE0Ids[cellId] = eids[0];
            triangleE1Ids[cellId] = eids[1];
            triangleE2Ids[cellId] = eids[2];
        }

        //flatten edgeNeighborVidsVec
        edgeNeighborVidsFlatVector.clear();
        edgeNeighborOffsetVector = std::vector<unsigned int>(edgesCnt + 1);

        for(auto eid = 0; eid < edgesCnt; ++eid)
        {
            const auto& neighborVids = edgeNeighborVidsVec[eid];
            for(const auto neighborVid: neighborVids)
            {
                edgeNeighborVidsFlatVector.push_back(neighborVid);
            }

            edgeNeighborOffsetVector[eid + 1] = edgeNeighborOffsetVector[eid] + neighborVids.size();
        }
    }

    InputMeshData GetInputMeshData(vtkPolyData* mesh)
    {
        std::vector<double> positionsX;
        std::vector<double> positionsY;
        std::vector<double> positionsZ;
        GetMeshPoints(mesh, positionsX, positionsY, positionsZ);

        std::vector<unsigned int> edgeV0Vids;
        std::vector<unsigned int> edgeV1Vids;
        std::vector<unsigned int> triangleV0Ids;
        std::vector<unsigned int> triangleV1Ids;
        std::vector<unsigned int> triangleV2Ids;
        std::unordered_map<std::pair<unsigned int, unsigned int>, unsigned int, PairHash> vidsToEdgeMap;
        std::vector<unsigned int> boundaryEdgeVidsFlatVector;
        std::unordered_set<unsigned int> boundaryEdgeIdsSet;
        InitializeEdgeTable(mesh, edgeV0Vids, edgeV1Vids, triangleV0Ids, triangleV1Ids, triangleV2Ids,
                            vidsToEdgeMap, boundaryEdgeVidsFlatVector, boundaryEdgeIdsSet);

        std::vector<unsigned int> vertexNeighborVidsFlatVector;
        std::vector<unsigned int> vertexNeighborOffsetVector;
        InitializeAdjacencyList(mesh->GetNumberOfPoints(), edgeV0Vids, edgeV1Vids, boundaryEdgeVidsFlatVector,
                                vertexNeighborVidsFlatVector, vertexNeighborOffsetVector);

        const auto edgesCnt = edgeV0Vids.size();
        std::vector<unsigned int> triangleE0Ids;
        std::vector<unsigned int> triangleE1Ids;
        std::vector<unsigned int> triangleE2Ids;
        std::vector<unsigned int> edgeNeighborVidsFlatVector;
        std::vector<unsigned int> edgeNeighborOffsetVector;
        ProcessTriangles(edgesCnt, triangleV0Ids, triangleV1Ids, triangleV2Ids,
                         boundaryEdgeIdsSet, vidsToEdgeMap,
                         triangleE0Ids, triangleE1Ids, triangleE2Ids,
                         edgeNeighborVidsFlatVector, edgeNeighborOffsetVector);

        InputMeshData meshData;
        meshData.positionsX = positionsX;
        meshData.positionsY = positionsY;
        meshData.positionsZ = positionsZ;
        meshData.edgeV0Vids = edgeV0Vids;
        meshData.edgeV1Vids = edgeV1Vids;
        meshData.triangleV0Ids = triangleV0Ids;
        meshData.triangleV1Ids = triangleV1Ids;
        meshData.triangleV2Ids = triangleV2Ids;
        meshData.triangleE0Ids = triangleE0Ids;
        meshData.triangleE1Ids = triangleE1Ids;
        meshData.triangleE2Ids = triangleE2Ids;
        meshData.vertexNeighborVidsFlatVector = vertexNeighborVidsFlatVector;
        meshData.vertexNeighborOffsetVector = vertexNeighborOffsetVector;
        meshData.edgeNeighborVidsFlatVector = edgeNeighborVidsFlatVector;
        meshData.edgeNeighborOffsetVector = edgeNeighborOffsetVector;

        return meshData;
    }

    void GetUpdatedPointsFlatVector(const InputMeshData& meshData,
                                    std::vector<double>& updatedPositionsX,
                                    std::vector<double>& updatedPositionsY,
                                    std::vector<double>& updatedPositionsZ)
    {
        //step1: iterate through edges and add new vertices
        const auto pointsCnt = meshData.positionsX.size();
        const auto edgesCnt = meshData.edgeV0Vids.size();

        constexpr double endpointsWeight = 0.375;//3/8
        constexpr double twoPI = 2.0 * std::numbers::pi;

        const auto newPointsCnt = pointsCnt + edgesCnt;
        updatedPositionsX = std::vector<double>(newPointsCnt);
        updatedPositionsY = std::vector<double>(newPointsCnt);
        updatedPositionsZ = std::vector<double>(newPointsCnt);

        for(auto edgeId = 0; edgeId < edgesCnt; ++edgeId)
        {
            const auto vid0 = meshData.edgeV0Vids[edgeId];
            const auto vid1 = meshData.edgeV1Vids[edgeId];
            const auto index = pointsCnt + edgeId;

            updatedPositionsX[index] = meshData.positionsX[vid0] + meshData.positionsX[vid1];
            updatedPositionsY[index] = meshData.positionsY[vid0] + meshData.positionsY[vid1];
            updatedPositionsZ[index] = meshData.positionsZ[vid0] + meshData.positionsZ[vid1];

            const auto neighborVidStartIndex = meshData.edgeNeighborOffsetVector[edgeId];
            const auto neighborVidEndIndex = meshData.edgeNeighborOffsetVector[edgeId + 1];

            if(neighborVidStartIndex == neighborVidEndIndex)//boundary edge
            {
                updatedPositionsX[index] *= 0.5;
                updatedPositionsY[index] *= 0.5;
                updatedPositionsZ[index] *= 0.5;
            }
            else//interior edge
            {
                updatedPositionsX[index] *= endpointsWeight;
                updatedPositionsY[index] *= endpointsWeight;
                updatedPositionsZ[index] *= endpointsWeight;

                for(auto n = neighborVidStartIndex; n < neighborVidEndIndex; ++n)
                {
                    constexpr double neighborWeight = 0.125;
                    const auto neighborVid = meshData.edgeNeighborVidsFlatVector[n];

                    updatedPositionsX[index] += neighborWeight * meshData.positionsX[neighborVid];
                    updatedPositionsY[index] += neighborWeight * meshData.positionsY[neighborVid];
                    updatedPositionsZ[index] += neighborWeight * meshData.positionsZ[neighborVid];
                }
            }
        }

        //step2: update old vertices pos
        for(auto vid = 0; vid < pointsCnt; vid++)
        {
            const auto startIndex = meshData.vertexNeighborOffsetVector[vid];
            const auto endIndex = meshData.vertexNeighborOffsetVector[vid + 1];
            const auto n = endIndex - startIndex;

            double beta = 0.0;
            if(n < 3)
            {
                //boundary vertex
                //	p0′=3/4 p0 + 1/8 (p1 + p2)
                updatedPositionsX[vid] = 0.75 * meshData.positionsX[vid];
                updatedPositionsY[vid] = 0.75 * meshData.positionsY[vid];
                updatedPositionsZ[vid] = 0.75 * meshData.positionsZ[vid];

                beta = 0.125;//1/8
            }
            else
            {
                //interior vertex
                //β=1/n {5/8−[3/8 + 1/4 cos(2π/n) ]^2 }
                //v_old' = (1−n⋅β)⋅v_old+Σ_(j=1)^n (β⋅v_j )
                constexpr double fiveOnEight = 0.625;

                const auto inner_bracket = endpointsWeight + (1.0 / 4.0) * cos(twoPI / n);
                const auto bracket_squared = inner_bracket * inner_bracket;
                beta = (1.0 / n) * (fiveOnEight - bracket_squared);

                const auto factor = (1.0 - n * beta);
                updatedPositionsX[vid] = factor * meshData.positionsX[vid];
                updatedPositionsY[vid] = factor * meshData.positionsY[vid];
                updatedPositionsZ[vid] = factor * meshData.positionsZ[vid];
            }

            for(auto index = startIndex; index < endIndex; ++index)
            {
                const auto neighborVid = meshData.vertexNeighborVidsFlatVector[index];

                updatedPositionsX[vid] += beta * meshData.positionsX[neighborVid];
                updatedPositionsY[vid] += beta * meshData.positionsY[neighborVid];
                updatedPositionsZ[vid] += beta * meshData.positionsZ[neighborVid];
            }
        }
    }

    void GetUpdatedTriangleVids(const InputMeshData& meshData,
                                std::vector<unsigned int>& updatedTriangleV0Ids,
                                std::vector<unsigned int>& updatedTriangleV1Ids,
                                std::vector<unsigned int>& updatedTriangleV2Ids)
    {
        constexpr int subTrianglesCnt = 4;
        const auto cellsCnt = meshData.triangleV0Ids.size();
        const auto updatedCellsCnt = cellsCnt * subTrianglesCnt;
        const unsigned int pointsCnt = meshData.positionsX.size();

        updatedTriangleV0Ids = std::vector<unsigned int>(updatedCellsCnt);
        updatedTriangleV1Ids = std::vector<unsigned int>(updatedCellsCnt);
        updatedTriangleV2Ids = std::vector<unsigned int>(updatedCellsCnt);

        for(auto cellID = 0; cellID < cellsCnt; ++cellID)
        {
            //newVertices internally maps eid->points. eg, eid 0->newVertices[0]
            //here we need to know the index of point in updatedPoints
            //the index would be eid+oldpoints.count
            const unsigned int newVids[3]
            {
                meshData.triangleE0Ids[cellID] + pointsCnt,
                meshData.triangleE1Ids[cellID] + pointsCnt,
                meshData.triangleE2Ids[cellID] + pointsCnt
            };

            auto index = cellID * subTrianglesCnt;

            //interior triangle: {N0, N1, N2}
            updatedTriangleV0Ids[index] = newVids[0];
            updatedTriangleV1Ids[index] = newVids[1];
            updatedTriangleV2Ids[index] = newVids[2];
            index++;

            //corner triangles: {V0, N0, N2}, {V1, N1, N0}, {V2, N2, N1}
            updatedTriangleV0Ids[index] = meshData.triangleV0Ids[cellID];
            updatedTriangleV1Ids[index] = newVids[0];
            updatedTriangleV2Ids[index] = newVids[2];
            index++;

            updatedTriangleV0Ids[index] = meshData.triangleV1Ids[cellID];
            updatedTriangleV1Ids[index] = newVids[1];
            updatedTriangleV2Ids[index] = newVids[0];
            index++;

            updatedTriangleV0Ids[index] = meshData.triangleV2Ids[cellID];
            updatedTriangleV1Ids[index] = newVids[2];
            updatedTriangleV2Ids[index] = newVids[1];
        }
    }

    vtkSmartPointer<vtkPolyData> GetLoopSubdivisionMesh(vtkPolyData* originalMesh, const int iteration)
    {
        auto currentMesh = vtkSmartPointer<vtkPolyData>::New();
        currentMesh->DeepCopy(originalMesh);

        for(auto iter = 0; iter < iteration; ++iter)
        {
            //pre-process
            const auto meshData = GetInputMeshData(currentMesh);

            //subdivision calculation
            std::vector<double> updatedPositionsX;
            std::vector<double> updatedPositionsY;
            std::vector<double> updatedPositionsZ;
            GetUpdatedPointsFlatVector(meshData, updatedPositionsX, updatedPositionsY, updatedPositionsZ);

            std::vector<unsigned int> updatedTriangleV0Ids;
            std::vector<unsigned int> updatedTriangleV1Ids;
            std::vector<unsigned int> updatedTriangleV2Ids;
            GetUpdatedTriangleVids(meshData, updatedTriangleV0Ids, updatedTriangleV1Ids, updatedTriangleV2Ids);

            //update mesh topology
            const auto newPointsCnt = meshData.positionsX.size() + meshData.edgeV0Vids.size();
            auto vtkPointsVec = vtkSmartPointer<vtkPoints>::New();
            for(auto vid = 0; vid < newPointsCnt; ++vid)
                vtkPointsVec->InsertNextPoint(updatedPositionsX[vid], updatedPositionsY[vid], updatedPositionsZ[vid]);

            auto cells = AlgorithmHelper::GetTriangleTopologyAsCellArray(updatedTriangleV0Ids, updatedTriangleV1Ids, updatedTriangleV2Ids);

            currentMesh->SetPoints(vtkPointsVec);
            currentMesh->SetPolys(cells);
            currentMesh->Modified();
        }

        return currentMesh;
    }
}
