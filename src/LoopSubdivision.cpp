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
    void ProcessEdges(vtkPolyData* mesh, std::vector<std::pair<int, int>>& edgeVidsVector,
                      std::vector<std::unordered_set<int>>& adjacencyMatrix)
    {
        edgeVidsVector.clear();
        adjacencyMatrix.clear();

        const auto pointsCnt = mesh->GetNumberOfPoints();
        adjacencyMatrix.resize(pointsCnt);

        //process all boundary edges
        std::vector<bool> isBoundaryVid(pointsCnt, false);
        const auto boundaryEdges = TopologyUtil::GetBoundaryEdgeVids(mesh);
        for(const auto& [startVid, endVid]: boundaryEdges)
        {
            isBoundaryVid[startVid] = true;
            isBoundaryVid[endVid] = true;

            adjacencyMatrix[startVid].insert(endVid);
            adjacencyMatrix[endVid].insert(startVid);
        }

        const auto allEdges = TopologyUtil::GetEdgeVids(mesh);
        const auto linesCnt = allEdges.size();

        //traverse all edges
        for(auto i = 0; i < linesCnt; i++)
        {
            const auto& [startVid, endVid] = allEdges[i];
            edgeVidsVector.push_back({startVid, endVid});

            //don't count boundary vids' neighbors
            if(!isBoundaryVid[startVid])
                adjacencyMatrix[startVid].insert(endVid);

            if(!isBoundaryVid[endVid])
                adjacencyMatrix[endVid].insert(startVid);
        }
    }

    void ProcessTriangles(vtkPolyData* mesh, const std::unordered_map<std::pair<int, int>, int, PairHash>& vidsToEdgeMap,
                          std::unordered_map<int, std::vector<int>>& eidToTidsMap, std::vector<std::vector<int>>& triangleEidsVec)
    {
        eidToTidsMap.clear();
        triangleEidsVec.clear();

        const auto cellsCnt = mesh->GetNumberOfCells();
        triangleEidsVec = std::vector<std::vector<int>>(cellsCnt, std::vector<int>(3, 0));

        //iterate through all triangles
        for(vtkIdType cellId = 0; cellId < cellsCnt; ++cellId)
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

                eidToTidsMap[eid].push_back(cellId);
                triangleEidsVec[cellId][i] = eid;
            }
        }
    }

    vtkSmartPointer<vtkPolyData> GetLoopSubdivisionMesh(vtkPolyData* originalMesh, const int iteration)
    {
        auto currentMesh = vtkSmartPointer<vtkPolyData>::New();
        currentMesh->DeepCopy(originalMesh);

        for(int i = 0; i < iteration; ++i)
        {
            //pre-process
            const auto originalPoints = AlgorithmHelper::GetPoints(currentMesh);

            std::vector<std::pair<int, int>> edgeVidsVector;
            std::vector<std::unordered_set<int>> adjacencyMatrix;
            ProcessEdges(currentMesh, edgeVidsVector, adjacencyMatrix);

            const auto vidsToEdgeMap = AlgorithmHelper::GetVidsToEdgeMap(edgeVidsVector);
            const auto triangleVidsVector = TopologyUtil::GetTriangleVids(currentMesh);

            std::unordered_map<int, std::vector<int>> eidToTidsMap;
            std::vector<std::vector<int>> triangleEidsVec;
            ProcessTriangles(currentMesh, vidsToEdgeMap, eidToTidsMap, triangleEidsVec);

            const auto edgesCnt = edgeVidsVector.size();

            //step1: iterate through edges and add new vertices
            const double endpointsWeight = 0.375;//3/8
            const double neighborWeight = 0.125;//1/8
            std::vector<vtkVector3d> newVertices(edgesCnt);
            for(int i = 0; i < edgesCnt; ++i)
            {
                const auto triangles = eidToTidsMap.at(i);

                const auto vid0 = edgeVidsVector[i].first;
                const auto vid1 = edgeVidsVector[i].second;

                const auto v0 = originalPoints[vid0];
                const auto v1 = originalPoints[vid1];

                vtkVector3d point;

                if(triangles.size() < 2)
                {
                    //boundary edge
                    point = 0.5 * (v0 + v1);
                }
                else
                {
                    //interior edge
                    //find the neighbors
                    std::unordered_set<int> neighborVidsSet;
                    for(const auto triangle: triangles)
                    {
                        const auto triVids = triangleVidsVector[triangle];
                        for(const auto vid: triVids)
                        {
                            if(vid == vid0 || vid == vid1)
                                continue;

                            neighborVidsSet.insert(vid);
                        }
                    }

                    point = endpointsWeight * (v0 + v1);
                    for(const auto& vid: neighborVidsSet)
                    {
                        point += neighborWeight * originalPoints[vid];
                    }
                }

                newVertices[i] = point;
            }

            //step2: update old vertices pos
            const auto pointsCnt = originalPoints.size();
            std::vector<vtkVector3d> updatedPoints(pointsCnt);

            for(auto i = 0; i < pointsCnt; i++)
            {
                const auto neighborVids = adjacencyMatrix[i];
                const auto n = neighborVids.size();

                double beta = 0.0;
                if(n < 3)
                {
                    //boundary vertex
                    //	p0′=3/4 p0 + 1/8 (p1 + p2)
                    updatedPoints[i] = 0.75 * originalPoints[i];
                    beta = 0.125;//1/8
                }
                else
                {
                    //interior vertex
                    //β=1/n {5/8−[3/8+1/4  cos(2π/n) ]^2 }
                    //v_old' = (1−n⋅β)⋅v_old+Σ_(j=1)^n (β⋅v_j )
                    const auto inner_bracket = (3.0 / 8.0) + (1.0 / 4.0) * cos((2.0 * std::numbers::pi) / n);
                    const auto bracket_squared = inner_bracket * inner_bracket;
                    beta = (1.0 / n) * ((5.0 / 8.0) - bracket_squared);

                    updatedPoints[i] = (1.0 - n * beta) * originalPoints[i];
                }

                for(const auto vid: neighborVids)
                {
                    updatedPoints[i] += beta * originalPoints[vid];
                }
            }

            //step 3: update topology and create new mesh
            //points
            for(const auto& point: newVertices)
            {
                updatedPoints.push_back(point);
            }

            //triangles
            std::vector<std::vector<int>> newTriangles;

            //newVertices internally maps eid->points
            //eg, eid 0->newVertices[0]
            //here we need to know the index of point in updatedPoints
            //the index would be eid+oldpoints.count
            for(int cellID = 0; cellID < currentMesh->GetNumberOfCells(); ++cellID)
            {
                const auto eids = triangleEidsVec[cellID];

                std::vector<int> newVids(3);
                for(int i = 0; i < 3; ++i)
                {
                    newVids[i] = eids[i] + pointsCnt;
                }

                //interior triangle: {N0, N1, N2}
                newTriangles.push_back({newVids[0], newVids[1], newVids[2]});

                //corner triangles: {p0, N0, N2}, {p1, N1, N0}, {p2, N2, N1}
                const auto triVids = triangleVidsVector[cellID];
                newTriangles.push_back({triVids[0], newVids[0], newVids[2]});
                newTriangles.push_back({triVids[1], newVids[1], newVids[0]});
                newTriangles.push_back({triVids[2], newVids[2], newVids[1]});
            }

            //set mesh topology
            auto vtkPointsVec = vtkSmartPointer<vtkPoints>::New();
            for(const auto& p: updatedPoints)
                vtkPointsVec->InsertNextPoint(p.GetData());

            auto cells = TopologyUtil::GetTriangleTopologyAsCellArray(newTriangles);

            auto newMesh = vtkSmartPointer<vtkPolyData>::New();
            newMesh->SetPoints(vtkPointsVec);
            newMesh->SetPolys(cells);
            newMesh->Modified();

            currentMesh->DeepCopy(newMesh);
        }

        return currentMesh;
    }
}
