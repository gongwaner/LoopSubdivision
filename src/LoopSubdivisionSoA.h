#pragma once

#include <vtkSmartPointer.h>

class vtkPolyData;

namespace AlgorithmSoA
{
    vtkSmartPointer<vtkPolyData> GetLoopSubdivisionMesh(vtkPolyData* mesh, int iteration);
}
