#pragma once

#include <vtkSmartPointer.h>

class vtkPolyData;

namespace AlgorithmAoS
{
    vtkSmartPointer<vtkPolyData> GetLoopSubdivisionMesh(vtkPolyData* mesh, int iteration);
}
