#pragma once

#include <vtkSmartPointer.h>

class vtkPolyData;

namespace Algorithm
{
    vtkSmartPointer<vtkPolyData> GetLoopSubdivisionMesh(vtkPolyData* mesh, int iteration);
}
