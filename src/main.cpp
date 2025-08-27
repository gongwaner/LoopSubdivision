#include <filesystem>

#include <vtkPolyData.h>

#include "IOUtil.h"

#include "AlgorithmHelper.h"
#include "LoopSubdivision.h"


int main()
{
    const auto dataDir = AlgorithmHelper::GetDataDir();

    const auto file = dataDir / "bunny.stl";
    auto mesh = IOUtil::ReadMesh(file.string());

    const int iterationCount = 2;

    const auto resultMesh = Algorithm::GetLoopSubdivisionMesh(mesh, iterationCount);

    const auto fileDir = dataDir / "result.stl";
    IOUtil::WriteMesh(fileDir.string(), resultMesh);


    return 0;
}
