#include <filesystem>

#include <vtkPolyData.h>

#include "IOUtil.h"

#include "AlgorithmHelper.h"
#include "LoopSubdivisionAoS.h"
#include "LoopSubdivisionSoA.h"


int main()
{
    const auto dataDir = AlgorithmHelper::GetDataDir();

    const auto file = dataDir / "bunny.stl";
    auto mesh = IOUtil::ReadMesh(file.string());

    const int iterationCount = 2;

    {
        auto start = std::chrono::high_resolution_clock::now();

        const auto resultMesh = AlgorithmAoS::GetLoopSubdivisionMesh(mesh, iterationCount);

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "GetLoopSubdivisionMesh() AoS version takes " << duration << " ms" << std::endl;

        const auto fileDir = dataDir / "resultAoS.stl";
        IOUtil::WriteMesh(fileDir.string(), resultMesh);
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        const auto resultMesh = AlgorithmSoA::GetLoopSubdivisionMesh(mesh, iterationCount);

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "GetLoopSubdivisionMesh() SoA version takes " << duration << " ms" << std::endl;

        const auto fileDir = dataDir / "resultSoA.stl";
        IOUtil::WriteMesh(fileDir.string(), resultMesh);
    }

    return 0;
}
