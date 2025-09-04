# Loop Subdivision

<img width="2029" height="935" alt="image" src="https://github.com/user-attachments/assets/40bf3a05-7736-451a-968e-06b010ce2e46" />
Implementation of loop subdivision algorithm in *[Smooth Subdivision Surfaces Based on Triangles](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/thesis-10.pdf)* by *Charles T. Loop* (2003)


## Prerequisite
VTK 9.4.0+  
CMake 3.12+


## Implementation
VTK has its own implementation of loop subdivision, see [vtkLoopSubdivisionFilter](https://vtk.org/doc/nightly/html/classvtkLoopSubdivisionFilter.html).
I re-implemented the algorithm using data-oriented design.

### Input
The program accepts a triangular mesh file as input. The following file formats are supported:
- .stl
- .ply
- .obj
- .vtk

### Output
Loop subdivision result mesh file, exported to result.stl in **data** folder


## Build
1.Clone the repo
2.Initialize the submodule  
3.Build the project
```
cd LoopSubdivision/build
cmake -DVTK_DIR=/dir/to/your/vtk/install ..
make
```
