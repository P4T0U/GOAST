<div class="panel panel-success"  align="center">
<div class="panel-body">

<img src="doc/GOAST.svg"  width="500" height="500">

</div>

**goast**, *n.* &nbsp;&nbsp;  **1.** An obsolet spelling of ghost.[^1] &ndash; **2.** A C++ library for Geometry Processing.

</div>


# The Geometric Optimization and Simulation Toolbox (GOAST)
The Geometric Optimization and Simulation Toolbox (short: GOAST) is a C++ library for research in Geometry Processing. 
It currently focuses on variational approaches and triangle meshes, especially methods for the Riemannian shape space of discrete shells. 

**Warning:** GOAST is still in a fairly early stage of development, and a lot of the code is not up to the standards we want to set for our selves yet. 
That means that, in the coming months, the library will likely undergo rapid changes. 
The included examples, however, should remain functional at all times.


[[_TOC_]]

## Dependencies
The GOAST has currently three dependencies:
 - [Eigen](http://eigen.tuxfamily.org) for everything related to linear algebra
 - [OpenMesh](https://www.openmesh.org/) as the underlying library for (triangle) meshes
 - [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html)

Furthermore, some parts are only usable if the following additional libaries are available
 - [Ipopt](https://coin-or.github.io/Ipopt/)
 - [MOSEK](https://www.mosek.com/)
 - [ADOL-C](https://projects.coin-or.org/ADOL-C)
 - [trlib](https://github.com/felixlen/trlib)
 - [VTK](https://vtk.org/)

They are enabled via CMake or preprocessor defines and we provide convenient CMake modules for their configuration, see below for more information.

## Using GOAST
Essentially, there are two different usage scenarios for GOAST.

### Standalone
If you just want to check out the examples, use the included tools, test the functionalities of GOAST, or work on the source code you can use GOAST, i.e. this repository, standalone.
To this end, you can use the common CMake routine
```bash
mkdir build
cd build
cmake ..
make
```
which will build the tests, examples, and tools by default.
Depending on how you installed the dependencies, this might not find all of them. 
In this case, you will need to tell CMake where they are explicitly by setting the options documented below.

### As dependency
GOAST is a header-only library, which means that to use it as an external dependency in your project you simply have to add GOAST's `include` directory to your include path and are ready to go.
However, in this case, you have to take care of the dependencies yourself and make sure that their headers are available and that they are linked when necessary.

If you use CMake, we provide a convenient module to directly include GOAST from its source folder.
It is located in `cmake` and includes the CMake scripts for finding dependencies provided with GOAST.
With it, you can include GOAST in the following way:
```cmake
list(APPEND CMAKE_MODULE_PATH <PATH_TO_GOAST>/cmake) # path to GOAST module
include(GOAST) # load GOAST module
```

Alternatively you can build and install GOAST to your system. For example, on a unix flavoured os
```bash
mkdir build 
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=<where to install to>
make
(sudo) make install
```
and use the CMake module provided in the directory where GOAST was installed.

In both cases, adding GOAST (and all its transitive dependancies) as a depency to your target works in both cases the same. Simply add the following line(s) to your CMakeList.txt 
```cmake
# for only using the main library:
target_link_libraries(<YOUR_TARGET> PUBLIC GOAST::GOAST)
# or if you also want to use all available additional modules/dependencies:
target_link_libraries(<YOUR_TARGET> PUBLIC GOAST::All)
```

Details on the involved CMake options for finding and using dependencies are provided below.

We also created an example project which demonstrates how you can use GOAST as a git submodule in your CMake project, you can find it [here](https://gitlab.com/numod/goast-example-project).

## Configuration
GOAST is configured via a range of CMake options and/or preprocessor defines. 
Most of them are related to external dependencies.

#### Turning optional external libraries on/off
These are both CMake options (set to `ON` or `OFF`) and preprocessor defines (simply created or not), typically named as `USE_<LIBRARY>`. 
The CMake options will enable searching for the library, which automatically stops the CMake creation if the library is not found.
If it is found, the necessary directories are added to the include paths, a CMake variable `<LIBRARY>_LIBRARIES` to link the library is defined, and the aforementioned preprocessor define is set.

| Option | Explanation | Default value |
| ---      |  ------  |---------:|
| `GOAST_WITH_IPOPT`   | Enable Ipopt as additional optimization method  | `OFF`   |
| `GOAST_WITH_ADOLC`   | Enable ADOL-C to allow automatially differentiating energies  | `OFF`   |
| `GOAST_WITH_VTK`  | Enable VTK for additional I/O functions  | `OFF`   |
| `GOAST_WITH_TRLIB`  | Enable trlib as additional trust-region subproblem solver  | `OFF`   |
| `GOAST_WITH_MOSEK`  | Enable MOSEK for solving quadratic optimization problems  | `OFF`   |
| `GOAST_WITH_MKL`  | Use Intel's MKL[^2]  | `OFF`   |

#### Search paths
These are only CMake options controlling where the search modules will look for dependencies additionally to some default paths.

| Option | Explanation  |
| ---      |  ------ |
| `SUITESPARSE_HOME`   | Directory containing the `include` and `lib` directories of SuiteSparse  |
| `EIGEN_INCLUDE_DIR`   | Directory containing Eigen includes  | 
| `OPENMESH_LIBRARY_DIR`  | Directory of your OpenMesh installation  | 
| `IPOPT_HOME`   | Directory containing the `include` and `lib` directories of Ipopt  |
| `VTK_DIR`   | Directory of the CMake modules of your VTK installation[^3]  |
| `TRLIB_ROOT`   | Directory containing the `include` and `lib` directories of trlib  |
| `MOSEK_HOME`   | Directory containing the `bin` and `h` directories of MOSEK  |

<!---
#### Standalone project
When you configure the main GOAST project, i.e. this repository, there are the following additional options:

| Option | Explanation | Default value |
| ---      |  ------  |---------:|
| `OPTIMIZE_FOR_NATIVE`   | Set `-march=native` if it is available , i.e. optimization for your specific processor architecture. | `OFF`   |
--->

## Authors
The library is primarily developed by 
[Behrend Heeren](https://ins.uni-bonn.de/staff/heeren) and 
[Josua Sassen](https://ins.uni-bonn.de/staff/sassen) in the group of 
[Martin Rumpf](https://ins.uni-bonn.de/staff/rumpf) at the 
Institute for Numerical Simulation at the University of Bonn. Other members of 
the group also contributed smaller parts of the code which is documented in the
respective files.

## License
GOAST is under the [Mozilla Public License 2.0](https://www.mozilla.org/en-US/MPL/2.0/),
see also its [FAQ](https://www.mozilla.org/en-US/MPL/2.0/FAQ/).

## Citation
If you use GOAST in a project, please cite the implemented papers as appropriate. 
To cite the library in general, please use this BibTeX entry:
```
@misc{GOAST,
	title = {{T}he {G}eometric {O}ptimization {A}nd {S}imulation {T}oolbox},
	author = {Heeren, Behrend and Sassen, Josua and others},
	url = {https://gitlab.com/numod/goast},
	year = {2020},
}
```

---------
[^1]: According to The Century Dictionary.
[^2]: To use the MKL, the environment variable `MKLROOT` needs to be set during configuration and building.
[^3]: We do not provide an own search script for VTK hence rely on the one provided by your system, which means you have to check the VTK documentation for more details.
