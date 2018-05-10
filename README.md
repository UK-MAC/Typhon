# Typhon

Typhon is a distributed communications library for unstructured mesh
applications.

Typhon is implemented in C++, but the API is additionally available in Fortran
and C.

## Installation through Spack

Typhon can be installed using [Spack](https://spack.io/).

```
spack install typhon
```

## Manual installation

### Dependencies

The only hard dependency is MPI.

[yaml-cpp](https://github.com/jbeder/yaml-cpp) is an optional dependency, used
for dumping serialised internal Typhon data structures using the `typhon_dump`
tool (built when `BUILD_TYPHON_DUMP` is specified). If you plan to develop
software using Typhon, you may find this functionality useful.

If Doxygen is present and `BUILD_DOCS` is specified, Typhon will build browsable
HTML documentation. This can be found at `html/index.html` in the build
directory.

### Building and installing

Typhon uses CMake, and tries to be idiomatic in doing so. The typical CMake
process looks something like the following (inside the top-level Typhon
directory):

```
mkdir build
cd build
cmake \
    -DCMAKE_INSTALL_PREFIX=$HOME \
    -DCMAKE_BUILD_TYPE="Release" \
    ..
make
make test
make install
```

If you are having trouble getting CMake to find yaml-cpp, you can pass it the
`YamlCpp_ROOT` variable specifying its location. For issues with MPI, please
refer to the CMake [FindMPI documentation](https://cmake.org/cmake/help/latest/module/FindMPI.html).

## Usage

To use Typhon from C, C++ or Fortran, simply `#include <typhon.h>` (in C/C++) or
`use typhon` (in Fortran) and link against `libtyphon`. Remember that Typhon
requires MPI in order to function.

To simplify usage, the C and C++ APIs for Typhon are identical. This means that
function overloading and default arguments are both not supported, as these are
specifically C++ features.

The Fortran API is implemented as a set of wrapper routines around bindings to
the C API. These wrapper routines handle marshalling between C and Fortran
types, and are included in the `libtyphon` binary, so no additional linking is
required.

*Note on boolean data: As the bool type in C++ has some unhelpful qualities, and
as there does not exist a clean way of interfacing between C and C++ code that
uses _Bool/bool, Typhon uses 'int' as the underlying type for all boolean
operations (mainly TYPH_Gather and TYPH_Reduce). The expectation is that '1' is
used to represent true, and '0' to represent false. The Fortran wrappers
correspondingly expect 'logical(kind=int32)', and the same true/false values
should be used.*

### Examples

For examples of using Typhon, please see the `examples` directory.

## Release history

* 02/05/2018, version 3.0
