# Typhon

Typhon is a distributed communications library for unstructured mesh
applications.

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

## Release history

* 02/05/2018, version 3.0
