# Comet
**Co**mputational **m**odeling of **e**volving **t**umors

this is a modification of tumopp, a single-cell based spatial tumor growth simulation framework.
[Project page on GitHub](https://github.com/heavywatal/tumopp)

The differences between Comet and the original tumopp packages is:
- [x] passenger mutations are generated and printed out for each cell
- [x] all dead cells and lineages are printed out (for lineage tracing)
- [x] Metastatic seeding cells are sampled randomly over time (removed from the growing tumor)

*tumopp* is a fast and flexible tumor growth simulator.
The core library is written in C++.

## Requirements

- Unix-like environment (macOS, Linux, WSL, MinGW on MSYS2, etc.)
- C++14 compiler (clang++ >= Apple LLVM 8.1, g++ >= 5.3)
- [CMake](https://cmake.org/) (>= 3.8.0)

The following libraries are optional or automatically installed:

- [clippson](https://github.com/heavywatal/clippson)
- [cxxwtl](https://github.com/heavywatal/cxxwtl)
- [sfmt-class](https://github.com/heavywatal/sfmt-class)
- [zlib](https://zlib.net)


## Installation

```sh
git clone https://github.umn.edu/sunpath/Comet.git
cd Comet/
mkdir build
cd build/
DESTINATION=${HOME}/local
cmake -DCMAKE_INSTALL_PREFIX=$DESTINATION ..
make -j2
make install
PATH=${DESTINATION}/bin:$PATH
```

Example:
```sh
tumopp -h
tumopp -N20000 -D3 -Chex -k100 -d0.1 -m0.5 -o OUTPUT_DIR
```

## Reference:
to be added

## Author:
Ruping Sun