# Comet
**Co**mputational **m**odeling of **e**volving **t**umors

Comet is built on [Tumopp](https://github.com/heavywatal/tumopp), a single-cell based, flexible and efficient spatial tumor growth simulator written in C++.

The differences between Comet and the original tumopp packages are:
- [x] passenger mutations are generated and printed out for each cell
- [x] all dead cells and lineages are printed out (for lineage tracing)
- [x] Metastatic seeding cells are sampled randomly over time (removed from the growing tumor)

Test

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
comet -h
comet -N1000000 -D3 -Cmoore -k1 -u 0.15 -b 0.25 -d 0.18 --ub 0.0005 --mb 0.05 --extinction 100000 --local linear --seedingSize 5000 -o $OUTPUT_DIR
```

## Reference:
Ruping Sun and Athanasios N. Nikolakopoulos\
**Elements and Evolutionary Determinants of Genomic Divergence Between Paired Primary and Metastatic Tumors**,\
__bioRxiv 2020.08.24.262378__,\
☕**DOI:** [10.1101/2020.08.24.262378](https://doi.org/10.1101/2020.08.24.262378)

## Author:
Ruping Sun, Athanasios Nikolakopoulos
