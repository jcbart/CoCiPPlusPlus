# CoCiP++

CoCiP++ is a basic C++11 implementation of the Contrail Cirrus Prediction tool (CoCiP).

It is designed to be used both as a standalone program and as a library (in the future).

## Compilation

### System requirements
- C++ compiler supporting C++11 or later
- CMake
- Git

### Bundled dependencies

The following are included as Git submodules and are compiled with the project:

- [yaml-cpp](https://github.com/jbeder/yaml-cpp)

### Optional external dependencies
- [Earth System Modeling Framework (ESMF)](https://earthsystemmodeling.org/) (in the future)

### Instructions

Clone the GitHub repo with
```bash
git clone https://github.com/jcbart/CoCiPPlusPlus.git
```

In the top-level directory, run the following commands:
```bash
git submodule update --init --recursive
mkdir build
cd build
cmake ..
cmake --build .
```

This will compile the project in `CoCiPPlusPlus/build`. The project can be compiled in any directory by instead running `cmake /path/to/CoCiPPlusPlus`.

The project can be compiled in debug mode by specifying `cmake .. -DCMAKE_BUILD_TYPE=Debug`. The default build type is `Release`.

Successful compilation will create an executable named `CoCiP`.

If they do not already exist, CMake will copy the default input files `CoCiP-config.yaml` and `CoCiP-params.yaml` to the build directory.

## Execution

Runtime configuration options can be specified in `CoCiP-config.yaml`. Parameters can be specified in `CoCiP-params.yaml`. These files must be in the same directory as the `CoCiP` executable.

Run the program with the command `./CoCiP`.

Outputs are (currently) written to `cocip.out`.

## Features

CoCiP++ is based on the [pycontrails](https://github.com/contrailcirrus/pycontrails/) implementation of CoCiP.

Currently, CoCiP++ is missing the following features:
- Extended K15 model
- Humidity scaling
- Contrail contrail overlap

## Acknowledgements

This project has been developed by Jack Bartlett.

The original CoCiP model was developed by Ulrich Schumann and detailed in [Schumann (2012)](https://doi.org/10.5194/gmd-5-543-2012).

Acknowledgement for the transliterated code is directed to [pycontrails](https://github.com/contrailcirrus/pycontrails/).