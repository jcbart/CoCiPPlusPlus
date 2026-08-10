# CoCiP++

CoCiP++ is a work-in-progress C++20 implementation of the Contrail Cirrus Prediction tool (CoCiP).

It is designed to be used both as a standalone program and as a library.

## Compilation

### System requirements
- C++ compiler supporting C++20
- CMake
- Git

### Bundled dependencies

The following are included as Git submodules and are compiled with the project:

- [yaml-cpp](https://github.com/jbeder/yaml-cpp)

### Instructions

Clone the GitHub repo with
```bash
git clone https://github.com/jcbart/CoCiPPlusPlus.git
```
and, in the top-level directory, run the following command to gather the Git submodules:
```bash
git submodule update --init --recursive
```

Use the following commands to configure the project in a directory named `build`:
```bash
mkdir build
cd build
cmake ..
```
The project can be compiled in any directory by instead running `cmake /path/to/CoCiPPlusPlus`.

Configuration options that can be used in the `cmake` command are listed below:
- `-DCMAKE_BUILD_TYPE=Debug`: Compile with debug flags and no optimisation (default is `-DCMAKE_BUILD_TYPE=Release`).
- `-DCOCIP_BUILD_EXECUTABLE=OFF`: Do not build an executable; only a static library will be built (default is `-DCOCIP_BUILD_EXECUTABLE=ON`).
- `-DCOCIP_USE_EXTERNAL_LIBRARIES=ON`: CMake will search for required libraries on your system instead of building the Git submodules as internal targets (default is `-DCOCIP_USE_EXTERNAL_LIBRARIES=OFF`).

Finally, run
```bash
cmake --build .
```
to compile and build.

Successful compilation will create a static library named `libcocip.a` and, unless turned off, an executable named `CoCiP`.

If they do not already exist, CMake will copy the default input files `CoCiP-config.yaml` and `CoCiP-params.yaml` to the build directory.

## Standalone execution

Runtime configuration options can be specified in `CoCiP-config.yaml`. Parameters can be specified in `CoCiP-params.yaml`. These files must be in the current working directory. Alternatively, a path to a YAML file containing configuration options can be specified in the command line.

Run the program with the command `./CoCiP /optional/path/to/config`.

Outputs are (currently) written to `cocip.out`.

## Use in another project

The static library can be used to include CoCiP++ in another project.

### Instructions

Include the headers in `include/CoCiP++`.

Declare a `CoCiP<MetType>` object where `MetType` is one of:
1. `SimpleMet`
2. `ArrayMet<arrayType>` where `arrayType` is `float` or `double`

Assign a shared pointer to a `Params` object to `CoCiP::params` and a unique pointer to a derived `IMet` object (e.g. `ArrayMet`) to `CoCiP::met` respectively.

Assign flight inputs and initial time and location to the `CoCiP` object, as well as initial meterological values to the `CoCiP::met` object if required.

Call `CoCiP::formation`.

If `CoCiP::sac` is true, call `CoCiP::simulate_wake_vortex_downwash`, then `CoCiP::initial_properties`.

If `CoCiP::persistent` is true, call `CoCiP::process_downwash_flight` passing the segment heading angle.

In a loop, advect the plume location, update the meterological variables, time, and location, then call `CoCiP::evolve` passing the segment length ratio, the segment heading angle, and the time step duration. End the loop if `CoCiP::persistent` is not true.

### Notes

CoCiP++ is location agnostic. It does not advect the plume internally; this behaviour must be driven externally. The time and location are used solely to estimate solar direct radiation.

The plume will sediment and alter `CoCiP::altitude` during `evolve`, so ensure that advection takes this alteration into account.

The parameters in `CoCiP-params.yaml` generally reflect those in pycontrails' [`cocip_params.py`](https://github.com/contrailcirrus/pycontrails/blob/main/pycontrails/models/cocip/cocip_params.py) except for:
- `interp_with_pressure`: if true, will interpolate temperature, specific humidity, and wind linearly against linear pressure (matching pycontrails) and, if false, interpolate against linear geopotential height; and
- `effective_vertical_resolution`: may be null in CoCiP++ in which case the actual vertical resolution of the met data will be used (only applicable with `ArrayMet`).

## Features

CoCiP++ is based on the [pycontrails](https://github.com/contrailcirrus/pycontrails/) implementation of CoCiP. It is up to date with pycontrails v0.63.1.

Currently, CoCiP++ is missing the following features:
- Extended K15 model
- Humidity scaling
- Contrail contrail overlap
- Schumann 2025 radiative forcing model

## Acknowledgements

This project has been developed by Jack Bartlett.

The original CoCiP model was developed by Ulrich Schumann and detailed in [Schumann (2012)](https://doi.org/10.5194/gmd-5-543-2012).

Most code in CoCiP++ has been translitered from [pycontrails](https://github.com/contrailcirrus/pycontrails/). Credit for this code is directed to them.