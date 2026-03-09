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

Declare a `CoCiP` object.

Assign a shared pointer to a `Params` object to `CoCiP::params` and a unique pointer to a derived `IMet` object (e.g. `ArrayMet`) to `CoCiP::met` respectively.

Assign flight inputs and initial time and location to the `CoCiP` object, as well as initial meterological values to the `CoCiP::met` object if required.

Call `CoCiP::formation`.

If `CoCiP::sac` is true, call `CoCiP::simulate_wake_vortex_downwash`, then `CoCiP::initial_properties`.

If `CoCiP::persistent` is true, call `CoCiP::process_downwash_flight` passing the segment heading angle.

In a loop, advect the plume location, update the meterological variables, time, and location, then call `CoCiP::evolve` passing the segment length ratio, the segment heading angle, and the time step duration. End the loop is `CoCiP::persistent` is not true.

### Notes

CoCiP++ is location agnostic. It does not advect the plume internally; this behaviour must be driven externally. The time and location are used solely to estimate solar direct radiation.

The plume will sediment and and alter `CoCiP::altitude`, so ensure that advection does not overwrite this alteration.

## Features

CoCiP++ is based on the [pycontrails](https://github.com/contrailcirrus/pycontrails/) implementation of CoCiP. It is up to date with pycontrails v0.60.3.

Currently, CoCiP++ is missing the following features:
- Extended K15 model
- Humidity scaling
- Contrail contrail overlap

## Acknowledgements

This project has been developed by Jack Bartlett.

The original CoCiP model was developed by Ulrich Schumann and detailed in [Schumann (2012)](https://doi.org/10.5194/gmd-5-543-2012).

Acknowledgement for the transliterated code is directed to [pycontrails](https://github.com/contrailcirrus/pycontrails/).