# LISA Simulation

This is a simulation built to analyze the effects of asteroids in the solar system on the LISA satellites (which will launch in the 2030's). Key data files are not included (Data is accessed from jpl small body database.)

## Getting started

There are three simulations and several Python scripts. There is the Spectrum Simulation, which uses some CUDA code to run on the GPU. This simulation has not been maintained with the file structure changes and new compilation format, so it will take some care to get compiled and running again. There is the Secular Simulation, which is written entirely in C++ and is up to date. And finally there is the Bolen Test Simulation, which is a "how to" simulation for a future project. This is also up to date. The Python scripts are for handling and graphing the the outputs of these simulations.

### Prerequisites

This was built and tested on linux, though it should be easy enough to run on windows.

#### Spectral Simulation Requirements

For use on linux, ppa drivers were used. This was to get the nvidia gpu working properly. This code will also require nvcc to compile, and will require a nvidia GPU with CUDA support.

#### Secular and Bolen Test Requirements

This was coded entirely in C++, so no exceptional requirements are needed. Though, a C++17 compatable compiler will be needed, and the std flag will need to be set to C++17. For convenience, a CMake file is included. For this, CMake 3.0 is needed, as well as Make. I'm uncertain what the minimum requirement for this project is, but I have been using Make 4.1.

### Building

#### Spectral Simulation Build

In "Spectrum_Simulation_Code/Sh_Scripts" there is a suite of compiling scripts. "comp_all.sh" is something of a quickstart compile for the main simulation. You'll need a file called "Executables" and "Obects" next to "Sh_Scripts" for compilation files to go into. If you run the executable, there should be a "file not found" type error.

#### Secular and Bolen Test Build

Both of these simulations are built with one CMake script. To build these, create a build directory separate from the git directory. Once navigated to the build directory, call cmake on the git directory. It will automatically call the CMake script and build the necessary build files. Then call make, and both simulations will be compiled in the build directory. At his point they will seg fault due to not having the correct reference files. See Data Files for more information.

### Data Files

#### Reference File

Due to file size, key data files were excluded. This data comes from the jpl small body database. For simplicity, if you do not have the Reference file, send me an email and I will send you a zip file of it.

#### Output Files

The CMake script automatically creates an output file for the secular simulation as well as one for the Bolen test simulation. These are labelled clearly in the build directory and will be the destination of all outputs of the respective simulations.
