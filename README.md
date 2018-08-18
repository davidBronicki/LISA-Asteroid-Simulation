# LISA Simulation

This is a simulation built to analyze the effects of asteroids in the solar system on the LISA satellites (which will launch in the 2030's). Key data files are not included, I'll put in some details on how to get the correct data in file later. (Data is accessed from jpl small body database.)

## Getting started

This is a guide for getting the key object files and executable files up.

### Prerequisites

This was built and tested on linux, though it should be easy enough to run on windows. Also, ppa drivers were used. There's only one place that provides these for linux, so if you find something you've got the right stuff. You'll need nvidia's cuda tools to compile this aswell as a nvidia gpu. If typing 'nvcc' in terminal returns an error, you've got a problem.

### Building

 In "LISA_Simulation_Code" there is a "Sh_Scripts" file which contains the core of the compiling scripts. "comp_all.sh" is something of a quickstart compile for the main simulation. You'll need a file called "Executables" and "Obects" next to "Sh_Scripts" for compilation files to go into. If you run the executable, there should be a "file not found" type error.

### Data files

Due to file size, key data files were excluded. This data comes from the jpl small body database. In a future release, I'll put in the full file structure as well as instructions for how to fill the files.
