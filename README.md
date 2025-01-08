This repository contains code to run and reproduce the simulations described in the following article:

Amelia W. Paine, Michael F. Hagan, Vinothan N. Manoharan, Disassembly of virus-like particles and the stabilizing role of nucleic acid cargo, _Journal of Physical Chemistry B_, 2025.

To regenerate the simulation files:

1. Set up the conda environments found in the "environments" folder. The hoomd-blue environment also has HOOMD version 2.9 installed.
2. In the "simulation_init" folder, run mk_0.py, mk_1.py etc., which will generate the simulation files and folders for each of the 9 initial conditions provided, with the correct random seeds. These scripts submit the simulations as Slurm jobs.

To regenerate the basic analysis of inter-subunit bonds and polymer reinforcements:

1. Copy the folder "analysis" and its contents into each of the simulation batch folders (named "0_PS_off", "1_PS_off" etc.).
2. Within the analysis folders, run mk_analysis.py, which will set up and run the analysis scripts within the folder for each simulation run in the batch. These analyses are also set up to run as Slurm jobs.

To perform higher-level analyses and recreate the plots found in the paper, use the two Jupyter notebooks found in the "analysis_notebooks" folder.
