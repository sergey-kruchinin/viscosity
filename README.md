# run_green_cubo.sbatch
**Usage.**
Run as `sbatch run_green_cubo.sbatch`. It is a simple script for SLURM for simultaneous submitting of several tasks with common parameters (job array). It submits a given number of MD simulations for creating a hierarchical structure of directories and files needed for running script `viscosity.py`. The script needs a start directory with necessary files for MD calculations. This directory will be copied by the script the required number of times and a job with given commands will be submitted for every created directory.

**Setup**
*	The start directory name is set in the script by variable `INIT_DIR` (default is `init` in the same directory with the script);
*	Number of copies is set as SBATCH variable `--array=N-M` (N – minimal job index, М – maximal job index);
*	The directories names for submitting jobs are constructed from two parts: first - name, which is set in the script by variable `WORK_DIR` (default is `traj_`) and second - index, which is taken from values of SBATCH variable `--array` (see before);
*	Desireable amount of CPU, GPU and another computational resources is set by common way for a SBATCH scripts. The indicated amount of computational resources is the amount for the one job;
*	The command sequence for running is written in the end of the script.

After successful script execution the (M-N) directories with the calculations data will be created. 

# viscosity.py
**Usage.**
For execution the script needs Python 3 (was tested on v.3.8) interpretator and NumPy (was tested on v.1.21.5) and SciPy(was tested on v.1.7.3) libraries. The script calculates the shear viscosity rate from the set of autocorrelation functions of the pressure tensor components by the algorithm from the paper [^Zhang2015]. For correct work the hierarchical structure of the input data must be created by the script  `run_green_cubo.sbatch` or analogous. 

`viscosity.py [-h] [--fraction FRACTION] [--amount AMOUNT] [--oa AVE_NAME] [--os STD_NAME] [--of FIT_NAME] [directory] [stem]`
* `directory` — The name of directory with the results of calculations. Default is current directory;
* `stem`  — The common name part for directories with the results of calculations. Default is `traj_`;

**Parameters:**
* `--fraction N` — Ratio of STD value to its average. Using for calculation t<sub>cut</sub>. Default is 0.4;
*	`--amount AMOUNT` — Number of directories using for analysis. If lesser than 0 or greater than number of available directories, all directories will be taken. Default is 0;
* `--oa AVE_NAME` — File name for saving the averaged viscosity time dependence [ps, P]. If not set this data will not be saved;
*	`--os STD_NAME` — File name for saving the viscosity STD time dependence [ps, P]. If not set this data will not be saved;
* `--of FIT_NAME` — File name for saving the fitted viscosity time dependence [ps, P]. If not set this data will not be saved;
* `-h`, `--help` — Usage help.

After successful execution the script will print the shear viscosity rate in poises [P]. 
 
# run_nemd.sbatch
**Usage.**
Run as `sbatch run_nemd.sbatch`. It is a script for SLURM for simultaneous submitting of several jobs with common parameters (job array). It submits a given number of MD simulations for calculating shear viscosity rate by spatially oscillating perturbations method. At first stage the script submits N simulations with the same modelling parameters. At second stage the script creates additional `.sbatch` files which submits a simulations with different values of `cos-acceleration` parameter (acceleration of system particles) [^AllenCSL] for every simulation from first stage. The script needs a directory with necessary files for running MD simulations. This directory will be copied by the script the required number of times and a job with given commands will be submitted for every created directory.

**Setup**
* The input directory is set in the script by the variable `INIT_DIR` (Default is `init` in the same directory with the script);
* Number of copies for first stage is set as SBATCH variable `--array=N-M` (N – minimal job index, М – maximal job index);
*	The directories names for the running tasks are constructed from two parts: first - name, which is set in the script by variable `WORK_DIR` (default is `traj_`) and second - index, which is taken from values of SBATCH variable `--array` (see before);
* The `cos-acceleration` values is set in the variable `ACCELERATIONS`;
* The name of the server which used for submitting jobs is set in the variable `LOGIN_SERVER`. It is needed for running second stage of calculation;
* The desireable amount of CPU, GPU and another computational resources is set by common way for a SBATCH scripts. The indicated amount of computational resources is the amount for the one job.

After successful script execution the (M-N) directories with the files containing shear rate viscosities for every given acceleration will be created.

[^Zhang2015]: Y. Zhang, A. Otani, E.J. Maginn, Reliable Viscosity Calculation from Equilibrium Molecular Dynamics Simulations: A Time Decomposition Method, J. Chem. Theory Comput. 11 (2015) 3537–3546. https://doi.org/10.1021/acs.jctc.5b00351.

[^AllenCSL]: M.P. Allen, D.J. Tildesley, Computer Simulation of Liquids, (2017). https://doi.org/10.1093/oso/9780198803195.001.0001.

