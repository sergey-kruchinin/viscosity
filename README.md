# run_green_cubo.sbatch
**Usage.**
Run as `sbatch run_green_cubo.sbatch`. It is a simple SLURM script, used to simultaneously submit several tasks with common parameters (job array). It submits a given number of MD simulations and creates a hierarchical structure of directories and files needed for subsequent run of the script `viscosity.py` in a current directory (project directory). The script needs an input directory with necessary files for MD calculations. It will copy this directory required number of times and submit the job with given commands for each created directory.

**Setup**
*	The input directory name is set in the script by variable `INIT_DIR` (by default it is `init` in the same directory with the script);
*	Number of copies is set by SBATCH variable `--array=N-M` (N – minimum job index, М – maximum job index);
*	The directories names of submitted jobs are constructed from two parts: first - the root name, which is set in the script by variable `WORK_DIR` (by default it is `traj_`), and second - index, which is taken from the values of SBATCH variable `--array` (see before);
*	Desireable amount of CPU, GPU and other computational resources is set in a way, standard for the SBATCH scripts. The indicated amount of computational resources is the amount for single job;
*	The command sequence for submittion is written below the setup part.

After successful script execution the output data will be organized in these directories for subsequent usage of the `viscosity.py` script. 

# viscosity.py
**Usage.**
For execution the script needs Python 3 (was tested on v.3.8) interpretator and NumPy (was tested on v.1.21.5) and SciPy (was tested on v.1.7.3) libraries. The script calculates the shear viscosity rate from the set of autocorrelation functions of the pressure tensor components by the algorithm from the paper [^Zhang2015]. For correct work the input data must be organized by the script `run_green_cubo.sbatch` (as written above) or by analogous one. 

`viscosity.py [-h] [--fraction FRACTION] [--amount AMOUNT] [--oa AVE_NAME] [--os STD_NAME] [--of FIT_NAME] [directory] [stem]`
* `directory` — The name of project directory. Default is current directory;
* `stem`  — The common name part for directories with the results of calculations. Default is `traj_`;

**Parameters:**
* `--fraction N` — Threshold value for STD to average ratio. Used for calculation of t<sub>cut</sub>. Default is 0.4;
*	`--amount AMOUNT` — Number of directories used for analysis. If lesser than 0 or greater than number of available directories, all directories will be used. Default is 0;
* `--oa AVE_NAME` — File name for saving the averaged viscosity time dependence [ps, P]. If not set this data will not be saved;
*	`--os STD_NAME` — File name for saving the viscosity STD time dependence [ps, P]. If not set this data will not be saved;
* `--of FIT_NAME` — File name for saving the fitted viscosity time dependence [ps, P]. If not set this data will not be saved;
* `-h`, `--help` — shows this help message and exit.

After successful execution the script will print the shear viscosity rate in poises [P]. 
 
# run_nemd.sbatch
**Usage.**
Run as `sbatch run_nemd.sbatch`. It is a SLURM script, used to simultaneously submit several tasks with common parameters (job array). It submits a given number of MD simulations for calculating shear viscosity rate by spatially oscillating perturbations method. At first stage the script submits N simulations with the same modelling parameters. At second stage the script creates additional `.sbatch` files which submit simulations with different values of `cos-acceleration` parameter (acceleration of system particles) [^AllenCSL] for each simulation from the first stage. The script needs a directory with necessary files for running MD simulations. It will copy this directory required number of times and submit the job with given commands for each created directory.

**Setup**
* The input directory is set in the script by the variable `INIT_DIR` (Default is `init` in the same directory with the script);
* Number of copies for the first stage is set by the SBATCH variable `--array=N-M` (N – minimum job index, М – maximum job index);
* The directories names of submitted jobs are constructed from two parts: the root name, which is set in the script by variable `WORK_DIR` (by default it is `traj_`), and second - index, which is taken from the values of SBATCH variable `--array` (see before);
* The `cos-acceleration` values is set in the variable `ACCELERATIONS`;
* The name of the server, which is used to submit the jobs, is set by the variable `LOGIN_SERVER`. It is needed for execution of the second stage of calculation;
* Desireable amount of CPU, GPU and other computational resources is set in a way, standard for the SBATCH scripts. The indicated amount of computational resources is the amount for single job.

After successful script execution the files containing shear rate viscosities for each given acceleration value will be placed in corresponding directories. 

[^Zhang2015]: Y. Zhang, A. Otani, E.J. Maginn, Reliable Viscosity Calculation from Equilibrium Molecular Dynamics Simulations: A Time Decomposition Method, J. Chem. Theory Comput. 11 (2015) 3537–3546. https://doi.org/10.1021/acs.jctc.5b00351.

[^AllenCSL]: M.P. Allen, D.J. Tildesley, Computer Simulation of Liquids, (2017). https://doi.org/10.1093/oso/9780198803195.001.0001.

