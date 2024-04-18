#!/bin/bash
#SBATCH --job-name=viscosity
#SBATCH --error=viscosity-%A-%a.err
#SBATCH --output=viscosity-%A-%a.log
#SBATCH --array=1-40
#SBATCH --time=16-00:10:00
#SBATCH --ntasks 1
#SBATCH --nodes 1


module load GROMACS/2022.2 

IDX=${SLURM_ARRAY_TASK_ID}
INIT_DIR=init
WORK_DIR=traj_${IDX}

rm -rf ${WORK_DIR}
cp -r ${INIT_DIR} ${WORK_DIR}
cd ${WORK_DIR}

gmx_mpi grompp -f em.mdp -o em.tpr -p mix.top -c BMIMBF4.gro
srun --mpi=pmix_v3 gmx_mpi mdrun -deffnm em -ntomp 1 -v 

gmx_mpi grompp -f eq.mdp -o eq.tpr -p mix.top -c em.gro
srun --mpi=pmix_v3 gmx_mpi mdrun -deffnm eq -ntomp 1 -v 

gmx_mpi grompp -f prod.mdp -o prod.tpr -p mix.top -c eq.gro
srun --mpi=pmix_v3 gmx_mpi mdrun -deffnm prod -ntomp 1 -v 
echo 25 26 29 0 | gmx_mpi energy -f prod.edr -s prod.tpr -o press_XY_XZ_YZ.xvg -b 1000
gmx_mpi analyze -f press_XY_XZ_YZ.xvg -ac acf.xvg -subav no -normalize no -oneacf yes
echo 12 0 | gmx_mpi energy -f prod.edr -s prod.tpr -b 1000 -aver yes | grep -Po '(?<=Temperature)[\-0-9.\s]+' > temperature.dat
gmx_mpi editconf -f prod.gro 2>/dev/null | grep -Po '(?<=Volume:)[0-9.\s]*(?=nm\^3)' > volume.dat

cd ../
