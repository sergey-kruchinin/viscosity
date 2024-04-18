#!/bin/bash
#SBATCH --job-name=nemd
#SBATCH --error=nemd-%A-%a.err
#SBATCH --output=nemd-%A-%a.log
#SBATCH --array=1-4
#SBATCH --time=3-00:00:00
#SBATCH --ntasks 8
#SBATCH --nodes 1


module load GROMACS/2022.2

IDX=${SLURM_ARRAY_TASK_ID}
INIT_DIR=init
WORK_DIR=traj_${IDX}
ACCELERATIONS=($(seq 0.010 0.005 0.040))
LOGIN_SERVER=sms

# STAGE 1
if $(squeue -u $(whoami) -o '.18%i %.Z' | grep -q "$(pwd)/${WORK_DIR}")
then
    squeue -u $(whoami) -o 'Cannot start the job in %Z due to conflict with job %i' | grep "$(pwd)/${WORK_DIR}"
    exit 1
fi    

if [ ! -d ${WORK_DIR} ]
then
    cp -r ${INIT_DIR} ${WORK_DIR}
fi

if $(diff ${INIT_DIR} ${WORK_DIR} | grep -q ${INIT_DIR}) 
then 
    rm -rf ${WORK_DIR}
    cp -r ${INIT_DIR} ${WORK_DIR}
fi
cd ${WORK_DIR}
#
if ! $(grep -q 'Finished mdrun' em.log) 
then
    gmx_mpi grompp -f em.mdp -o em.tpr -p mix.top -c emimCG_NVT.gro -maxwarn 10
    srun --mpi=pmix_v3 gmx_mpi mdrun -ntomp 1 -v -deffnm em 
fi

if ! $(grep -q 'Finished mdrun' eq.log) 
then
    gmx_mpi grompp -f eq.mdp -o eq.tpr -p mix.top -c em.gro -maxwarn 10
    srun --mpi=pmix_v3 gmx_mpi mdrun -ntomp 1 -v -deffnm eq
fi

if ! $(grep -q 'Finished mdrun' prod.log) 
then
    gmx_mpi grompp -f prod.mdp -o prod.tpr -p mix.top -c eq.gro -maxwarn 10
    srun --mpi=pmix_v3 gmx_mpi mdrun -ntomp 1 -v -deffnm prod
fi

# STAGE 2
cat > run_nemd.sbatch << EOF
#!/bin/bash
#SBATCH --job-name=${SLURM_JOB_NAME}-${IDX}
#SBATCH --error=${SLURM_JOB_NAME}-%A-%a.err
#SBATCH --output=${SLURM_JOB_NAME}-%A-%a.log
#SBATCH --array=0-$((${#ACCELERATIONS[@]}-1))
$(grep -P '^#SBATCH --time' $0)
$(grep -P '^#SBATCH --ntasks' $0)
$(grep -P '^#SBATCH --nodes' $0)


module swap gnu8 gnu10/10.1
module load GROMACS/2022.2

IDX=\${SLURM_ARRAY_TASK_ID}
ACC=(${ACCELERATIONS[@]}) 

sed 's|cos-acceleration.*|cos-acceleration = '"\${ACC[\${IDX}]}"'|' nemd_template.mdp > nemd_"\${ACC[\${IDX}]}".mdp
gmx_mpi grompp -f nemd_"\${ACC[\${IDX}]}".mdp -o nemd_"\${ACC[\${IDX}]}".tpr -p mix.top -c prod.gro -maxwarn 10
srun --mpi=pmix_v3 gmx_mpi mdrun -ntomp 1 -v -deffnm nemd_"\${ACC[\${IDX}]}"
echo 33 | gmx_mpi energy -f nemd_"\${ACC[\${IDX}]}".edr -o inverse_viscosity_"\${ACC[\${IDX}]}".xvg | grep -Po '(?<=1/Viscosity)[\-0-9.\s]+' > inverse_viscosity_average_"\${ACC[\${IDX}]}".dat


#  1  Angle            2  LJ-(SR)          3  Coulomb-(SR)     4  Coul.-recip.
#  5  LJ-recip.        6  Potential        7  Kinetic-En.      8  Total-Energy
#  9  Conserved-En.   10  Temperature     11  Pressure        12  Constr.-rmsd
# 13  Vir-XX          14  Vir-XY          15  Vir-XZ          16  Vir-YX
# 17  Vir-YY          18  Vir-YZ          19  Vir-ZX          20  Vir-ZY
# 21  Vir-ZZ          22  Pres-XX         23  Pres-XY         24  Pres-XZ
# 25  Pres-YX         26  Pres-YY         27  Pres-YZ         28  Pres-ZX
# 29  Pres-ZY         30  Pres-ZZ         31  #Surf*SurfTen   32  2CosZ*Vel-X
# 33  1/Viscosity     34  Coul-SR:AP-AP   35  LJ-SR:AP-AP     36  Coul-SR:AP-AN
# 37  LJ-SR:AP-AN                         38  Coul-SR:AP-rest
# 39  LJ-SR:AP-rest                       40  Coul-SR:AN-AN
# 41  LJ-SR:AN-AN                         42  Coul-SR:AN-rest
# 43  LJ-SR:AN-rest                       44  Coul-SR:rest-rest
# 45  LJ-SR:rest-rest                     46  T-GRA
# 47  T-GRN           48  T-AP            49  T-AN            50  Lamb-GRA
# 51  Lamb-GRN        52  Lamb-AP         53  Lamb-AN
EOF

ssh $LOGIN_SERVER "cd $(pwd); sbatch run_nemd.sbatch"
cd ../

