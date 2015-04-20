#!/bin/bash
#SBATCH --job-name=ticket18203
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --ntasks-per-node=12
#SBATCH --exclusive
#SBATCH --time=00:15:00
#SBATCH --output=DORA_%j.txt

# DEFINE SIMULATION NUMBER
nTasks=12
nTasksPerNode=12

# EXTRA SETTINGS
ulimit -s unlimited
ulimit -c unlimited

# update counters and launch
for i in {1..10}
do
    echo ' ****** START OF JOB ******'
    date
    aprun -N $nTasksPerNode -n $nTasks -cc none ./program.exe
    date
    echo ' ****** END OF JOB ****** '
done
