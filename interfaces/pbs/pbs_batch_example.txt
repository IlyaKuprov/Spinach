#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=59:00:00
#PBS -V
#PBS -m n
#PBS -t 1-48
#PBS -q your_queue_name

cd $PBS_O_WORKDIR

# Replace this command as appropriate
module load matlab/2016b

matlab -nodesktop -r "your_function_name(`expr ${PBS_ARRAYID}`); exit"   \
                  1 > your_function_name_`expr ${PBS_ARRAYID}`.out       \
                  2 > your_function_name_`expr ${PBS_ARRAYID}`.err       \
                    < /dev/null



