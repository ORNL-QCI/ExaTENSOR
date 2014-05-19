#It is crucial to launch MPI processes consecutively within a node!
#Environment variable QF_PROCS_PER_NODE must be set appropriately!
export OMP_NUM_THREADS=2       #mandatory
export QF_NUM_PROCS=1          #mandatory
export QF_PROCS_PER_NODE=1     #mandatory
export QF_GPUS_PER_PROCESS=1   #optional
export QF_MICS_PER_PROCESS=0   #optional
export MIC_ENV_PREFIX=MIC_     #mandatory when using MIC
export MIC_OMP_PREFIX=MIC_     #mandatory when using MIC
export MIC_OMP_NUM_THREADS=224 #mandatory when using MIC
export MIC_KMP_PLACETHREADS="56c,4t" #optional
export MIC_KMP_AFFINITY="granularity=fine,compact" #optional
export OFFLOAD_REPORT=2        #optional
export QFORCE_PATH=/home/dima/Projects/QFORCE
export KMP_AFFINITY=verbose,compact
export MKL_NUM_THREADS=$OMP_NUM_THREADS

rm *.tmp *.log *.out
cp $QFORCE_PATH/qforce.v13.01.x ./
nvprof --log-file nv_profile.log --print-gpu-trace ./qforce.v13.01.x #&> qforce.log &
#nvprof --log-file nv_profile.log --print-gpu-trace --metrics branch_efficiency,gld_efficiency,gst_efficiency ./qforce.v13.01.x #&> qforce.log &
#/usr/bin/mpiexec.mpich2 -n $QF_NUM_PROCS ./qforce.v13.01.x #&> qforce.log &
#gprof ./qforce.v13.01.x gmon.out > profile.log
