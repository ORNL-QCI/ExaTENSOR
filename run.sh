#It is crucial to launch MPI processes consecutively within a node!
#Environment variable QF_PROCS_PER_NODE must be set appropriately!
export QFORCE_PATH=/autofs/na3_home1/div/src/ExaTensor #mandatory
export OMP_NUM_THREADS=8                            #mandatory
export QF_NUM_PROCS=1                               #mandatory
export QF_PROCS_PER_NODE=1                          #mandatory
export QF_GPUS_PER_PROCESS=1                        #optional (Nvidia GPU)
export QF_MICS_PER_PROCESS=0                        #optional (Intel Xeon Phi)
export QF_AMDS_PER_PROCESS=0                        #optional (AMD GPU)
export MIC_PREFIX=MIC                               #mandatory when using MIC
export MIC_ENV_PREFIX=MIC                           #mandatory when using MIC
export MIC_OMP_PREFIX=MIC                           #mandatory when using MIC
export MIC_OMP_NUM_THREADS=224                      #mandatory when using MIC
export MIC_MKL_NUM_THREADS=$MIC_OMP_NUM_THREADS     #mandatory when using MIC
export MIC_KMP_PLACE_THREADS="56c,4t"               #optional (MIC only)
export MIC_KMP_AFFINITY="granularity=fine,compact"  #optional (MIC only)
export MKL_MIC_ENABLE=1                             #optional (MIC only)
export OFFLOAD_REPORT=2                             #optional (MIC only)
export KMP_AFFINITY=compact                         #optional (CPU only)
#export MKL_NUM_THREADS=$OMP_NUM_THREADS            #optional (CPU only: MKL)

rm *.tmp *.log *.out *.x
cp $QFORCE_PATH/qforce.v13.01.x ./
aprun -n $QF_NUM_PROCS -N $QF_PROCS_PER_NODE -d $OMP_NUM_THREADS -m 16384 ./qforce.v13.01.x
#nvprof --log-file nv_profile.log --print-gpu-trace ./qforce.v13.01.x # &> qforce.log &
#nvprof --log-file nv_profile.log --print-gpu-trace --metrics branch_efficiency,gld_efficiency,gst_efficiency ./qforce.v13.01.x # &> qforce.log &
#gprof ./qforce.v13.01.x gmon.out > profile.log
