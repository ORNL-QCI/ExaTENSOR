#It is crucial to launch MPI processes consecutively within a node
#if multiple MPI processes reside on the same node. In this case
#the environment variable QF_PROCS_PER_NODE must be set appropriately!
export QF_PATH=/ccs/home/div/src/ExaTensor          #full path to ExaTensor
#export QF_PATH=/home/dima/Projects/QFORCE
export QF_NUM_PROCS=1                               #number of MPI processes
export QF_CORES_PER_PROC=16                         #number of cores per MPI process
export QF_PROCS_PER_NODE=1                          #number of MPI processes per node
export OMP_NUM_THREADS=8                            #max number of OpenMP threads per process
export QF_GPUS_PER_PROCESS=1                        #number of Nvidia GPU's per process (optional)
export QF_MICS_PER_PROCESS=0                        #number of Intel Xeon Phi's per process (optional)
export QF_AMDS_PER_PROCESS=0                        #number of AMD GPU's per process (optional)
export MKL_NUM_THREADS=$OMP_NUM_THREADS             #number of Intel MKL threads per process (optional)
export OMP_MAX_ACTIVE_LEVELS=3                      #max number of OpenMP nesting levels (at least 3)
export OMP_THREAD_LIMIT=64                          #max total number of OpenMP threads per process
export OMP_DYNAMIC=FALSE                            #no OpenMP dynamic threading
export OMP_NESTED=TRUE                              #OpenMP nested parallelism is mandatory
export OMP_WAIT_POLICY=PASSIVE                      #idle thread behavior (optional)
#export KMP_AFFINITY=compact                        #Intel CPU thread affinity (optional)
export MIC_PREFIX=MIC                               #mandatory when using MIC
export MIC_ENV_PREFIX=MIC                           #mandatory when using MIC
export MIC_OMP_PREFIX=MIC                           #mandatory when using MIC
export MIC_OMP_NUM_THREADS=236                      #mandatory when using MIC
export MIC_MKL_NUM_THREADS=$MIC_OMP_NUM_THREADS     #mandatory when using MIC (Intel MIC MKL)
export MIC_KMP_PLACE_THREADS="59c,4t"               #optional (MIC only)
export MIC_KMP_AFFINITY="granularity=fine,compact"  #optional (MIC only)
export MIC_USE_2MB_BUFFERS=64K                      #optional (MIC only)
export MKL_MIC_ENABLE=0                             #optional (MIC only: MKL MIC auto-offloading)
export OFFLOAD_REPORT=2                             #optional (MIC only)

#CRAY specific:
export CRAY_OMP_CHECK_AFFINITY=TRUE                 #CRAY: Shows thread placement
export MPICH_NEMESIS_ASYNC_PROGRESS="SC"            #CRAY: Activates MPI asynchronous progress thread
export MPICH_MAX_THREAD_SAFETY=multiple             #CRAY: Required for MPI asynchronous progress
export MPICH_GNI_ASYNC_PROGRESS_TIMEOUT=0           #CRAY:
export MPICH_GNI_MALLOC_FALLBACK=enabled            #CRAY:
export MPICH_RMA_OVER_DMAPP=1                       #CRAY: DMAPP backend for CRAY-MPICH
export _DMAPPI_NDREG_ENTRIES=16384                  #CRAY: Max number of entries in UDREG memory registration cache
#export MPICH_ALLOC_MEM_HUGE_PAGES=1
#export MPICH_ALLOC_MEM_HUGEPG_SZ=2M
#export MPICH_ENV_DISPLAY=1
#export MPICH_GNI_MEM_DEBUG_FNAME=MPICH.memdebug
#export MPICH_RANK_REORDER_DISPLAY=1

rm *.tmp *.log *.out *.x
cp $QF_PATH/qforce.v13.01.x ./
aprun -n $QF_NUM_PROCS -N $QF_PROCS_PER_NODE -d $QF_CORES_PER_PROC -cc 0,2,4,6,8,10,12,14,1,3,5,7,9,11,13 -r1 ./Qforce.x #> qforce.log
#mpiexec -n $QF_NUM_PROCS -npernode $QF_PROCS_PER_NODE ./Qforce.x #> qforce.log
#nvprof --log-file nv_profile.log --print-gpu-trace ./Qforce.x # &> qforce.log &
#nvprof --log-file nv_profile.log --print-gpu-trace --metrics branch_efficiency,gld_efficiency,gst_efficiency ./Qforce.x # &> qforce.log &
#gprof ./Qforce.x gmon.out > profile.log
