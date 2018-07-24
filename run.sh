#It is crucial to launch MPI processes consecutively within a node
#if multiple MPI processes reside on the same node. In this case
#the environment variable QF_PROCS_PER_NODE must be set appropriately!

export QF_PATH=/home/div/src/ExaTensor #full path to ExaTensor
export QF_NUM_PROCS=8                  #number of MPI processes
export QF_PROCS_PER_NODE=8             #number of MPI processes per node
export QF_CORES_PER_PROC=1             #number of cores per MPI process
export QF_NUM_THREADS=8                #number of threads per MPI process (at least 8)
export QF_GPUS_PER_PROCESS=0           #number of discrete NVIDIA GPU's per process (optional)
export QF_MICS_PER_PROCESS=0           #number of discrete Intel Xeon Phi's per process (optional)
export QF_AMDS_PER_PROCESS=0           #number of dsicrete AMD GPU's per process (optional)

#OpenMP:
export OMP_NUM_THREADS=$QF_NUM_THREADS #initial number of OpenMP threads per MPI process
export OMP_MAX_ACTIVE_LEVELS=3         #max number of OpenMP nesting levels (at least 3)
export OMP_THREAD_LIMIT=256            #max total number of OpenMP threads per process
export OMP_DYNAMIC=FALSE               #no OpenMP dynamic threading
export OMP_NESTED=TRUE                 #OpenMP nested parallelism is mandatory
export OMP_WAIT_POLICY=PASSIVE         #idle thread behavior

#Intel specific:
export KMP_AFFINITY="verbose,granularity=core,compact" #Intel CPU thread affinity
export MIC_PREFIX=MIC                                  #mandatory when using MIC
export MIC_ENV_PREFIX=MIC                              #mandatory when using MIC
export MIC_OMP_PREFIX=MIC                              #mandatory when using MIC
export MIC_OMP_NUM_THREADS=236                         #mandatory when using MIC
export MIC_MKL_NUM_THREADS=$MIC_OMP_NUM_THREADS        #mandatory when using MIC (Intel MIC MKL)
export MIC_KMP_PLACE_THREADS="59c,4t"                  #Intel MIC thread placement
export MIC_KMP_AFFINITY="granularity=fine,compact"     #Intel MIC thread affinity
export MIC_USE_2MB_BUFFERS=64K                         #Intel MIC only
export MKL_MIC_ENABLE=0                                #Intel MIC MKL auto-offloading
export MKL_NUM_THREADS=$OMP_NUM_THREADS                #number of Intel MKL threads per process
export OFFLOAD_REPORT=2                                #Intel MIC offload reporting level

#Cray specific:
export CRAY_OMP_CHECK_AFFINITY=TRUE          #CRAY: Show thread placement
export MPICH_NEMESIS_ASYNC_PROGRESS="SC"     #CRAY: Activate MPI asynchronous progress thread {"SC","MC"}
export MPICH_MAX_THREAD_SAFETY=multiple      #CRAY: Required for MPI asynchronous progress
export MPICH_GNI_ASYNC_PROGRESS_TIMEOUT=0    #CRAY:
export MPICH_GNI_MALLOC_FALLBACK=enabled     #CRAY:
export MPICH_RMA_OVER_DMAPP=1                #CRAY: DMAPP backend for CRAY-MPICH
#export _DMAPPI_NDREG_ENTRIES=16384          #CRAY: Max number of entries in UDREG memory registration cache
#export MPICH_ALLOC_MEM_HUGE_PAGES=1
#export MPICH_ALLOC_MEM_HUGEPG_SZ=2M
#export MPICH_ENV_DISPLAY=1
#export MPICH_GNI_MEM_DEBUG_FNAME=MPICH.memdebug
#export MPICH_RANK_REORDER_DISPLAY=1

#Summit (jsrun) specific:
unset PAMI_IBV_ENABLE_DCT

rm *.tmp *.log *.out *.x
cp $QF_PATH/Qforce.x ./

#/usr/local/mpi/openmpi-3.1.0/bin/mpiexec -n $QF_NUM_PROCS -npernode $QF_PROCS_PER_NODE -oversubscribe ./Qforce.x #>& qforce.log

#/usr/local/mpi/mpich-3.2/bin/mpiexec -n $QF_NUM_PROCS ./Qforce.x #>& qforce.log

#aprun -n $QF_NUM_PROCS -N $QF_PROCS_PER_NODE -d $QF_CORES_PER_PROC -cc 0,2,4,6,8,10,12,14,1,3,5,7,9,11,13 -r1 ./Qforce.x #>& qforce.log

#aprun -n $QF_NUM_PROCS -N $QF_PROCS_PER_NODE -cc none ./Qforce.x #>& qforce.log

#nvprof --log-file nv_profile.log --print-gpu-trace ./Qforce.x #>& qforce.log

#nvprof --log-file nv_profile.log --print-gpu-trace --metrics branch_efficiency,gld_efficiency,gst_efficiency ./Qforce.x #>& qforce.log

#gprof ./Qforce.x gmon.out > profile.log
