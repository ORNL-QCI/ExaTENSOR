#!/bin/bash

#BSUB -nnodes 129
#BSUB -W 00:20
#BSUB -P STF006
#BSUB -J exatensor
#BSUB -o exa.out.%J
#BSUB -e exa.err.%J
###BSUB -alloc_flags "gpumps smt1"

# Load modules:
module load gcc/8.1.1
module load essl
module load cuda
module load git
module load cmake

# Set full path to the executable:
export EXA_PATH=/ccs/home/div/src/ExaTensor_Summit_gcc8

#ExaTENSOR specific:
export QF_NUM_PROCS=258           #total number of MPI processes
export QF_PROCS_PER_NODE=2        #number of MPI processes per logical node (logical nodes are created by node resource isolation)
export QF_CORES_PER_PROCESS=21    #number of physical CPU cores per MPI process (no less than 1)
export QF_MEM_PER_PROCESS=130000  #host RAM memory limit per MPI process in MB
export QF_NVMEM_PER_PROCESS=0     #non-volatile memory limit per MPI process in MB
export QF_HOST_BUFFER_SIZE=180000 #host buffer size per MPI process in MB (must be less than QF_MEM_PER_PROCESS)
export QF_GPUS_PER_PROCESS=3      #number of discrete NVIDIA GPU's per MPI process (optional)
export QF_MICS_PER_PROCESS=0      #number of discrete Intel Xeon Phi's per MPI process (optional)
export QF_AMDS_PER_PROCESS=0      #number of discrete AMD GPU's per MPI process (optional)
export QF_NUM_THREADS=21          #initial number of CPU threads per MPI process (irrelevant, keep it 8)

#OpenMP generic:
export OMP_NUM_THREADS=$QF_NUM_THREADS #initial number of OpenMP threads per MPI process
export OMP_DYNAMIC=false               #no OpenMP dynamic threading
export OMP_NESTED=true                 #OpenMP nested parallelism is mandatory
export OMP_MAX_ACTIVE_LEVELS=3         #max number of OpenMP nesting levels (at least 3)
export OMP_THREAD_LIMIT=256            #max total number of OpenMP threads per process
export OMP_WAIT_POLICY=PASSIVE         #idle thread behavior
#export OMP_STACKSIZE=200M             #stack size per thread
#export OMP_DISPLAY_ENV=VERBOSE        #display OpenMP environment variables
#export GOMP_DEBUG=1                   #GNU OpenMP debugging
#export LOMP_DEBUG=1                   #IBM XL OpenMP debugging

#OpenMP thread binding:
export OMP_PLACES_DEFAULT=cores                                      #default thread binding to CPU logical cores
export OMP_PLACES_EOS="{1},{3},{5},{7,9},{0:16:2},{11},{13},{15}"    #Eos 16-core hyperthreaded Intel Xeon thread binding (even logical cores do computing)
export OMP_PLACES_TITAN="{1},{3},{5},{7,9},{0:8:2},{11},{13},{15}"   #Titan 16-core 8-FPU AMD thread binding (even logical cores do computing)
#export OMP_PLACES_POWER9="{0:4},{4:4},{8:4},{12:4},{28:56},{16:4},{20:4},{24:4}" #Summit 21-core SMT4 Power9 socket thread binding (even logical cores do computing)
export OMP_PLACES_POWER9="{0},{4},{8},{12},{28:56},{16},{20},{24}"   #Summit 21-core SMT4 Power9 socket thread binding (even logical cores do computing)
export OMP_PLACES_KNL="{1},{3},{5},{7,9},{0:128:2},{11},{13},{15}"   #Percival 64-core SMT4 KNL thread binding (even logical cores do computing)
#export OMP_PLACES=$OMP_PLACES_POWER9
export OMP_PROC_BIND="close,spread,spread" #nest1: Functional threads (DSVU)
                                           #nest2: TAVP-WRK:Dispatcher spawns coarse-grain Executors
                                           #nest3: TAVP-WRK:Dispatcher:Executor spawns execution threads in CP-TAL kernels
#MKL specific:
export MKL_NUM_THREADS_DEFAULT=1                #keep consistent with chosen OMP_PLACES!
export MKL_NUM_THREADS_EOS=16                   #keep consistent with chosen OMP_PLACES!
export MKL_NUM_THREADS_TITAN=8                  #keep consistent with chosen OMP_PLACES!
export MKL_NUM_THREADS_POWER9=56                #keep consistent with chosen OMP_PLACES!
export MKL_NUM_THREADS_KNL=128                  #keep consistent with chosen OMP_PLACES!
export MKL_NUM_THREADS=$MKL_NUM_THREADS_DEFAULT #number of Intel MKL threads per process
export MKL_DYNAMIC=false

#Intel MIC specific:
#export KMP_AFFINITY="verbose,granularity=core,compact"     #Intel CPU thread affinity
#export MIC_PREFIX=MIC                                      #mandatory when using MIC
#export MIC_ENV_PREFIX=MIC                                  #mandatory when using MIC
#export MIC_OMP_PREFIX=MIC                                  #mandatory when using MIC
#export MIC_OMP_NUM_THREADS=256                             #mandatory when using MIC
#export MIC_MKL_NUM_THREADS=$MIC_OMP_NUM_THREADS            #mandatory when using MIC (Intel MIC MKL)
#export MIC_KMP_PLACE_THREADS="64c,4t"                      #Intel MIC thread placement
#export MIC_KMP_AFFINITY="verbose,granularity=fine,compact" #Intel MIC thread affinity
#export MIC_USE_2MB_BUFFERS=64K                             #Intel MIC only
#export MKL_MIC_ENABLE=0                                    #Intel MIC MKL auto-offloading
#export OFFLOAD_REPORT=2

#Cray/MPICH specific:
#export CRAY_OMP_CHECK_AFFINITY=TRUE         #CRAY: Show thread placement
export MPICH_MAX_THREAD_SAFETY=multiple      #CRAY: Required for MPI asynchronous progress
export MPICH_NEMESIS_ASYNC_PROGRESS="MC"     #CRAY: Activate MPI asynchronous progress thread {"SC","MC"}
#export MPICH_RMA_OVER_DMAPP=1               #CRAY: DMAPP backend for CRAY-MPICH
#export MPICH_GNI_ASYNC_PROGRESS_TIMEOUT=0   #CRAY:
#export MPICH_GNI_MALLOC_FALLBACK=enabled    #CRAY:
#export MPICH_ALLOC_MEM_HUGE_PAGES=1         #CRAY: Huge pages
#export MPICH_ALLOC_MEM_HUGEPG_SZ=2M         #CRAY: Huge page size
#export _DMAPPI_NDREG_ENTRIES=16384          #CRAY: Max number of entries in UDREG memory registration cache
#export MPICH_ENV_DISPLAY=1
#export MPICH_GNI_MEM_DEBUG_FNAME=MPICH.memdebug
#export MPICH_RANK_REORDER_DISPLAY=1

#Summit specific:
export PAMI_IBV_ADAPTER_AFFINITY=1
export PAMI_IBV_DEVICE_NAME="mlx5_0:1,mlx5_3:1"
export PAMI_IBV_DEVICE_NAME_1="mlx5_3:1,mlx5_0:1"
export PAMI_IBV_ENABLE_OOO_AR=1                   #adaptive routing is default
export PAMI_ENABLE_STRIPING=1                     #increases network bandwidth, also increases latency
export PAMI_IBV_DISABLE_ODP=0                     #ODP (requires CAPI for performance)
#export PAMI_IBV_ENABLE_TAG_MATCHING=1            #hardware tag matching
export PAMI_IBV_ENABLE_DCT=1                      #reduces MPI_Init() time at large scale
#unset PAMI_IBV_ENABLE_DCT
#export PAMI_IBV_QP_SERVICE_LEVEL=8
#export PAMI_PMIX_DATACACHE=1
#export PAMI_IBV_DEBUG_CQE=1                      #CQE error debugging
#export PAMI_IBV_DEBUG_QP_TIMEOUT=22
#export PAMI_IBV_DEBUG_RNR_RETRY=9
#export OMPI_LD_PRELOAD_POSTPEND=$OLCF_SPECTRUM_MPI_ROOT/lib/libmpitrace.so

ulimit -s unlimited

rm core.* *.tmp *.log *.out *.x
cp $EXA_PATH/Qforce.x ./

#One MPI process per node (on only one socket):
#jsrun --smpiargs='-mca common_pami_use_odp 1' -D PAMI_IBV_DISABLE_ODP=0 -n $QF_NUM_PROCS -r $QF_PROCS_PER_NODE -a 1 -c $QF_CORES_PER_PROCESS -g $QF_GPUS_PER_PROCESS -bnone ./Qforce.x

#Two MPI processes per node (one MPI process per socket):
#jsrun --smpiargs='-mca common_pami_use_odp 1' -D PAMI_IBV_DISABLE_ODP=0 -n $QF_NUM_PROCS -r $QF_PROCS_PER_NODE -a 1 -c $QF_CORES_PER_PROCESS -g $QF_GPUS_PER_PROCESS -bnone ./exec_exa.sh
#jsrun --smpiargs='-mca pml_pami_use_2sided 1' -E PAMI_IBV_ENABLE_TAG_MATCHING=1 --smpiargs='-mca common_pami_use_odp 1' -D PAMI_IBV_DISABLE_ODP=0 -n $QF_NUM_PROCS -r $QF_PROCS_PER_NODE -a 1 -c $QF_CORES_PER_PROCESS -g $QF_GPUS_PER_PROCESS -bnone ./exec_exa.sh
jsrun --smpiargs='-async' --smpiargs='-mca common_pami_use_odp 1' -D PAMI_IBV_DISABLE_ODP=0 -n $QF_NUM_PROCS -r $QF_PROCS_PER_NODE -a 1 -c $QF_CORES_PER_PROCESS -g $QF_GPUS_PER_PROCESS -bnone ./exec_exa.sh
