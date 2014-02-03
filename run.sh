#It is crucial to launch MPI processes consecutively within a node!
#Environment variable QF_PROCS_PER_NODE must be set appropriately!
export OMP_NUM_THREADS=2     #mandatory: max number of threads per MPI process
export QF_PROCS_NUM=3        #mandatory: total number of MPI processes
export QF_PROCS_PER_NODE=1   #mandatory: number of MPI processes per node
export QF_GPUS_PER_PROCESS=1 #optional: max number of GPUs per MPI process
export QF_MICS_PER_PROCESS=0 #optional: max number of MICs per MPI process
export QFORCE_PATH=/home/dima/Projects/QFORCE #QFORCE path
export KMP_AFFINITY=verbose,compact #set thread affinity for Intel processors

export MKL_NUM_THREADS=$OMP_NUM_THREADS
rm *.log *.tmp
mpiexec -n $QF_PROCS_NUM $QFORCE_PATH/qforce.v13.01.x
