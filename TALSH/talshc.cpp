/** ExaTensor::TAL-SH: Device-unified user-level API.
REVISION: 2016/05/02

Copyright (C) 2014-2016 Dmitry I. Lyakh (Liakh)
Copyright (C) 2014-2016 Oak Ridge National Laboratory (UT-Battelle)

This file is part of ExaTensor.

ExaTensor is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ExaTensor is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with ExaTensor. If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------------------
**/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "talsh.h"

//GLOBALS:
// General:
static int talsh_on=0;           //TAL-SH initialization flag (1:initalized; 0:not)
static clock_t talsh_begin_time; //TAL-SH begin time (zero time reference)
// Accelerator configuration (`Needs modification for non-contiguous subranges):
static int talsh_gpu_beg;        //first Nvidia GPU in the range `Obsolete
static int talsh_gpu_end;        //last Nvidia GPU in the range `Obsolete
// Device status:
static int talsh_cpu=DEV_OFF;
static int talsh_gpu[MAX_GPUS_PER_NODE]={DEV_OFF}; //current GPU status: {DEV_OFF,DEV_ON,DEV_ON_BLAS}
static int talsh_mic[MAX_MICS_PER_NODE]={DEV_OFF}; //current MIC status: {DEV_OFF,DEV_ON,DEV_ON_BLAS}
static int talsh_amd[MAX_AMDS_PER_NODE]={DEV_OFF}; //current AMD status: {DEV_OFF,DEV_ON,DEV_ON_BLAS}
// Failure statistics:
static unsigned long long int not_clean_count=0LL; //number of times a NOT_CLEAN status was returned (possible indication of a memory leak)

//INTERNAL TYPES:
// Host task:
typedef struct{
 int task_error; //task error code (-1:empty or in progress; 0:success; >0:error code)
} host_task_t;

//PROTOTYPES OF IMPORTED FUNCTIONS:
#ifdef __cplusplus
extern "C"{
#endif
// Fortran tensor block aliasing:
int talsh_tensor_f_assoc(talsh_tens_t * talsh_tens, int image_id);
int talsh_tensor_f_dissoc(talsh_tens_t * talsh_tens);
#ifdef __cplusplus
}
#endif

//PROTOTYPES OF INTERNAL FUNCTIONS:
#ifdef __cplusplus
extern "C"{
#endif
// Tensor block image info (exported to talshf.F90):
int talsh_tensor_image_info(const talsh_tens_t * talsh_tens, int image_id, int * dev_id, int * data_kind, void ** gmem_p, int * buf_entry);
#ifdef __cplusplus
}
#endif
// Error counters:
static void talsh_raise_not_clean();
// Host task API:
static int host_task_create(host_task_t ** host_task);
static int host_task_clean(host_task_t * host_task);
static int host_task_destroy(host_task_t * host_task);
// C tensor block aliasing:
static int talsh_tensor_c_assoc(talsh_tens_t * talsh_tens, int image_id);
static int talsh_tensor_c_dissoc(talsh_tens_t * talsh_tens);
// Construct a TAL-SH task:
static int talshTaskConstruct(talsh_task_t * talsh_task, int dev_kind, int data_kind = NO_TYPE);

//INTERNAL FUNCTIONS:
// Error counters:
static void talsh_raise_not_clean(){++not_clean_count;}
// Host task API:
static int host_task_create(host_task_t ** host_task)
{
 *host_task=(host_task_t*)malloc(sizeof(host_task_t));
 if(*host_task == NULL) return TRY_LATER;
 return host_task_clean(*host_task);
}

static int host_task_clean(host_task_t * host_task)
{
 if(host_task == NULL) return TALSH_INVALID_ARGS;
 host_task->task_error=-1;
 return TALSH_SUCCESS;
}

static int host_task_destroy(host_task_t * host_task)
{
 if(host_task == NULL) return TALSH_INVALID_ARGS;
 free(host_task);
 return TALSH_SUCCESS;
}

int talsh_tensor_image_info(const talsh_tens_t * talsh_tens, int image_id,
                            int * dev_id, int * data_kind, void ** gmem_p, int * buf_entry)
{
 talsh_dev_rsc_t *drsc;

 if(talsh_tens == NULL) return TALSH_INVALID_ARGS;
 if(talshTensorIsEmpty(talsh_tens) == YEP) return TALSH_OBJECT_IS_EMPTY;
 if(talsh_tens->ndev <= 0 || talsh_tens->ndev > talsh_tens->dev_rsc_len ||
    talsh_tens->dev_rsc == NULL || talsh_tens->data_kind == NULL) return TALSH_FAILURE;
 if(image_id < 0 || image_id >= talsh_tens->ndev) return TALSH_INVALID_ARGS;
 drsc=&(talsh_tens->dev_rsc[image_id]);
 if(tensDevRsc_is_empty(drsc) == YEP) return TALSH_FAILURE;
 *data_kind=talsh_tens->data_kind[image_id];
 *dev_id=drsc->dev_id;
 *gmem_p=drsc->gmem_p;
 *buf_entry=drsc->buf_entry;
 return TALSH_SUCCESS;
}

static int talsh_tensor_c_assoc(talsh_tens_t * talsh_tens, //inout: TAL-SH tensor
                                int image_id)              //in: id of the tensor body image to be used
/** Fills in the <tensBlck_t> component inside the <talsh_tens_t> object by importing the information
    from <talsh_tens>. A return status TRY_LATER indicates temporary shortage in available resources. **/
{
 int i,errc;
 tensBlck_t *ctens;
 talsh_dev_rsc_t *src_rsc_p;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 if(talsh_tens == NULL) return TALSH_INVALID_ARGS;
 if(talshTensorIsEmpty(talsh_tens) == YEP) return TALSH_INVALID_ARGS;
 if(image_id < 0 || image_id >= talsh_tens->ndev) return TALSH_INVALID_ARGS;
 if(talsh_tens->tensC != NULL) return TALSH_OBJECT_NOT_EMPTY;
 if(talsh_tens->dev_rsc == NULL || talsh_tens->data_kind == NULL) return TALSH_FAILURE;
 if(tens_valid_data_kind(talsh_tens->data_kind[image_id]) != YEP) return TALSH_FAILURE;
 src_rsc_p=&(talsh_tens->dev_rsc[image_id]);
 errc=tensBlck_create(&ctens); if(errc){if(errc != TRY_LATER) errc=TALSH_FAILURE; return errc;}
 errc=tensBlck_construct(ctens,YEP,talsh_tens->shape_p->num_dim,talsh_tens->shape_p->dims,talsh_tens->shape_p->divs,talsh_tens->shape_p->grps);
 if(errc){if(errc != TRY_LATER) errc=TALSH_FAILURE; i=tensBlck_destroy(ctens); return errc;}
 errc=tensBlck_attach_body(ctens,talsh_tens->data_kind[image_id],src_rsc_p->dev_id,src_rsc_p->gmem_p,src_rsc_p->buf_entry);
 if(errc){if(errc != TRY_LATER) errc=TALSH_FAILURE; i=tensBlck_destroy(ctens); return errc;}
 talsh_tens->tensC=(void*)ctens; //.tensC has the right shape, data_kind, and source data
 return TALSH_SUCCESS;
}

static int talsh_tensor_c_dissoc(talsh_tens_t * talsh_tens) //inout: TAL-SH tensor
/** Destroys the <tensBlck_t> component in the <talsh_tens_t> object. **/
{
 int errc;
 tensBlck_t *ctens;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 if(talsh_tens == NULL) return TALSH_INVALID_ARGS;
 if(talshTensorIsEmpty(talsh_tens) == YEP) return TALSH_INVALID_ARGS;
 if(talsh_tens->tensC == NULL) return TALSH_OBJECT_IS_EMPTY;
 ctens=(tensBlck_t*)talsh_tens->tensC; errc=tensBlck_destroy(ctens);
 if(errc){if(errc != NOT_CLEAN) errc=TALSH_FAILURE;}
 return errc;
}

//EXPORTED FUNCTIONS:
// TAL-SH control API:
int talshInit(size_t * host_buf_size,    //inout: Host Argument Buffer size in bytes (in: suggested; out: actual)
              int * host_arg_max,        //out: Max number of arguments that can fit into the Host Argument Buffer
              int ngpus, int gpu_list[], //in: number of Nvidia GPU(s) to use and the list of Nvidia GPU(s) to use
              int nmics, int mic_list[], //in: number of Intel Xeon Phi(s) to use and the list of Intel Xeon Phi(s) to use
              int namds, int amd_list[]) //in: number of AMD GPU(s) to use and the list of AMD GPU(s) to use
/** Initializes the TAL-SH runtime. **/
{
 int i,j,gpu_beg,gpu_end,errc;

 if(talsh_on) return TALSH_ALREADY_INITIALIZED;
//CPU Host:
#ifndef NO_BLAS
 talsh_cpu=DEV_ON_BLAS;
#else
 talsh_cpu=DEV_ON;
#endif
//NVidia GPU accelerators:
#ifndef NO_GPU
 if(ngpus > 0){
  if(ngpus > MAX_GPUS_PER_NODE) return TALSH_INVALID_ARGS;
  gpu_beg=gpu_list[0]; gpu_end=gpu_list[ngpus-1]; //`Allow for non-consecutive GPU ranges in arg_buf_allocate()
  if(gpu_beg < 0 || gpu_beg >= MAX_GPUS_PER_NODE) return TALSH_INVALID_ARGS;
  if(gpu_end < 0 || gpu_end >= MAX_GPUS_PER_NODE) return TALSH_INVALID_ARGS;
  for(i=1;i<ngpus;i++){
   if(gpu_list[i] != gpu_list[i-1]+1){
    printf("#FATAL(TALSH::talshInit): The current version only supports consecutive GPU ranges!");
    return TALSH_FAILURE;
   }
  }
 }else{
#endif
  gpu_beg=0; gpu_end=-1;
#ifndef NO_GPU
 }
#endif
//Intel Xeon Phi accelerators:
#ifndef NO_MIC
 if(nmics > 0){
  printf("#FATAL(TALSH::talshInit): Intel Xeon Phi is not fully supported yet!");
  return TALSH_NOT_IMPLEMENTED; //`Future
 }
#endif
//AMD GPU accelerators:
#ifndef NO_AMD
 if(namds > 0){
  printf("#FATAL(TALSH::talshInit): AMD GPU is not supported yet!");
  return TALSH_NOT_IMPLEMENTED; //`Future
 }
#endif
 errc=arg_buf_allocate(host_buf_size,host_arg_max,gpu_beg,gpu_end); if(errc) return TALSH_FAILURE;
#ifndef NO_GPU
 for(i=0;i<ngpus;i++){
  j=gpu_list[i]; if(j < 0 || j >= MAX_GPUS_PER_NODE) return TALSH_INVALID_ARGS;
  talsh_gpu[j]=gpu_is_mine(j);
 }
#endif
 talsh_gpu_beg=gpu_beg; talsh_gpu_end=gpu_end;
 talsh_on=1; talsh_begin_time=clock();
 return TALSH_SUCCESS;
}

int talshShutdown()
/** Shuts down the TAL-SH runtime. **/
{
 int i,errc;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 errc=arg_buf_deallocate(talsh_gpu_beg,talsh_gpu_end);
 talsh_gpu_beg=0; talsh_gpu_end=-1; talsh_on=0;
 talsh_cpu=DEV_OFF;
 for(i=0;i<MAX_GPUS_PER_NODE;i++) talsh_gpu[i]=DEV_OFF;
 for(i=0;i<MAX_MICS_PER_NODE;i++) talsh_mic[i]=DEV_OFF;
 for(i=0;i<MAX_AMDS_PER_NODE;i++) talsh_amd[i]=DEV_OFF;
 if(errc) return TALSH_FAILURE;
 return TALSH_SUCCESS;
}

int talshFlatDevId(int dev_kind, //in: device kind
                   int dev_num)  //in: device Id within its kind (0..MAX)
/** Converts a kind-specific device Id into the flat device Id.
    DEV_MAX return status indicates invalidity of the arguments. **/
{
 return encode_device_id(dev_kind,dev_num); //Success:[0..DEV_MAX-1]; Failure: DEV_MAX
}

int talshKindDevId(int dev_id,     //in: flat device Id: [0:DEV_MAX-1]
                   int * dev_kind) //out: device kind
/** Converts a flat device Id into the kind specific device Id.
    A negative return value indicates an invalid flat device Id. **/
{
 return decode_device_id(dev_id,dev_kind);
}

int talshDeviceState(int dev_num,  //in: either a flat or kind specific (when <dev_kind> is present) device id
                     int dev_kind) //in: device kind (note that it changes the meaning of the <dev_num> argument)
/** Returns device state (Success:[DEV_OFF,DEV_ON,DEV_ON_BLAS]) **/
{
 int devk,i,sts;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 if(dev_kind == DEV_NULL){
  i=talshKindDevId(dev_num,&devk);
  if(i < 0) return TALSH_INVALID_ARGS;
 }else{
  devk=dev_kind;
 }
 switch(devk){
  case DEV_HOST:
   sts=talsh_cpu;
   break;
  case DEV_NVIDIA_GPU:
   sts=talsh_gpu[i];
   break;
  case DEV_INTEL_MIC:
   sts=talsh_mic[i];
   break;
  case DEV_AMD_GPU:
   sts=talsh_amd[i];
   break;
  default:
   return TALSH_INVALID_ARGS;
 }
 return sts;
}

int talshDeviceState_(int dev_num, int dev_kind) //Fortran wrapper
{
 return talshDeviceState(dev_num,dev_kind);
}

int talshDeviceBusyLeast(int dev_kind) //in: device kind (defaults to any kind)
/** Returns the least busy device id. **/
{
 int i;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 switch(dev_kind){
  case DEV_NULL:
   return talshFlatDevId(DEV_HOST,0); //`if device kind not specified, return CPU Host for simplicity
  case DEV_HOST:
   return talshFlatDevId(DEV_HOST,0);
  case DEV_NVIDIA_GPU:
#ifndef NO_GPU
   i=gpu_busy_least();
   if(i < 0 || i >= MAX_GPUS_PER_NODE) return TALSH_FAILURE;
   return i;
#else
   return TALSH_NOT_AVAILABLE;
#endif
  case DEV_INTEL_MIC:
#ifndef NO_PHI
   return TALSH_NOT_IMPLEMENTED; //`Implement in future
#else
   return TALSH_NOT_AVAILABLE;
#endif
  case DEV_AMD_GPU:
#ifndef NO_AMD
   return TALSH_NOT_IMPLEMENTED; //`Implement in future
#else
   return TALSH_NOT_AVAILABLE;
#endif
 }
 return TALSH_INVALID_ARGS;
}

int talshDeviceBusyLeast_(int dev_kind) //Fortran wrapper
{
 return talshDeviceBusyLeast(dev_kind);
}

int talshStats(int dev_id,   //in: device id (either flat or kind specific device id, see below)
               int dev_kind) //in: device kind (if present, <dev_id> will be interpreted as kind specific)
/** Prints the run-time statistics for devices of interest. **/
{
 int rc=TALSH_SUCCESS,devk,devn;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 switch(dev_kind){
  case DEV_NULL:
   if(dev_id < 0){ //print stats for all active devices
    rc=talshStats(-1,DEV_HOST);
    rc=talshStats(-1,DEV_NVIDIA_GPU);
    rc=talshStats(-1,DEV_INTEL_MIC);
    rc=talshStats(-1,DEV_AMD_GPU);
    rc=TALSH_SUCCESS;
   }else{
    devn=talshKindDevId(dev_id,&devk);
    rc=talshStats(devn,devk);
   }
   break;
  case DEV_HOST:
   rc=TALSH_NOT_IMPLEMENTED; //`Implement
   break;
  case DEV_NVIDIA_GPU:
#ifndef NO_GPU
   rc=gpu_print_stats(dev_id);
#else
   rc=TALSH_NOT_AVAILABLE;
#endif
   break;
  case DEV_INTEL_MIC:
#ifndef NO_PHI
   rc=TALSH_NOT_IMPLEMENTED; //`Implement in future
#else
   rc=TALSH_NOT_AVAILABLE;
#endif
   break;
  case DEV_AMD_GPU:
#ifndef NO_AMD
   rc=TALSH_NOT_IMPLEMENTED; //`Implement in future
#else
   rc=TALSH_NOT_AVAILABLE;
#endif
   break;
  default:
   rc=TALSH_INVALID_ARGS;
 }
 return rc;
}

int talshStats_(int dev_id, int dev_kind) //Fortran wrapper
{
 return talshStats(dev_id,dev_kind);
}

// TAL-SH tensor block API:
int talshTensorCreate(talsh_tens_t ** tens_block) //out: pointer to a newly created empty tensor block
/** Returns a pointer to a newly created empty tensor block (0:Success; TRY_LATER:Short on memory). **/
{
 *tens_block=(talsh_tens_t*)malloc(sizeof(talsh_tens_t));
 if(*tens_block == NULL) return TRY_LATER;
 return talshTensorClean(*tens_block);
}

int talshTensorClean(talsh_tens_t * tens_block)
/** Cleans an undefined tensor block (default ctor) making it defined-empty. **/
{
 if(tens_block == NULL) return TALSH_INVALID_ARGS;
 tens_block->shape_p=NULL;
 tens_block->dev_rsc=NULL;
 tens_block->data_kind=NULL;
 tens_block->dev_rsc_len=0;
 tens_block->ndev=0;
 tens_block->tensF=NULL;
 tens_block->tensC=NULL;
 return TALSH_SUCCESS;
}

int talshTensorIsEmpty(const talsh_tens_t * tens_block)
/** Returns YEP if the tensor block is empty, NOPE otherwise. **/
{
 if(tens_block == NULL) return TALSH_INVALID_ARGS;
 if(tens_block->shape_p == NULL) return YEP;
 return NOPE;
}

int talshTensorConstruct(talsh_tens_t * tens_block,     //inout: empty tensor block on entrance, constructed tensor block on exit
                         int data_kind,                 //in: data kind: {R4,R8,C4,C8,NO_TYPE}
                         int tens_rank,                 //in: tensor block rank (number of dimensions)
                         const int tens_dims[],         //in: tensor block dimension extents
                         int dev_id,                    //in: flat device ID on which the tensor block will reside
                         void * ext_mem,                //in: pointer to externally provided memory for tensor elements
                         int in_hab,                    //in: if >=0, a non-NULL <ext_mem> points to the HAB entry #<in_hab>
                         talsh_tens_init_i init_method, //in: user-defined initialization method (function pointer)
                         double init_val_real,          //in: initialization value (real part), defaults to 0.0
                         double init_val_imag)          //in: initialization value (imaginary part), defaults to 0.0
/** Constructs a tensor block: {0: success; TRY_LATER: no enough free memory available; DEVICE_UNABLE: device unable}.
    If <data_kind> == NO_TYPE, the tensor body will not be allocated (only the tensor shape),
    unless an external storage is provided (<ext_mem>). In case the tensor body storage space
    is provided externally (<ext_mem> != NULL), the initialization step is skipped. In other cases,
    unless <data_kind>=NO_TYPE, the newly allocated tensor body will be initialized by a user-defined
    method, or, if no method is provided (NULL), by a user-defined value, which defaults to zero.
    If the tensor body initialization failed, a status NOT_CLEAN is returned but
    the tensor block is ready for use (except its body value is still undefined). **/
{
 int i,j,dev_num,dev_kind,dksize,errc,already_allocated,use_hab;
 size_t tvol,tsize;
 float fval;
 float *fp;
 double *dp;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 errc=TALSH_SUCCESS;
 //Check arguments:
 if(tens_block == NULL) return TALSH_INVALID_ARGS; //tensor block must have been preallocated
 if(talshTensorIsEmpty(tens_block) != YEP) return TALSH_OBJECT_NOT_EMPTY; //tensor block is not empty (destruct it first)
 if(tens_valid_data_kind(data_kind,&dksize) != YEP) return TALSH_INVALID_ARGS; //unknown data kind (NO_TYPE is a valid type here)
 dev_num=talshKindDevId(dev_id,&dev_kind); if(dev_num < 0) return TALSH_INVALID_ARGS; //invalid device id
 already_allocated=0; if(ext_mem != NULL) already_allocated=1; //check whether an external memory space is provided for the tensor body
 if(in_hab >= 0){use_hab=YEP;}else{in_hab=-1; use_hab=NOPE;}
 //Tensor shape:
 errc=tensShape_create(&(tens_block->shape_p)); if(errc == TRY_LATER || errc == DEVICE_UNABLE) return errc;
 if(errc != 0 || tens_block->shape_p == NULL) return TALSH_FAILURE;
 errc=tensShape_construct(tens_block->shape_p,NOPE,tens_rank,tens_dims); //NOPE = not pinned
 if(errc != 0 && errc != TRY_LATER && errc != DEVICE_UNABLE) errc=TALSH_FAILURE;
 if(errc != 0){i=talshTensorDestruct(tens_block); return errc;}
 //Device resource storage:
 if(tens_block->dev_rsc_len == 0 && tens_block->dev_rsc == NULL && tens_block->data_kind == NULL){ //tensor block must be defined-empty
  tens_block->dev_rsc=(talsh_dev_rsc_t*)malloc(TALSH_MAX_DEV_PRESENT*sizeof(talsh_dev_rsc_t));
  if(tens_block->dev_rsc != NULL){
   tens_block->dev_rsc_len=TALSH_MAX_DEV_PRESENT; tens_block->ndev=0;
   for(j=0;j<TALSH_MAX_DEV_PRESENT;++j){i=tensDevRsc_clean(&(tens_block->dev_rsc[j]));}
   tens_block->data_kind=(int*)malloc(TALSH_MAX_DEV_PRESENT*sizeof(int));
   if(tens_block->data_kind == NULL){i=talshTensorDestruct(tens_block); return TRY_LATER;}
   for(j=0;j<TALSH_MAX_DEV_PRESENT;++j){tens_block->data_kind[j]=NO_TYPE;}
  }else{
   i=talshTensorDestruct(tens_block);
   return TRY_LATER;
  }
 }else{
  i=talshTensorDestruct(tens_block);
  return TALSH_INVALID_ARGS;
 }
 //Tensor body:
 if(already_allocated){ //tensor body storage has been allocated outside (no initialization will be performed)
  errc=tensDevRsc_attach_mem(&(tens_block->dev_rsc[0]),dev_id,ext_mem,in_hab);
  if(errc){i=talshTensorDestruct(tens_block); return TALSH_FAILURE;}
  tens_block->data_kind[0]=data_kind; tens_block->ndev=1; //present on one device, even if NO_TYPE
 }else{ //tensor body storage needs to be allocated here (will also be initialized), unless NO_TYPE
  if(data_kind != NO_TYPE){
   tvol=talshTensorVolume(tens_block);
   if(tvol > 0){
    tsize=tvol*dksize;
    if(tsize <= 0){i=talshTensorDestruct(tens_block); return TALSH_INTEGER_OVERFLOW;}
    errc=tensDevRsc_allocate_mem(&(tens_block->dev_rsc[0]),dev_id,tsize,use_hab);
    if(errc != 0 && errc != TRY_LATER && errc != DEVICE_UNABLE) errc=TALSH_FAILURE;
    if(errc != 0){i=talshTensorDestruct(tens_block); return errc;}
    tens_block->data_kind[0]=data_kind; tens_block->ndev=1; //present on one device
   }else{
    i=talshTensorDestruct(tens_block);
    return TALSH_FAILURE;
   }
   //Initialization:
   if(tens_block->ndev > 0){
    if(dev_kind == DEV_HOST){ //`Currently supported only on Host
     if(init_method != NULL){
      init_method(tens_block->dev_rsc[0].gmem_p,data_kind,tens_rank,tens_dims,&errc);
      if(errc) errc=NOT_CLEAN; //initialization failed, tensor block value is undefined, but one may continue
     }else{
      switch(data_kind){
       case R4:
        fval = (float)init_val_real;
        fp = (float*)(tens_block->dev_rsc[0].gmem_p);
#pragma omp parallel for schedule(guided)
        for(size_t l=0; l < tvol; l++) fp[l]=fval;
        break;
       case R8:
        dp = (double*)(tens_block->dev_rsc[0].gmem_p);
#pragma omp parallel for schedule(guided)
        for(size_t l=0; l < tvol; l++) dp[l]=init_val_real;
        break;
       default:
        return NOT_CLEAN; //`Enable initialization for complex data kinds C4 and C8
      }
     }
    }else{
     errc=TALSH_NOT_IMPLEMENTED; //`Initialization on other device kinds should be enabled
    }
   }
  }
 }
 return errc;
}

int talshTensorConstruct_(talsh_tens_t * tens_block, int data_kind, int tens_rank, const int tens_dims[], int dev_id, //Fortran wrapper
                          void * ext_mem, int in_hab, talsh_tens_init_i init_method, double init_val_real, double init_val_imag)
{
 return talshTensorConstruct(tens_block, data_kind, tens_rank, tens_dims, dev_id,
                             ext_mem, in_hab, init_method, init_val_real, init_val_imag);
}

int talshTensorDestruct(talsh_tens_t * tens_block) //in: non-NULL pointer to a tensor block (empty tensor block on exit)
/** Destructs a tensor block and sets its status to empty. **/
{
 int i,j,errc;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 errc=TALSH_SUCCESS;
 if(tens_block == NULL) return TALSH_INVALID_ARGS;
 if(tens_block->tensF != NULL){free(tens_block->tensF); tens_block->tensF=NULL;}
 if(tens_block->tensC != NULL){free(tens_block->tensC); tens_block->tensC=NULL;}
 if(tens_block->shape_p != NULL){
  i=tensShape_destroy(tens_block->shape_p); tens_block->shape_p=NULL;
  if(i == 0 || i == NOT_CLEAN){if(errc == 0) errc=i;}else{errc=TALSH_FAILURE;}
 }
 if(tens_block->ndev > tens_block->dev_rsc_len){tens_block->ndev=tens_block->dev_rsc_len; errc=TALSH_FAILURE;}
 if(tens_block->dev_rsc != NULL){
  for(j=0;j<tens_block->ndev;j++){
   i=tensDevRsc_release_all(&(tens_block->dev_rsc[j]));
   if(i == 0 || i == NOT_CLEAN){if(errc == 0) errc=i;}else{errc=TALSH_FAILURE;}
  }
  free(tens_block->dev_rsc); tens_block->dev_rsc=NULL;
 }
 if(tens_block->data_kind != NULL){free(tens_block->data_kind); tens_block->data_kind=NULL;}
 i=talshTensorClean(tens_block); //set to an empty status
 return errc;
}

int talshTensorDestroy(talsh_tens_t * tens_block) //in: non-NULL pointer to a tensor block
/** Completely destroys a talsh_tens_t object. **/
{
 int errc;
 if(tens_block == NULL) return TALSH_INVALID_ARGS;
 errc=talshTensorDestruct(tens_block);
 free(tens_block);
 return errc;
}

size_t talshTensorVolume(const talsh_tens_t * tens_block) //in: tensor block
/** Returns the total number of elements in the tensor block.
    0 on return means the tensor block is empty. **/
{
 if(tens_block == NULL) return TALSH_INVALID_ARGS;
 if(talshTensorIsEmpty(tens_block) == YEP) return 0;
 return tensShape_volume(tens_block->shape_p);
}

int talshTensorShape(const talsh_tens_t * tens_block, talsh_tens_shape_t * tens_shape)
/** Returns the shape of the tensor block. The tensor shape object <tens_shape>
    passed here must either be either empty defined or value defined. It is errorneous
    to pass an undefined <tens_shape> object. **/
{
 int errc;

 if(tens_block == NULL || tens_shape == NULL) return TALSH_INVALID_ARGS;
 if(talshTensorIsEmpty(tens_block) == YEP) return TALSH_OBJECT_IS_EMPTY;
 errc=tensShape_construct(tens_shape,NOPE,tens_block->shape_p->num_dim,tens_block->shape_p->dims, //NOPE: not pinned
                          tens_block->shape_p->divs,tens_block->shape_p->grps);
 if(errc) errc=TALSH_FAILURE;
 return errc;
}

int talshTensorPresence(const talsh_tens_t * tens_block, int * ncopies, int copies[], int data_kinds[], int dev_kind, int dev_id)
/** Returns the list of devices on which a copy of the tensor block resides, together with the data kind.
    The presence of optional <dev_kind> and <dev_id> arguments further customizes the search. **/
{
 int i,j,m,devnum,devk,specific_kind,specific_device;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 *ncopies=0; devk=DEV_NULL; devnum=-1;
 if(tens_block == NULL) return TALSH_INVALID_ARGS;
 if(talshTensorIsEmpty(tens_block) == YEP) return TALSH_OBJECT_IS_EMPTY;
 if(valid_device_kind(dev_kind) != YEP) return TALSH_INVALID_ARGS;
 if(dev_kind == DEV_NULL){
  if(dev_id >= 0){
   devnum=talshKindDevId(dev_id,&devk); if(devnum < 0) return TALSH_INVALID_ARGS;
   specific_kind=1; specific_device=1;
  }else{
   specific_kind=0; specific_device=0;
  }
 }else{
  specific_kind=1; devk=dev_kind;
  if(dev_id >= 0){
   devnum=talshFlatDevId(dev_kind,dev_id); if(devnum >= DEV_MAX) return TALSH_INVALID_ARGS;
   specific_device=1; devnum=dev_id;
  }else{
   specific_device=0;
  }
 }
 if(tens_block->ndev > 0){
  if(tens_block->ndev > TALSH_MAX_DEV_PRESENT || tens_block->ndev > tens_block->dev_rsc_len) return TALSH_FAILURE;
  if(tens_block->dev_rsc == NULL || tens_block->data_kind == NULL) return TALSH_FAILURE;
  for(i=0;i<tens_block->ndev;i++){
   j=talshKindDevId(tens_block->dev_rsc[i].dev_id,&m); if(j < 0) return TALSH_FAILURE;
   if((m == devk || specific_kind == 0) && (j == devnum || specific_device == 0)){
    (*ncopies)++; copies[*ncopies]=tens_block->dev_rsc[i].dev_id; data_kinds[*ncopies]=tens_block->data_kind[i];
   }
  }
 }
 return TALSH_SUCCESS;
}

int talshTensorPresence_(const talsh_tens_t * tens_block, int * ncopies, int copies[], int data_kinds[], int dev_kind, int dev_id) //Fortran wrapper
{
 return talshTensorPresence(tens_block, ncopies, copies, data_kinds, dev_kind, dev_id);
}

// TAL-SH task API:
int talshTaskCreate(talsh_task_t ** talsh_task)
/** Creates a clean <talsh_task_t> object on heap. **/
{
 *talsh_task=(talsh_task_t*)malloc(sizeof(talsh_task_t));
 if(*talsh_task == NULL) return TRY_LATER;
 return talshTaskClean(*talsh_task);
}

int talshTaskClean(talsh_task_t * talsh_task)
/** Cleans an undefined (statically allocated) <talsh_task_t> object making it defined-empty.
    Never call this function on value-defined <talsh_task_t> objects. **/
{
 talsh_task->task_p=NULL;
 talsh_task->dev_kind=DEV_NULL;
 talsh_task->data_kind=NO_TYPE;
 talsh_task->data_vol=0.0;
 talsh_task->flops=0.0;
 talsh_task->exec_time=0.0;
 return TALSH_SUCCESS;
}

static int talshTaskConstruct(talsh_task_t * talsh_task, int dev_kind, int data_kind)
/** Constructs a TAL-SH task. It is errorneous to pass undefined <talsh_task_t> objects here!
    At the same time, it is fine to pass a value-defined <talsh_task_t> object here because
    it will be destructed before the new construction. **/
{
 int i,errc;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 errc=TALSH_SUCCESS;
 if(talsh_task == NULL) return TALSH_INVALID_ARGS;
 if(valid_device_kind(dev_kind) != YEP) return TALSH_INVALID_ARGS;
 if(data_kind != NO_TYPE){if(tens_valid_data_kind(data_kind) != YEP) return TALSH_INVALID_ARGS;}
 if(talsh_task->dev_kind != DEV_NULL) errc=talshTaskDestruct(talsh_task); //destruct value-defined tasks first
 if(errc != TALSH_SUCCESS && errc != NOT_CLEAN) return TALSH_FAILURE; if(errc == NOT_CLEAN) talsh_raise_not_clean();
 switch(dev_kind){
  case DEV_HOST:
   i=host_task_create((host_task_t**)(&(talsh_task->task_p)));
   if(i != 0){
    errc=talshTaskClean(talsh_task);
    if(i == TRY_LATER || i == DEVICE_UNABLE){return i;}else{return TALSH_FAILURE;} //overrules previous NOT_CLEAN status
   }
   break;
  case DEV_NVIDIA_GPU:
#ifndef NO_GPU
   i=cuda_task_create((cudaTask_t**)(&(talsh_task->task_p)));
   if(i != 0){
    errc=talshTaskClean(talsh_task);
    if(i == TRY_LATER || i == DEVICE_UNABLE){return i;}else{return TALSH_FAILURE;} //overrules previous NOT_CLEAN status
   }
#else
   return TALSH_NOT_AVAILABLE;
#endif
   break;
  case DEV_INTEL_MIC:
#ifndef NO_PHI
   return TALSH_NOT_IMPLEMENTED; //`Future
#else
   return TALSH_NOT_AVAILABLE;
#endif
   break;
  case DEV_AMD_GPU:
#ifndef NO_AMD
   return TALSH_NOT_IMPLEMENTED; //`Future
#else
   return TALSH_NOT_AVAILABLE;
#endif
   break;
  default:
   return TALSH_INVALID_ARGS;
 }
 talsh_task->dev_kind=dev_kind;
 talsh_task->data_kind=data_kind;
 return errc;
}

int talshTaskDestruct(talsh_task_t * talsh_task)
/** Destructs a TAL-SH task, putting it back into the defined-empty (clean) state. **/
{
 int i,errc;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 errc=TALSH_SUCCESS;
 if(talsh_task == NULL) return TALSH_INVALID_ARGS;
 switch(talsh_task->dev_kind){
  case DEV_HOST:
   if(talsh_task->task_p == NULL) return TALSH_INVALID_ARGS;
   errc=host_task_destroy((host_task_t*)(talsh_task->task_p));
   if(errc != 0 && errc != TRY_LATER && errc != NOT_CLEAN) errc=TALSH_FAILURE;
   break;
  case DEV_NVIDIA_GPU:
#ifndef NO_GPU
   if(talsh_task->task_p == NULL) return TALSH_INVALID_ARGS;
   errc=cuda_task_destroy((cudaTask_t*)(talsh_task->task_p));
   if(errc != 0 && errc != TRY_LATER && errc != NOT_CLEAN) errc=TALSH_FAILURE;
#else
   return TALSH_NOT_AVAILABLE;
#endif
   break;
  case DEV_INTEL_MIC:
#ifndef NO_PHI
   return TALSH_NOT_IMPLEMENTED; //`Future
#else
   return TALSH_NOT_AVAILABLE;
#endif
   break;
  case DEV_AMD_GPU:
#ifndef NO_AMD
   return TALSH_NOT_IMPLEMENTED; //`Future
#else
   return TALSH_NOT_AVAILABLE;
#endif
   break;
  case DEV_NULL: //defined-empty task
   break;
  default:
   return TALSH_INVALID_ARGS;
 }
 i=talshTaskClean(talsh_task);
 return errc;
}

int talshTaskDestroy(talsh_task_t * talsh_task)
/** Completely destroys a <talsh_task_t> object. **/
{
 int errc;

 if(talsh_task == NULL) return TALSH_INVALID_ARGS;
 errc=talshTaskDestruct(talsh_task);
 free(talsh_task);
 return errc;
}

int talshTaskDevId(talsh_task_t * talsh_task, int * dev_kind)
/** Returns either a flat (<dev_kind> is absent) or kind-specific (<dev_kind> is present)
    device id on which the TAL-SH task is scheduled. DEV_NULL on return means an error. **/
{
 int devid,errc;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 if(talsh_task == NULL) return DEV_NULL;
 errc=talshTaskStatus(talsh_task);
 if(errc == TALSH_FAILURE || errc == TALSH_TASK_EMPTY) return DEV_NULL;
 switch(talsh_task->dev_kind){
  case DEV_HOST:
   devid=0; //Host is always single
   break;
  case DEV_NVIDIA_GPU:
#ifndef NO_GPU
   devid=cuda_task_gpu_id((cudaTask_t*)(talsh_task->task_p));
   if(devid < 0) return DEV_NULL;
#else
   return DEV_NULL;
#endif
   break;
  case DEV_INTEL_MIC:
#ifndef NO_PHI
   return DEV_NULL; //`Future
#else
   return DEV_NULL;
#endif
   break;
  case DEV_AMD_GPU:
#ifndef NO_AMD
   return DEV_NULL; //`Future
#else
   return DEV_NULL;
#endif
   break;
  default:
   return DEV_NULL;
 }
 if(devid < 0) return DEV_NULL;
 if(dev_kind != NULL){
  *dev_kind=talsh_task->dev_kind;
 }else{
  devid=talshFlatDevId(talsh_task->dev_kind,devid); //convert to flat device id
  if(devid < 0 || devid >= DEV_MAX) devid=DEV_NULL;
 }
 return devid;
}

int talshTaskDevId_(talsh_task_t * talsh_task, int * dev_kind) //Fortran wrapper
{
 return talshTaskDevId(talsh_task,dev_kind);
}

int talshTaskStatus(talsh_task_t * talsh_task)
/** Returns the current status of the TAL-SH task or an error status. **/
{
 int i,errc;
 host_task_t *host_task_p;
 cudaTask_t *cuda_task_p;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 if(talsh_task == NULL) return TALSH_INVALID_ARGS;
 if(talsh_task->dev_kind == DEV_NULL) return TALSH_TASK_EMPTY;
 if(talsh_task->task_p == NULL) return TALSH_INVALID_ARGS;
 switch(talsh_task->dev_kind){
  case DEV_HOST:
   host_task_p=((host_task_t*)(talsh_task->task_p));
   if(host_task_p->task_error == 0){errc=TALSH_TASK_COMPLETED;}else{errc=TALSH_TASK_ERROR;}
   break;
  case DEV_NVIDIA_GPU:
#ifndef NO_GPU
   cuda_task_p=(cudaTask_t*)(talsh_task->task_p);
   i=cuda_task_status(cuda_task_p);
   switch(i){
    case CUDA_TASK_ERROR: errc=TALSH_TASK_ERROR; break;
    case CUDA_TASK_EMPTY: errc=TALSH_TASK_EMPTY; break;
    case CUDA_TASK_SCHEDULED: errc=TALSH_TASK_SCHEDULED; break;
    case CUDA_TASK_STARTED: errc=TALSH_TASK_STARTED; break;
    case CUDA_TASK_INPUT_THERE: errc=TALSH_TASK_INPUT_READY; break;
    case CUDA_TASK_OUTPUT_THERE: errc=TALSH_TASK_OUTPUT_READY; break;
    case CUDA_TASK_COMPLETED: errc=TALSH_TASK_COMPLETED; break;
    default:
     errc=TALSH_FAILURE;
   }
#else
   return TALSH_NOT_AVAILABLE;
#endif
   break;
  case DEV_INTEL_MIC:
#ifndef NO_PHI
   return TALSH_NOT_IMPLEMENTED; //`Future
#else
   return TALSH_NOT_AVAILABLE;
#endif
   break;
  case DEV_AMD_GPU:
#ifndef NO_AMD
   return TALSH_NOT_IMPLEMENTED; //`Future
#else
   return TALSH_NOT_AVAILABLE;
#endif
   break;
  default:
   return TALSH_INVALID_ARGS;
 }
 return errc;
}

int talshTaskCompleted(talsh_task_t * talsh_task, int * stats, int * ierr)
/** Returns YEP if the TAL-SH has completed, NOPE otherwise.
    The TAL-SH task status will be returned in <stats>. **/
{
 int errc;
 host_task_t *host_task_p;
 cudaTask_t *cuda_task_p;

 errc=NOPE;
 if(talsh_on == 0){*ierr=TALSH_NOT_INITIALIZED; return errc;}
 if(ierr == NULL) return TALSH_INVALID_ARGS;
 if(talsh_task == NULL || stats == NULL){*ierr=TALSH_INVALID_ARGS; return errc;}
 if(talsh_task->task_p == NULL){*ierr=TALSH_OBJECT_IS_EMPTY; return errc;}
 *ierr=TALSH_SUCCESS;
 switch(talsh_task->dev_kind){
  case DEV_HOST:
   host_task_p=((host_task_t*)(talsh_task->task_p));
   if(host_task_p->task_error == 0){*stats=TALSH_TASK_COMPLETED;}else{*stats=TALSH_TASK_ERROR;}
   return YEP;
  case DEV_NVIDIA_GPU:
#ifndef NO_GPU
   cuda_task_p=(cudaTask_t*)(talsh_task->task_p);
   *stats=cuda_task_completed(cuda_task_p);
   switch(*stats){
    case CUDA_TASK_ERROR: *stats=TALSH_TASK_ERROR; errc=YEP; break;
    case CUDA_TASK_EMPTY: *stats=TALSH_TASK_EMPTY; break;
    case CUDA_TASK_SCHEDULED: *stats=TALSH_TASK_SCHEDULED; break;
    case CUDA_TASK_STARTED: *stats=TALSH_TASK_STARTED; break;
    case CUDA_TASK_INPUT_THERE: *stats=TALSH_TASK_INPUT_READY; break;
    case CUDA_TASK_OUTPUT_THERE: *stats=TALSH_TASK_OUTPUT_READY; break;
    case CUDA_TASK_COMPLETED: *stats=TALSH_TASK_COMPLETED; errc=YEP; break;
    default:
     *stats=TALSH_FAILURE; *ierr=TALSH_FAILURE;
   }
#else
   *ierr=TALSH_NOT_AVAILABLE;
#endif
   return errc;
  case DEV_INTEL_MIC:
#ifndef NO_PHI
   *ierr=TALSH_NOT_IMPLEMENTED; //`Future
#else
   *ierr=TALSH_NOT_AVAILABLE;
#endif
   return errc;
  case DEV_AMD_GPU:
#ifndef NO_AMD
   *ierr=TALSH_NOT_IMPLEMENTED; //`Future
#else
   *ierr=TALSH_NOT_AVAILABLE;
#endif
   return errc;
  default:
   *ierr=TALSH_INVALID_ARGS;
   return errc;
 }
 return errc;
}

int talshTaskWait(talsh_task_t * talsh_task, int * stats)
/** Returns upon completion of a TAL-SH task. **/
{
 int errc;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 if(talsh_task == NULL || stats == NULL) return TALSH_INVALID_ARGS;
 errc=TALSH_SUCCESS;
 while(talshTaskCompleted(talsh_task,stats,&errc) == NOPE){if(errc != TALSH_SUCCESS) break;};
 return errc;
}

int talshTasksWait(int ntasks, talsh_task_t talsh_tasks[], int stats[])
/** Returns upon completion of a number of TAL-SH tasks. **/
{
 int i,tc,errc;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 if(ntasks <= 0 || talsh_tasks == NULL || stats == NULL) return TALSH_INVALID_ARGS;
 for(i=0;i<ntasks;++i) stats[i]=TALSH_TASK_EMPTY;
 tc=ntasks; errc=TALSH_SUCCESS;
 while(tc > 0){
  for(i=0;i<ntasks;++i){
   if(talsh_tasks[i].task_p == NULL || talsh_tasks[i].dev_kind == DEV_NULL) return TALSH_OBJECT_IS_EMPTY;
   if(stats[i] == TALSH_TASK_EMPTY){if(talshTaskCompleted(&(talsh_tasks[i]),&(stats[i]),&errc) == YEP) --tc;}
   if(errc != TALSH_SUCCESS) return TALSH_FAILURE;
  };
 };
 return TALSH_SUCCESS;
}

int talshTaskTime(talsh_task_t * talsh_task, double * total, double * comput, double * input, double * output)
/** Returns the timing information for a given TAL-SH task. **/
{
 int sts,errc;
 float tot_tm,in_tm,out_tm,comp_tm;
 cudaTask_t *cuda_task_p;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 if(talsh_task == NULL || total == NULL) return TALSH_INVALID_ARGS;
 if(talsh_task->task_p == NULL) return TALSH_OBJECT_IS_EMPTY;
 if(talshTaskCompleted(talsh_task,&sts,&errc) == NOPE){
  if(errc != TALSH_SUCCESS) return TALSH_FAILURE;
  return TALSH_IN_PROGRESS;
 }
 switch(talsh_task->dev_kind){
  case DEV_HOST:
   tot_tm=(float)talsh_task->exec_time; in_tm=-1.0f; out_tm=-1.0f; comp_tm=-1.0f;
   if(tot_tm < 0.0f) errc=TALSH_FAILURE;
  case DEV_NVIDIA_GPU:
#ifndef NO_GPU
   cuda_task_p=(cudaTask_t*)(talsh_task->task_p);
   tot_tm=cuda_task_time(cuda_task_p,&in_tm,&out_tm,&comp_tm);
   if(tot_tm < 0.0f) errc=TALSH_FAILURE;
#else
   return TALSH_NOT_AVAILABLE;
#endif
   break;
  case DEV_INTEL_MIC:
#ifndef NO_PHI
   return TALSH_NOT_IMPLEMENTED; //`Future
#else
   return TALSH_NOT_AVAILABLE;
#endif
   break;
  case DEV_AMD_GPU:
#ifndef NO_AMD
   return TALSH_NOT_IMPLEMENTED; //`Future
#else
   return TALSH_NOT_AVAILABLE;
#endif
   break;
  default:
   return TALSH_INVALID_ARGS;
 }
 *total=(double)tot_tm;
 if(comput != NULL) *comput=(double)comp_tm;
 if(input != NULL) *input=(double)in_tm;
 if(output != NULL) *output=(double)out_tm;
 return errc;
}

int talshTaskTime_(talsh_task_t * talsh_task, double * total, double * comput, double * input, double * output) //Fortran wrapper
{
 return talshTaskTime(talsh_task,total,comput,input,output);
}

// TAL-SH tensor operations API:
int talshTensorPlace(talsh_tens_t * tens, int dev_id, int dev_kind, int copy_ctrl, talsh_task_t * talsh_task)
/** Places a tensor block on a specific device. **/
{
 int errc;

 errc=TALSH_SUCCESS;
 return errc;
}

int talshTensorPlace_(talsh_tens_t * tens, int dev_id, int dev_kind, int copy_ctrl, talsh_task_t * talsh_task) //Fortran wrapper
{
 return talshTensorPlace(tens,dev_id,dev_kind,copy_ctrl,talsh_task);
}

int talshTensorDiscard(talsh_tens_t * tens, int dev_id, int dev_kind)
/** Discards a tensor block on a specific device. **/
{
 int i,j,k,errc,devid;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 if(tens == NULL) return TALSH_INVALID_ARGS;
 if(talshTensorIsEmpty(tens) == YEP) return TALSH_OBJECT_IS_EMPTY;
 if(tens->ndev <= 0 || tens->ndev > tens->dev_rsc_len) return TALSH_FAILURE;
 errc=TALSH_SUCCESS;
 if(dev_kind == DEV_NULL){devid=dev_id;}else{devid=talshFlatDevId(dev_kind,dev_id);}
 if(devid < 0 || devid >= DEV_MAX) return TALSH_INVALID_ARGS;
 k=0;
 for(i=0;i<tens->ndev;i++){
  if(tens->dev_rsc[i].dev_id == devid){
   j=tensDevRsc_release_all(&(tens->dev_rsc[i]));
   if(j != 0 && errc != TALSH_FAILURE){if(j == NOT_CLEAN){errc=NOT_CLEAN;}else{errc=TALSH_FAILURE;}}
  }else{
   if(i != k){tens->dev_rsc[k]=tens->dev_rsc[i]; tens->data_kind[k]=tens->data_kind[i];}
   k++;
  }
 }
 tens->ndev=k;
 return errc;
}

int talshTensorDiscard_(talsh_tens_t * tens, int dev_id, int dev_kind) //Fortran wrapper
{
 return talshTensorDiscard(tens,dev_id,dev_kind);
}

int talshTensorContract(const char * cptrn,        //in: C-string: symbolic contraction pattern, e.g. "D(a,b,c,d)+=L(c,i,j,a)*R(b,j,d,i)"
                        talsh_tens_t * dtens,      //inout: destination tensor block
                        talsh_tens_t * ltens,      //inout: left source tensor block
                        talsh_tens_t * rtens,      //inout: right source tensor block
                        int copy_ctrl,             //in: copy control (COPY_XXX), defaults to COPY_MTT
                        double scale_real,         //in: scaling value (real part), defaults to 1
                        double scale_imag,         //in: scaling value (imaginary part), defaults to 0
                        int dev_id,                //in: device id (flat or kind-specific)
                        int dev_kind,              //in: device kind (if present, <dev_id> is kind-specific)
                        talsh_task_t * talsh_task) //inout: TAL-SH task (must be clean)
/** Tensor contraction. **/
{
 int errc;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 errc=TALSH_SUCCESS;

 return errc;
}

int talshTensorContract_(const char * cptrn, talsh_tens_t * dtens, talsh_tens_t * ltens, talsh_tens_t * rtens, int copy_ctrl,
                         double scale_real, double scale_imag, int dev_id, int dev_kind, talsh_task_t * talsh_task) //Fortran wrapper
{
 return talshTensorContract(cptrn,dtens,ltens,rtens,copy_ctrl,scale_real,scale_imag,dev_id,dev_kind,talsh_task);
}
