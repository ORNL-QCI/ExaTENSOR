/** ExaTensor::TAL-SH: Device-unified user-level API.
REVISION: 2016/04/18

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

#include "talsh.h"

static int talsh_on=0;    //TAL-SH initialization flag (1:initalized; 0:not)
static int talsh_gpu_beg; //first Nvidia GPU in the range `Obsolete
static int talsh_gpu_end; //last Nvidia GPU in the range `Obsolete

static int talsh_cpu=DEV_OFF;
static int talsh_gpu[MAX_GPUS_PER_NODE]={DEV_OFF}; //current GPU status: {DEV_OFF,DEV_ON,DEV_ON_BLAS}
static int talsh_mic[MAX_MICS_PER_NODE]={DEV_OFF}; //current MIC status: {DEV_OFF,DEV_ON,DEV_ON_BLAS}
static int talsh_amd[MAX_AMDS_PER_NODE]={DEV_OFF}; //current AMD status: {DEV_OFF,DEV_ON,DEV_ON_BLAS}

static talsh_task_t talsh_tasks[TALSH_MAX_ACTIVE_TASKS]; //reusable TAL-SH tasks

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
 if(ngpus > 0){
  if(ngpus > MAX_GPUS_PER_NODE) return TALSH_INVALID_ARGS;
  gpu_beg=gpu_list[0]; gpu_end=gpu_list[ngpus-1];
  if(gpu_beg < 0 || gpu_beg >= MAX_GPUS_PER_NODE) return TALSH_INVALID_ARGS;
  if(gpu_end < 0 || gpu_end >= MAX_GPUS_PER_NODE) return TALSH_INVALID_ARGS;
  for(i=1;i<ngpus;i++){
   if(gpu_list[i] != gpu_list[i-1]+1){
    printf("#FATAL(TALSH::talshInit): The current version only supports consecutive GPU ranges!");
    return TALSH_FAILURE;
   }
  }
  errc=arg_buf_allocate(host_buf_size,host_arg_max,gpu_beg,gpu_end); if(errc) return TALSH_FAILURE;
  talsh_gpu_beg=gpu_beg; talsh_gpu_end=gpu_end;
  for(i=0;i<ngpus;i++){
   j=gpu_list[i]; if(j < 0 || j >= MAX_GPUS_PER_NODE) return TALSH_INVALID_ARGS;
   talsh_gpu[j]=gpu_is_mine(j);
  }
 }else{
  talsh_gpu_beg=0; talsh_gpu_end=-1;
 }
//Intel Xeon Phi accelerators:
 if(nmics > 0){
  printf("#FATAL(TALSH::talshInit): Intel Xeon Phi is not fully supported yet!");
  return TALSH_NOT_IMPLEMENTED; //`Future
 }
//AMD GPU accelerators:
 if(namds > 0){
  printf("#FATAL(TALSH::talshInit): AMD GPU is not supported yet!");
  return TALSH_NOT_IMPLEMENTED; //`Future
 }
 talsh_on=1;
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
   i=gpu_busy_least();
   if(i < 0 || i >= MAX_GPUS_PER_NODE) return TALSH_FAILURE;
   return i;
  case DEV_INTEL_MIC:
   return TALSH_NOT_IMPLEMENTED; //`Implement in future
  case DEV_AMD_GPU:
   return TALSH_NOT_IMPLEMENTED; //`Implement in future
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
   rc=gpu_print_stats(dev_id);
   break;
  case DEV_INTEL_MIC:
   rc=TALSH_NOT_IMPLEMENTED; //`Implement in future
   break;
  case DEV_AMD_GPU:
   rc=TALSH_NOT_IMPLEMENTED; //`Implement in future
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
 int errc;
 *tens_block=(talsh_tens_t*)malloc(sizeof(talsh_tens_t));
 if(*tens_block == NULL) return TRY_LATER;
 errc=talshTensorClean(*tens_block);
 return errc;
}

int talshTensorClean(talsh_tens_t * tens_block)
/** Cleans a tensor block (default ctor). **/
{
 if(tens_block == NULL) return TALSH_INVALID_ARGS;
 tens_block->shape_p=NULL;
 tens_block->dev_rsc=NULL;
 tens_block->data_type=NULL;
 tens_block->dev_rsc_len=0;
 tens_block->ndev=0;
 tens_block->last_write=DEV_NULL;
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
                         int data_type,                 //in: data type: {R4,R8,C4,C8,NO_TYPE}
                         int tens_rank,                 //in: tensor block rank (number of dimensions)
                         const int tens_dims[],         //in: tensor block dimension extents
                         int dev_id,                    //in: flat device ID on which the tensor block will reside
                         void * ext_mem,                //in: pointer to externally provided memory for tensor elements
                         int in_hab,                    //in: if >=0, a non-NULL <ext_mem> points to the HAB entry #<in_hab>
                         talsh_tens_init_i init_method, //in: user-defined initialization method (function pointer)
                         double init_val_real,          //in: initialization value (real part), defaults to 0.0
                         double init_val_imag)          //in: initialization value (imaginary part), defaults to 0.0
/** Constructs a tensor block: {0: success; TRY_LATER: no enough free memory available; DEVICE_UNABLE: device unable}.
    If <data_type> == NO_TYPE, the tensor body will not be allocated (only the tensor shape),
    unless an external storage is provided (<ext_mem>). In case the tensor body storage space
    is provided externally (<ext_mem> != NULL), the initialization step is skipped. In other cases,
    unless <data_type>=NO_TYPE, the newly allocated tensor body will be initialized by a user-defined
    method, or, if no method is provided (NULL), by a user-defined value, which defaults to zero.
    If the tensor body initialization failed, a status NOT_CLEAN is returned but
    the tensor block is ready for use (except its body value is undefined). **/
{
 int i,dev_num,dev_kind,dksize,errc,already_allocated,use_hab;
 size_t tvol,tsize;
 float fval;
 float *fp;
 double *dp;

 if(talsh_on == 0) return TALSH_NOT_INITIALIZED;
 errc=TALSH_SUCCESS;
 //Check arguments:
 if(tens_block == NULL) return TALSH_INVALID_ARGS; //tensor block must have been preallocated
 if(talshTensorIsEmpty(tens_block) != YEP) return TALSH_OBJECT_NOT_EMPTY; //tensor block is not empty (destruct it first)
 if(tens_valid_data_kind(data_type,&dksize) != YEP) return TALSH_INVALID_ARGS; //unknown data type (NO_TYPE is a valid type here)
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
 if(tens_block->dev_rsc_len == 0 && tens_block->dev_rsc == NULL){
  tens_block->dev_rsc=(talsh_dev_rsc_t*)malloc(TALSH_MAX_DEV_PRESENT*sizeof(talsh_dev_rsc_t));
  if(tens_block->dev_rsc != NULL){
   tens_block->dev_rsc_len=TALSH_MAX_DEV_PRESENT; tens_block->ndev=0; tens_block->last_write=DEV_NULL;
   tens_block->data_type=(int*)malloc(TALSH_MAX_DEV_PRESENT*sizeof(int));
   if(tens_block->data_type == NULL){i=talshTensorDestruct(tens_block); return TRY_LATER;}
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
  tens_block->data_type[0]=data_type; tens_block->ndev=1; //present on one device, even if NO_TYPE
 }else{ //tensor body storage needs to be allocated here (will also be initialized), unless NO_TYPE
  if(data_type != NO_TYPE){
   tvol=talshTensorVolume(tens_block);
   if(tvol > 0){
    tsize=tvol*dksize;
    if(tsize <= 0){i=talshTensorDestruct(tens_block); return TALSH_INTEGER_OVERFLOW;}
    errc=tensDevRsc_allocate_mem(&(tens_block->dev_rsc[0]),dev_id,tsize,use_hab);
    if(errc != 0 && errc != TRY_LATER && errc != DEVICE_UNABLE) errc=TALSH_FAILURE;
    if(errc != 0){i=talshTensorDestruct(tens_block); return errc;}
    tens_block->data_type[0]=data_type; tens_block->ndev=1; //present on one device
   }else{
    i=talshTensorDestruct(tens_block);
    return TALSH_FAILURE;
   }
   //Initialization:
   if(tens_block->ndev > 0){
    if(dev_kind == DEV_HOST){ //`Currently supported only on Host
     if(init_method != NULL){
      init_method(tens_block->dev_rsc[0].gmem_p,data_type,tens_rank,tens_dims,&errc);
      if(errc) errc=NOT_CLEAN; //initialization failed, tensor block value is undefined, but one may continue
     }else{
      switch(data_type){
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
        return NOT_CLEAN; //`Enable initialization for complex data types C4 and C8
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

int talshTensorConstruct_(talsh_tens_t * tens_block, int data_type, int tens_rank, int tens_dims[], int dev_id, //Fortran wrapper
                          void * ext_mem, int in_hab, talsh_tens_init_i init_method, double init_val_real, double init_val_imag)
{
 return talshTensorConstruct_(tens_block, data_type, tens_rank, tens_dims, dev_id,
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
 if(tens_block->data_type != NULL){free(tens_block->data_type); tens_block->data_type=NULL;}
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
    passed here must either be empty defined or value defined. It is errorneous
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
