/** ExaTensor::TAL-SH: Device-unified user-level API.
REVISION: 2016/04/06

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
static int talsh_gpu[MAX_GPUS_PER_NODE]={DEV_OFF}; //current GPU status
static int talsh_mic[MAX_MICS_PER_NODE]={DEV_OFF}; //current MIC status
static int talsh_amd[MAX_AMDS_PER_NODE]={DEV_OFF}; //current AMD status

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
  return TALSH_NOT_IMPLEMENTED;
 }
//AMD GPU accelerators:
 if(namds > 0){
  printf("#FATAL(TALSH::talshInit): AMD GPU is not supported yet!");
  return TALSH_NOT_IMPLEMENTED;
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
   return TALSH_NOT_IMPLEMENTED; //`Implement
  case DEV_AMD_GPU:
   return TALSH_NOT_IMPLEMENTED; //`Implement
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
   rc=TALSH_NOT_IMPLEMENTED;
   break;
  case DEV_NVIDIA_GPU:
   rc=gpu_print_stats(dev_id);
   break;
  case DEV_INTEL_MIC:
   rc=TALSH_NOT_IMPLEMENTED;
   break;
  case DEV_AMD_GPU:
   rc=TALSH_NOT_IMPLEMENTED;
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
 (*tens_block)=NULL;
 (*tens_block)=(talsh_tens_t*)malloc(sizeof(talsh_tens_t));
 if(*tens_block == NULL) return TRY_LATER;
 errc=talshTensorClean(*tens_block);
 return errc;
}

int talshTensorClean(talsh_tens_t * tens_block)
/** Clean a tensor block. **/
{
 if(tens_block == NULL) return TALSH_INVALID_ARGS;
 tens_block->shape_p=NULL;
 tens_block->ndev=0;
 tens_block->last_write=DEV_NULL;
 tens_block->dev_rsc_len=0;
 tens_block->dev_rsc=NULL;
 tens_block->tensF=NULL;
 tens_block->tensC=NULL;
 return TALSH_SUCCESS;
}

int talshTensorConstruct(talsh_tens_t * tens_block,     //inout: constructed tensor block (must be empty on entrance)
                         int data_type,                 //in: data type: {R4,R8,C4,C8,NO_TYPE}
                         int tens_rank,                 //in: tensor block rank (number of dimensions)
                         const int tens_dims[],         //in: tensor block dimension extents
                         int dev_id,                    //in: flat device ID on which the tensor block will reside
                         void * ext_mem,                //in: pointer to externally provided memory for tensor elements
                         int in_hab,                    //in: if >=0, <ext_mem> points to the HAB entry #<in_hab>
                         talsh_tens_init_i init_method, //in: user-defined initialization method (function pointer)
                         double init_val_real,          //in: initialization value (real part), defaults to 0.0
                         double init_val_imag)          //in: initialization value (imaginary part), defaults to 0.0
/** Constructs a tensor block: {0: success; TRY_LATER: currently no enough memory available}.
    If <data_type> == NO_TYPE, the tensor body is not allocated (only the tensor shape). **/
{
 int dev_num,dksize,errc,already_allocated;
 size_t tvol,tsize;

 errc=TALSH_SUCCESS;
 //Check arguments:
 if(tens_block == NULL) return TALSH_INVALID_ARGS; //tensor block object must have been preallocated
 if(tens_block->shape_p != NULL) return TALSH_ALREADY_INITIALIZED; //tensor block is not empty (destruct it first)
 if(tens_valid_data_kind(data_type,&dksize) != YEP) return TALSH_INVALID_ARGS; //unknown data type (NO_TYPE is a valid type here)
 dev_num=decode_device_id(dev_id); if(dev_num < 0) return TALSH_INVALID_ARGS; //invalid device id
 already_allocated=0; if(ext_mem != NULL) already_allocated=1; //check whether an external memory space is provided for the tensor body
 if(in_hab >= 0){ //HAB entry number must be accompanied with the external memory pointer <ext_mem>
  if(already_allocated){already_allocated=2;}else{return TALSH_INVALID_ARGS;};
 }else{
  in_hab=-1; //no HAB entry in use
 }
 //Tensor shape:
 errc=tensShape_create(&(tens_block->shape_p)); if(errc == TRY_LATER || errc == DEVICE_UNABLE) return errc;
 if(errc != 0 || tens_block->shape_p == NULL) return TALSH_FAILURE;
 errc=tensShape_construct(tens_block->shape_p,NOPE,tens_rank,tens_dims); //NOPE = not pinned
 if(errc != 0 && errc != TRY_LATER && errc != DEVICE_UNABLE) errc=TALSH_FAILURE;
 if(errc != 0){free(tens_block->shape_p); tens_block->shape_p=NULL; return errc;};
 //Device resource storage:
 if(tens_block->dev_rsc_len == 0 || tens_block->dev_rsc == NULL){
  tens_block->dev_rsc=(talsh_dev_rsc_t*)malloc(TALSH_MAX_DEV_PRESENT*sizeof(talsh_dev_rsc_t));
  if(tens_block->dev_rsc != NULL){
   tens_block->dev_rsc_len=TALSH_MAX_DEV_PRESENT; tens_block->ndev=0; tens_block->last_write=DEV_NULL;
  }else{
   free(tens_block->shape_p); tens_block->shape_p=NULL; tens_block->dev_rsc=NULL; tens_block->dev_rsc_len=0;
   return TRY_LATER;
  }
 }else{
  free(tens_block->shape_p); tens_block->shape_p=NULL;
  return TALSH_INVALID_ARGS;
 }
 //Tensor body:
 if(already_allocated != 0){
  errc=tensDevRsc_attach_mem(&(tens_block->dev_rsc[0]),dev_id,ext_mem,in_hab);
  if(errc != 0){
   free(tens_block->shape_p); tens_block->shape_p=NULL;
   free(tens_block->dev_rsc); tens_block->dev_rsc=NULL; tens_block->dev_rsc_len=0;
   return TALSH_FAILURE;
  }
  tens_block->ndev=1;
 }else{
  if(data_type != NO_TYPE){
   
  }
 }
 return errc;
}

int talshTensorDestruct(talsh_tens_t * tens_block) //in: non-NULL pointer to a tensor block (empty tensor block on exit)
/** Destructs a tensor block and sets its status to empty. **/
{
 int errc;

 errc=TALSH_SUCCESS;
 if(tens_block == NULL) return TALSH_INVALID_ARGS;
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
