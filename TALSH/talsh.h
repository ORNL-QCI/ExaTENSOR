/** ExaTensor::TAL-SH: Device-unified user-level API header.
REVISION: 2016/04/15

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
------------------------------------------------------------------------
**/

#ifndef _TALSH_H
#define _TALSH_H

#include "tensor_algebra.h"

//TAL-SH PARAMETERS:
#define TALSH_MAX_ACTIVE_TASKS 4096 //max number of active tasks on all devices on a node
#define TALSH_MAX_DEV_PRESENT 16 //max number of on-node devices the tensor block can be present on

//TAL-SH ERROR CODES:
#define TALSH_SUCCESS 0
#define TALSH_FAILURE -666
#define TALSH_NOT_AVAILABLE -888
#define TALSH_NOT_IMPLEMENTED -999
#define TALSH_NOT_INITIALIZED 1000000
#define TALSH_ALREADY_INITIALIZED 1000001
#define TALSH_INVALID_ARGS 1000002
#define TALSH_INTEGER_OVERFLOW 1000003

//TAL-SH DATA TYPES:
// Interoperable tensor block:
typedef struct{
 talsh_tens_shape_t * shape_p; //shape of the tensor block
 int ndev;                     //number of devices the tensor block resides on
 int last_write;               //flat device id where the last write happened, -1 means coherence on all devices where the tensor block resides
 int dev_rsc_len;              //capacity of dev_rsc[]: ndev <= dev_rsc_len
 talsh_dev_rsc_t * dev_rsc;    //list of device resources occupied by the tensor block on each device
 void * tensF;                 //pointer to Fortran <tensor_block_t> (CPU,Phi)
 void * tensC;                 //pointer to C <tensBlck_t> (Nvidia GPU)
} talsh_tens_t;

// Interoperable TAL-SH task handle:
typedef struct{
 void * task_p;    //pointer to the corresponding task object
 int dev_kind;     //device kind (DEV_NULL: uninitalized)
 double flops;     //number of floating point operations
 double exec_time; //execution time in seconds
} talsh_task_t;

//EXPORTED FUNCTIONS:
#ifdef __cplusplus
extern "C"{
#endif
// TAL-SH control API:
//  Initialize TAL-SH:
 int talshInit(size_t * host_buf_size,
               int * host_arg_max,
               int ngpus,
               int gpu_list[],
               int nmics,
               int mic_list[],
               int namds,
               int amd_list[]);
//  Shutdown TAL-SH:
 int talshShutdown();
//  Get the flat device Id:
 int talshFlatDevId(int dev_kind,
                    int dev_num);
//  Get the kind-specific device Id:
 int talshKindDevId(int dev_id,
                    int * dev_kind);
//  Query the state of a device:
 int talshDeviceState(int dev_num,
                      int dev_kind = DEV_NULL);
 int talshDeviceState_(int dev_num, int dev_kind);
//  Find the least busy device:
 int talshDeviceBusyLeast(int dev_kind = DEV_NULL);
 int talshDeviceBusyLeast_(int dev_kind);
//  Print TAL-SH statistics for specific devices:
 int talshStats(int dev_id = -1,
                int dev_kind = DEV_NULL);
 int talshStats_(int dev_id, int dev_kind);
// TAL-SH tensor block API:
//  Create an empty tensor block:
 int talshTensorCreate(talsh_tens_t ** tens_block);
//  Clean a tensor block (default constructor):
 int talshTensorClean(talsh_tens_t * tens_block);
//  Check whether a tensor block is empty:
 int talshTensorIsEmpty(const talsh_tens_t * tens_block);
//  Construct a tensor block:
 int talshTensorConstruct(talsh_tens_t * tens_block,
                          int data_type,
                          int tens_rank,
                          int tens_dims[],
                          int dev_id = 0,
                          void * ext_mem = NULL,
                          int in_hab = -1,
                          talsh_tens_init_i init_method = NULL,
                          double init_val_real = 0.0,
                          double init_val_imag = 0.0);
 int talshTensorConstruct_(talsh_tens_t * tens_block, int data_type, int tens_rank, int tens_dims[], int dev_id,
                           void * ext_mem, int in_hab, talsh_tens_init_i init_method, double init_val_real, double init_val_imag);
//  Destruct a tensor block:
 int talshTensorDestruct(talsh_tens_t * tens_block);
//  Destroy a tensor block:
 int talshTensorDestroy(talsh_tens_t * tens_block);
//  Get the tensor block volume (number of elements):
 size_t talshTensorVolume(const talsh_tens_t * tens_block);

#ifdef __cplusplus
}
#endif

//HEADER GUARD
#endif
