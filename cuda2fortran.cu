//This file contains C wrappers for CUDA run-time functions to be called from Fortran.
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

#ifndef NO_GPU
//Protect the C function names from name mangling:
#ifdef __cplusplus
extern "C"{
#endif
 void cudagetdevicecount(int* count, int* err_code);
 void cudasetdevice(int device, int* err_code);
 void cudagetdeviceproperties(int device, size_t *totalGlobalMem_, size_t *sharedMemPerBlock_,
       int *regsPerBlock_, int *warpSize_, int *maxThreadsPerBlock_, int *maxThreadsDim_, int *maxGridSize_,
       int *clockRate_, size_t *totalConstMem_, int *major_, int *minor_, int *deviceOverlap_,
       int *multiProcessorCount_, int *concurrentKernels_, int *ECCEnabled_, int *asyncEngineCount_,
       int *memoryClockRate_, int *memoryBusWidth_, int *maxThreadsPerMultiProcessor_, int* err_code);
 void cudadevicesynchronize(int *err_code);
#ifdef __cplusplus
}
#endif
#endif
#ifdef __cplusplus
extern "C"{
#endif
 void print_c_ptr(void *c_ptr);
#ifdef __cplusplus
}
#endif
//-------------------------------------------------------------------
#ifndef NO_GPU
//C Wrappers (called from Fortran to invoke CUDA run-time functions):
__host__ void cudagetdevicecount(int* count, int* err_code){
 cudaError_t err; const char* err_msg;
 *err_code=0;
 err=cudaGetDeviceCount(count); if(err!=cudaSuccess){
  err_msg=cudaGetErrorString(err);
  printf("#ERROR(cudagetdevicecount): %s \n",err_msg);
  *err_code=1;
 };
 return;
}

__host__ void cudasetdevice(int device, int* err_code){
 cudaError_t err; const char* err_msg;
 *err_code=0;
 err=cudaSetDevice(device); if(err!=cudaSuccess){
  err_msg=cudaGetErrorString(err);
  printf("#ERROR(cudasetdevice): %s \n",err_msg);
  *err_code=1;
 };
 return;
}

__host__ void cudagetdeviceproperties(int device, size_t *totalGlobalMem_, size_t *sharedMemPerBlock_,
               int *regsPerBlock_, int *warpSize_, int *maxThreadsPerBlock_, int *maxThreadsDim_, int *maxGridSize_,
               int *clockRate_, size_t *totalConstMem_, int *major_, int *minor_, int *deviceOverlap_,
               int *multiProcessorCount_, int *concurrentKernels_, int *ECCEnabled_, int *asyncEngineCount_,
               int *memoryClockRate_, int *memoryBusWidth_, int *maxThreadsPerMultiProcessor_, int* err_code){
 cudaError_t err; const char* err_msg; cudaDeviceProp prop;
 *err_code=0;
 err=cudaGetDeviceProperties(&prop,device);
 if(err!=cudaSuccess){
  err_msg=cudaGetErrorString(err);
  printf("#ERROR(cudagetdeviceproperties): %s \n",err_msg);
  *err_code=1;
 }else{
  *totalGlobalMem_=prop.totalGlobalMem;
  *sharedMemPerBlock_=prop.sharedMemPerBlock;
  *regsPerBlock_=prop.regsPerBlock;
  *warpSize_=prop.warpSize;
  *maxThreadsPerBlock_=prop.maxThreadsPerBlock;
  maxThreadsDim_[0]=prop.maxThreadsDim[0]; maxThreadsDim_[1]=prop.maxThreadsDim[1]; maxThreadsDim_[2]=prop.maxThreadsDim[2];
  maxGridSize_[0]=prop.maxGridSize[0]; maxGridSize_[1]=prop.maxGridSize[1]; maxGridSize_[2]=prop.maxGridSize[2];
  *clockRate_=prop.clockRate;
  *totalConstMem_=prop.totalConstMem;
  *major_=prop.major; *minor_=prop.minor;
  *deviceOverlap_=prop.deviceOverlap;
  *multiProcessorCount_=prop.multiProcessorCount;
  *concurrentKernels_=prop.concurrentKernels;
  *ECCEnabled_=prop.ECCEnabled;
  *asyncEngineCount_=prop.asyncEngineCount;
  *memoryClockRate_=prop.memoryClockRate;
  *memoryBusWidth_=prop.memoryBusWidth;
  *maxThreadsPerMultiProcessor_=prop.maxThreadsPerMultiProcessor;
 };
 return;
}
__host__ void cudadevicesynchronize(int *err_code)
{
 *err_code=0;
 cudaError_t err=cudaDeviceSynchronize(); if(err != cudaSuccess){*err_code=1;}
 return;
}
#endif
//---------------------------------------------------------------
//Auxiliary functions:
int string_len(const char* str){
 const int max_string_len=2147483647;
 int i;
 for(i=0;i<max_string_len;i++){if(str[i]==0) break;};
 return i;
}
void print_c_ptr(void *c_ptr){
 printf("%p ",c_ptr);
 return;
}
