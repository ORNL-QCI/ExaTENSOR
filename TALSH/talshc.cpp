/** ExaTensor::TAL-SH: Device-unified user-level API. **/

#include <stdio.h>
#include <stdlib.h>

#ifndef NO_GPU
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#include "tensor_algebra.h"
#include "talsh.h"
