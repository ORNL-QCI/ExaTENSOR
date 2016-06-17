const int TILEDIM = 32;
const int TILEROWS = 8;

// Arguments for transposeTensorKernel
// Enough to support tensors of rank 200 or so.
const int transposeArgSize = 2048 - 5 - 2*2;
__constant__ int transposeArg[transposeArgSize];

template <typename T>
__global__ void transposeTensorKernelArg(
  const int volMbar, const int volMm, const int sizeMmk, const int sizeMbar,
  const int cuDimMk, const int cuDimMm,
  const int3* __restrict__ gl_MbarOut,
  const T* __restrict__ dataIn, T* __restrict__ dataOut) {

  // Shared memory
  __shared__ T shTile[TILEDIM][TILEDIM+1];

  const int warpLane = threadIdx.x & (warpSize - 1);
  int3 MbarOut;
  if (warpLane < sizeMbar) {
    MbarOut = gl_MbarOut[warpLane];
  }

  int* dimMmkIn  = &transposeArg[0];
  int* dimMmkOut = &transposeArg[sizeMmk];

  const int xin = blockIdx.x * TILEDIM + threadIdx.x;
  const int yin = blockIdx.y * TILEDIM + threadIdx.y;

  const int xout = blockIdx.x * TILEDIM + threadIdx.y;
  const int yout = blockIdx.y * TILEDIM + threadIdx.x;

  for (int blockz=blockIdx.z;blockz < volMbar;blockz += blockDim.z*gridDim.z)
  {

    // Read from global memory
    {
      int pos0 = xin + yin*cuDimMk + blockz*volMm;

      __syncthreads();

      // Read data into shared memory tile
      for (int j=0;j < TILEDIM;j += TILEROWS) {
        int pos = pos0 + j*cuDimMk;
        if ((xin < dimMmkIn[0]) && (yin + j < dimMmkIn[1])) {
          shTile[threadIdx.y + j][threadIdx.x] = dataIn[pos];
        }
      }
    }

    // Write to global memory
    {

      int pos0 = 0;
      if (warpLane < sizeMbar) {
        int z = blockz/MbarOut.x;
        pos0 = (z % MbarOut.y)*MbarOut.z;
      }
#pragma unroll
      for (int i=16;i >= 1;i/=2) {
        pos0 += __shfl_xor(pos0, i);
      }
      pos0 += yout + xout*cuDimMm;

      __syncthreads();

      // Write data into global memory
      for (int j=0;j < TILEDIM;j += TILEROWS) {
        int pos = pos0 + j*cuDimMm;
        if ((yout < dimMmkOut[0]) && (xout + j < dimMmkOut[1])) {
          dataOut[pos] = shTile[threadIdx.x][threadIdx.y + j];
        }
      }
    }

  }

}
