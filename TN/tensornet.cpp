/** ExaTENSOR header for the tensor network infrastructure

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/12/07

!Copyright (C) 2014-2017 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2017 Oak Ridge National Laboratory (UT-Battelle)

!This file is part of ExaTensor.

!ExaTensor is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published
!by the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.

!ExaTensor is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!GNU Lesser General Public License for more details.

!You should have received a copy of the GNU Lesser General Public License
!along with ExaTensor. If not, see <http://www.gnu.org/licenses/>.

**/

#include "tensornet.hpp"

namespace exatensor {

/** Starts ExaTENSOR numerical runtime. **/
int start(std::size_t hostMemBufferSize){
 int errc, hostArgMax, nGPU, listGPU[MAX_GPUS_PER_NODE];
 errc = talshGetDeviceCount(DEV_NVIDIA_GPU,&nGPU); if(errc != TALSH_SUCCESS) return -2;
 for(int i = 0; i < nGPU; ++i) listGPU[i]=i;
 errc = talshInit(&hostMemBufferSize,&hostArgMax,nGPU,listGPU,0,NULL,0,NULL);
 if(errc != TALSH_SUCCESS) return -1;
 return 0;
}

/** Stops ExaTENSOR numerical runtime. **/
int stop(){
 return talshShutdown();
}

} //end namespace exatensor
