/** C++ adapters for ExaTENSOR: Tensor leg (tensor connection to other tensors)

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

#ifndef _EXA_TENSOR_LEG_H
#define _EXA_TENSOR_LEG_H

#include <assert.h>
#include <iostream>

#include "type_deduct.hpp"

#define _DEBUG_DIL

namespace exatensor {

/** Tensor leg: Connection to another tensor **/
class TensorLeg{

private:

 unsigned int TensorId; //connected tensor id: 0 is the output tensor (lhs), >0 is an input tensor (rhs)
 unsigned int DimesnId; //connected tensor dimension: [0..rank-1], where "rank" is the rank of the connected tensor

public:

//Life cycle:
 /** Leg (connection) constructor. **/
 explicit TensorLeg(const unsigned int tensorId,  //connected tensor id in the tensor network
                    const unsigned int dimesnId); //connected tensor dimension

//Accesors:
 /** Returns the connected tensor id: [0..max] **/
 unsigned int getTensorId() const;
 /** Returns the connected tensor dimension: [0..rank-1] **/
 unsigned int getDimensionId() const;
 /** Print. **/
 void printIt() const;

//Mutators:
 /** Resets the tensor leg to another connection. **/
 void resetConnection(const unsigned int tensorId, const unsigned int dimesnId);
 /** Resets only the tensor id in the tensor leg. **/
 void resetTensorId(const unsigned int tensorId);
 /** Resets only the tensor dimension id in the tensor leg. **/
 void resetDimensionId(const unsigned int dimesnId);

};

} //end namespace exatensor

#endif //_EXA_TENSOR_LEG_H
