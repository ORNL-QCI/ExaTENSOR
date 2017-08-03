/** C++ adapters for ExaTENSOR: Tensor leg (connection)

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/08/03

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

#include "tensor_leg.hpp"

namespace exatensor {

//Life cycle:

/** Leg (connection) constructor. **/
TensorLeg::TensorLeg(const unsigned int tensorId,  //connected tensor id in the tensor network
                     const unsigned int dimesnId): //connected tensor dimension
 TensorId(tensorId), DimesnId(dimesnId)
{
}

//Accesors:

/** Returns the connected tensor id: [0..max] **/
unsigned int TensorLeg::getTensorId() const
{
 return TensorId;
}

/** Returns the connected tensor dimension: [0..rank-1] **/
unsigned int TensorLeg::getDimensionId() const
{
 return DimesnId;
}

/** Print. **/
void TensorLeg::printIt() const
{
 std::cout << "{" << TensorId << ":" << DimesnId << "}";
 return;
}

//Mutators:

/** Resets the tensor leg to another connection. **/
void TensorLeg::resetConnection(const unsigned int tensorId, const unsigned int dimesnId)
{
 TensorId=tensorId; DimesnId=dimesnId;
 return;
}

/** Resets the tensor id in the tensor leg. **/
void TensorLeg::resetTensorId(const unsigned int tensorId)
{
 TensorId=tensorId;
 return;
}

/** Resets the tensor dimension id in the tensor leg. **/
void TensorLeg::resetDimensionId(const unsigned int dimesnId)
{
 DimesnId=dimesnId;
 return;
}

} //end namespace exatensor
