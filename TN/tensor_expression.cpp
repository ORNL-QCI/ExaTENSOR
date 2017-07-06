/** C++ adapters for ExaTENSOR: Implementation

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/07/06

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

#include "tensor_expression.hpp"

namespace exatensor{

template<typename T>
void TensorNetwork<T>::appendTensor(const TensorDenseAdpt<T> & tensor, const std::vector<TensorLeg> & connections){

 auto num_tens = Tensors.size(); //current total number of tensors set in the tensor network
 //Check the consistency of the new tensor candidate:
#ifdef _DEBUG_DIL
 assert(num_tens < (1 + NumInputTensors));
 assert(tensor.getRank() == connections.size());
 unsigned int i=0;
 for(auto it=connections.cbegin(); it != connections.cend(); ++it){
  TensorLeg & leg = *it; //new tensor leg
  auto tens_id = leg.getTensorId(); //tensor to which the new leg is connected
  assert(tens_id <= NumInputTensors);
  if(tens_id < num_tens){ //that tensor has already been appended into the tensor network
   TensorConn<T> & tensconn = Tensors[tens_id]; //reference to that tensor
   auto dimsn = leg.getDimensionId(); //specific dimension of that tensor
   TensorLeg & other_leg = tensconn.getTensorLeg(dimsn); //leg on the other side
   assert(other_leg.getTensorId() == num_tens && other_leg.getDimensionId() == i); //legs connectivity must match
   assert(tensor.getDimExtent(i) == tensconn.getDimExtent(dimsn)); //dimension extents must match as well
  }else if(tens_id == num_tens){ //self-contraction
   auto dimsn = leg.getDimensionId(); //specific dimension of the same tensor
   assert(dimsn != i); //dimension of a tensor cannot be contracted with itself
   TensorLeg & other_leg = connections.at(dimsn); //other leg of the same tensor (loop)
   assert(other_leg.getTensorId() == num_tens && other_leg.getDimensionId() == i); //legs connectivity must match
   assert(tensor.getDimExtent(i) == tensor.getDimExtent(dimsn)); //dimension extents must match as well
  }
  ++i;
 }
#endif
 //append the new tensor into the tensor network:
 Tensors.push_back(TensorConn<T>(tensor,connections));
 return;
}

} //end namespace exatensor
