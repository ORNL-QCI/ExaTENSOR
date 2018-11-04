/** C++ adapters for ExaTENSOR: Tensor network solver

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

#ifndef EXA_TENSOR_SOLVER_H_
#define EXA_TENSOR_SOLVER_H_

#include <vector>

#include "tensor_network.hpp"

namespace exatensor {

/** Optimizes a given subset of tensors in a closed tensor network
    in order to maximize its scalar value. If tensorNorms is not empty,
    its length must match that of optimizedTensIds. In this case
    the optimized tensors will be constrained to the given norms.
    If tensorNorms is empty, the input norms of the optimized
    tensors will be kept intact (they must not be zero). **/
template<typename T>
int optimizeOverlapMax(TensorNetwork<T> tensNet,                         //inout: closed tensor network
                       const std::vector<unsigned int> optimizedTensIds, //in: IDs of the r.h.s. tensors to be optimized
                       const std::vector<double> tensorNorms);           //in: imposed tensor norms

//Template definition:
#include "tensor_solver.cpp"

} //end namespace exatensor

#endif //EXA_TENSOR_SOLVER_H_
