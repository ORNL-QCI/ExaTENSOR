/** C++ adapters for ExaTENSOR: Tensor connected to other tensors

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

#ifndef _EXA_TENSOR_CONN_H
#define _EXA_TENSOR_CONN_H

#include <assert.h>
#include <iostream>
#include <memory>
#include <vector>

#include "type_deduct.hpp"

#include "tensor_leg.hpp"
#include "tensor_dense_adpt.hpp"

#define _DEBUG_DIL

namespace exatensor {

/** Tensor connected to other tensors via tensor legs **/
template<typename T>
class TensorConn{

private:

 TensorDenseAdpt<T> Tensor;   //tensor
 std::vector<TensorLeg> Legs; //tensor legs (connections to other tensors)

public:

//Life cycle:
 /** Constructor of a connected tensor. **/
 TensorConn(const TensorDenseAdpt<T> & tensor,           //tensor
            const std::vector<TensorLeg> & connections); //tensor connections (legs) to other tensors in a tensor network
 /** Destructor. **/
 virtual ~TensorConn();

//Accessors:
 /** Returns a const reference to the tensor. **/
 const TensorDenseAdpt<T> & getTensor() const;
 /** Returns the tensor rank. **/
 unsigned int getTensorRank() const;
 /** Returns the tensor volume (number of tensor elements). **/
 std::size_t getVolume() const;
 /** Returns the extent of a specific tensor dimension. **/
 std::size_t getDimExtent(const unsigned int dimension) const;
 /** Returns a const-reference to the specific tensor leg (connection to other tensors). **/
 const TensorLeg & getTensorLeg(const unsigned int leg) const;
 /** Returns the total number of tensor legs (connections). **/
 unsigned int getNumLegs() const;
 /** Returns true if the tensor has body. **/
 bool hasBody() const;
 /** Prints. **/
 void printIt() const;

//Mutators:
 /** Associates the tensor with an externally provided tensor body.
     Will fail if the tensor body is already present (defined).
     The new body may be null. **/
 void setBody(const std::shared_ptr<T> body);
 /** Reassociates the tensor with another body. The new body may be null. **/
 void resetBody(const std::shared_ptr<T> body);
 /** Allocates tensor body. **/
 void allocateBody();
 /** Sets tensor body to zero. **/
 void nullifyBody();
 /** Resets connection (leg). **/
 void resetConnection(const unsigned int legId, const TensorLeg & tensorLeg);
 /** Deletes the specified tensor dimension. **/
 void deleteDimension(const unsigned int dimesn);
 /** Appends a new dimension to the connected tensor as the last dimension. **/
 void appendDimension(const std::size_t dimExtent, //in: new dimension extent
                      const TensorLeg & leg);      //in: new tensor leg for the new dimension

};

//Template definition:
#include "tensor_conn.cpp"

} //end namespace exatensor

#endif //_EXA_TENSOR_CONN_H
