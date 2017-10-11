/** C++ adapters for ExaTENSOR: Wrapper for importing dense tensor blocks from clients

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/10/11

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

#ifndef _EXA_TENSOR_DENSE_ADPT_H
#define _EXA_TENSOR_DENSE_ADPT_H

#include <memory>
#include <vector>
#include <assert.h>
#include <iostream>

#include "type_deduct.hpp"

#define _DEBUG_DIL

namespace exatensor {

/** Simple dense tensor wrapper with imported (external) body. **/
template <typename T>
class TensorDenseAdpt{

private:

 unsigned int Rank;                        //VAL: tensor rank (number of dimensions)
 std::unique_ptr<std::size_t[]> DimExtent; //VAL: tensor dimension extents
 std::shared_ptr<T> Body;                  //REF: pointer to the imported tensor body (tensor elements)

public:

//Life cycle:
 /** Constructs TensorDenseAdpt without a body (shape only). **/
 TensorDenseAdpt(const unsigned int rank, const std::size_t dimExtent[]);
 /** Constructs TensorDensAdpt with an externally provided body. **/
 TensorDenseAdpt(const unsigned int rank, const std::size_t dimExtent[], const std::shared_ptr<T> data);
 /** Copy constructor. **/
 TensorDenseAdpt(const TensorDenseAdpt & tensor);
 /** Copy assignment. **/
 TensorDenseAdpt & operator=(const TensorDenseAdpt & tensor);
 /** Destructor. **/
 virtual ~TensorDenseAdpt();

//Accessors:
 /** Returns the tensor rank. **/
 unsigned int getRank() const;
 /** Returns the extent of the specific tensor dimension. **/
 std::size_t getDimExtent(const unsigned int dimension) const;
 /** Returns a pointer to the tensor dimension extents. **/
 const std::size_t * getDimExtents() const;
 /** Returns a shared pointer to the tensor body (NULL inside if there is no body). **/
 std::shared_ptr<T> getBodyAccess() const;
 /** Returns the tensor volume (total number of tensor elements). **/
 std::size_t getVolume() const;
 /** Returns the tensor size in bytes. **/
 std::size_t getSize() const;
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
 /** Reshapes the tensor to a different shape. If the tensor has a body,
     it will be nullified until the new body is supplied.  **/
 void reshape(const unsigned int rank,        //in: new tensor rank
              const std::size_t dimExtent[]); //in: new tensor dimension extents

};

//Template definition:
#include "tensor_dense_adpt.cpp"

} //end namespace exatensor

#endif //_EXA_TENSOR_DENSE_ADPT_H
