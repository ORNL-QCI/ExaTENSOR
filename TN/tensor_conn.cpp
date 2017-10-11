/** C++ adapters for ExaTENSOR: Tensor connected to other tensors

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

//Life cycle:

/** Constructor of a connected tensor. **/
template <typename T>
TensorConn<T>::TensorConn(const TensorDenseAdpt<T> & tensor,           //tensor
                          const std::vector<TensorLeg> & connections): //tensor connections (legs) to other tensors in a tensor network
 Tensor(tensor), Legs(connections)
{
#ifdef _DEBUG_DIL
 assert(tensor.getRank() == connections.size());
#endif
}

/** Destructor. **/
template <typename T>
TensorConn<T>::~TensorConn()
{
}

//Accessors:

/** Returns a const reference to the tensor. **/
template <typename T>
const TensorDenseAdpt<T> & TensorConn<T>::getTensor() const
{
 return Tensor;
}

/** Returns the tensor rank. **/
template <typename T>
unsigned int TensorConn<T>::getTensorRank() const
{
 return Tensor.getRank();
}

/** Returns the tensor volume (number of tensor elements). **/
template <typename T>
std::size_t TensorConn<T>::getVolume() const
{
 return Tensor.getVolume();
}

/** Returns the extent of a specific tensor dimension. **/
template <typename T>
std::size_t TensorConn<T>::getDimExtent(const unsigned int dimension) const
{
 return Tensor.getDimExtent(dimension);
}

/** Returns a const-reference to the specific tensor leg (connection to other tensors). **/
template <typename T>
const TensorLeg & TensorConn<T>::getTensorLeg(const unsigned int leg) const
{
#ifdef _DEBUG_DIL
 assert(leg < Tensor.getRank());
#endif
 return Legs.at(leg);
}

/** Returns the total number of tensor legs (connections). **/
template <typename T>
unsigned int TensorConn<T>::getNumLegs() const
{
 return Legs.size();
}

/** Returns true if the tensor has body. **/
template <typename T>
bool TensorConn<T>::hasBody() const
{
 return Tensor.hasBody();
}

/** Prints. **/
template <typename T>
void TensorConn<T>::printIt() const
{
 //std::cout << std::endl;
 std::cout << "TensorConn{" << std::endl;
 Tensor.printIt();
 std::cout << "Legs:";
 for(unsigned int i=0; i<Tensor.getRank(); ++i){std::cout << " "; Legs.at(i).printIt();}
 std::cout << std::endl << "}" << std::endl;
 return;
}

//Mutators:

/** Associates the tensor with an externally provided tensor body.
    Will fail if the tensor body is already present (defined).
    The new body may be null. **/
template <typename T>
void TensorConn<T>::setBody(const std::shared_ptr<T> body)
{
 Tensor.setBody(body);
 return;
}

/** Reassociates the tensor with another body. The new body may be null. **/
template <typename T>
void TensorConn<T>::resetBody(const std::shared_ptr<T> body)
{
 Tensor.resetBody(body);
 return;
}

/** Allocates tensor body. **/
template <typename T>
void TensorConn<T>::allocateBody()
{
 Tensor.allocateBody();
 return;
}

/** Sets tensor body to zero. **/
template <typename T>
void TensorConn<T>::nullifyBody()
{
 Tensor.nullifyBody();
 return;
}

/** Resets connection (leg). **/
template <typename T>
void TensorConn<T>::resetConnection(const unsigned int legId, const TensorLeg & tensorLeg)
{
#ifdef _DEBUG_DIL
 assert(legId < Legs.size());
#endif
 Legs[legId].resetConnection(tensorLeg.getTensorId(),tensorLeg.getDimensionId());
 return;
}

/** Deletes the specified tensor dimension. **/
template <typename T>
void TensorConn<T>::deleteDimension(const unsigned int dimesn)
{
 auto oldTensRank = Tensor.getRank();
#ifdef _DEBUG_DIL
 assert(dimesn < oldTensRank);
#endif
 auto newTensRank = oldTensRank - 1;
 const std::size_t * oldDims = Tensor.getDimExtents();
 std:size_t newDims[newTensRank];
 unsigned int j=0;
 for(unsigned int i = 0; i < oldTensRank; ++i){
  if(i != dimesn) newDims[j++] = oldDims[i];
 }
 Tensor.reshape(newTensRank,newDims);
 Legs.erase(Legs.begin()+dimesn);
 return;
}

/** Appends a new dimension to the connected tensor as the last dimension. **/
template <typename T>
void TensorConn<T>::appendDimension(const std::size_t dimExtent, //in: new dimension extent
                                    const TensorLeg & leg)       //in: new tensor leg for the new dimension
{
 auto oldTensRank = Tensor.getRank();
 auto newTensRank = oldTensRank + 1;
 const std::size_t * oldDims = Tensor.getDimExtents();
 std:size_t newDims[newTensRank];
 for(unsigned int i = 0; i < oldTensRank; ++i) newDims[i] = oldDims[i];
 newDims[newTensRank-1]=dimExtent;
 Tensor.reshape(newTensRank,newDims);
 Legs.push_back(leg);
 return;
}
