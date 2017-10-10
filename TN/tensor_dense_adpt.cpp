/** C++ adapters for ExaTENSOR: Wrapper for importing dense tensor blocks from clients

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/10/10

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

/** Constructs TensorDenseAdpt without a body (shape only). **/
template <typename T>
TensorDenseAdpt<T>::TensorDenseAdpt(const unsigned int rank, const std::size_t dimExtent[]):
 Rank(rank), DimExtent(new std::size_t[rank]), Body(nullptr)
{
 for(unsigned int i=0; i<rank; ++i) DimExtent[i]=dimExtent[i];
}

/** Constructs TensorDensAdpt with an externally provided body. **/
template <typename T>
TensorDenseAdpt<T>::TensorDenseAdpt(const unsigned int rank, const std::size_t dimExtent[], const std::shared_ptr<T> data):
 Rank(rank), DimExtent(new std::size_t[rank]), Body(data)
{
 for(unsigned int i=0; i<rank; ++i) DimExtent[i]=dimExtent[i];
}

/** Copy constructor. **/
template <typename T>
TensorDenseAdpt<T>::TensorDenseAdpt(const TensorDenseAdpt<T> & tensor):
 Rank(tensor.Rank), DimExtent(new std::size_t[tensor.Rank]), Body(tensor.Body)
{
 for(unsigned int i=0; i<tensor.Rank; ++i) DimExtent[i]=tensor.DimExtent[i];
}

/** Copy assignment. **/
template <typename T>
TensorDenseAdpt<T> & TensorDenseAdpt<T>::operator=(const TensorDenseAdpt<T> & tensor)
{
 if(&tensor == this) return *this;
 if(tensor.Rank != Rank){
  DimExtent.reset(new std::size_t[tensor.Rank]);
  Rank=tensor.Rank;
 }
 std::copy(&tensor.DimExtent[0],&tensor.DimExtent[0]+tensor.Rank,&DimExtent[0]);
 Body=tensor.Body;
 return *this;
}

/** Destructor. **/
template <typename T>
TensorDenseAdpt<T>::~TensorDenseAdpt()
{
 //std::cout << std::endl << "exatensor::~TensorDenseAdpt: Body reference count = " << Body.use_count(); //debug
}

//Accessors:

/** Returns the tensor rank. **/
template <typename T>
unsigned int TensorDenseAdpt<T>::getRank() const
{
 return Rank;
}

/** Returns the extent of the specific tensor dimension. **/
template <typename T>
std::size_t TensorDenseAdpt<T>::getDimExtent(const unsigned int dimension) const
{
#ifdef _DEBUG_DIL
 assert(dimension < Rank);
#endif
 return DimExtent[dimension];
}

/** Returns a pointer to the tensor dimension extents. **/
template <typename T>
const std::size_t * TensorDenseAdpt<T>::getDimExtents() const
{
 return DimExtent.get();
}

/** Returns a shared pointer to the tensor body (NULL inside if there is no body). **/
template <typename T>
std::shared_ptr<T> TensorDenseAdpt<T>::getBodyAccess() const
{
 return Body;
}

/** Returns the tensor volume (total number of tensor elements). **/
template <typename T>
std::size_t TensorDenseAdpt<T>::getVolume() const
{
 std::size_t vol=1;
 for(unsigned int i=0; i<Rank; ++i) vol*=DimExtent[i];
 return vol;
}

/** Returns the tensor size in bytes. **/
template <typename T>
std::size_t TensorDenseAdpt<T>::getSize() const
{
 return (this->getVolume())*sizeof(T);
}

/** Returns true if the tensor has body. **/
template <typename T>
bool TensorDenseAdpt<T>::hasBody() const
{
 return (Body != nullptr);
}

/** Prints. **/
template <typename T>
void TensorDenseAdpt<T>::printIt() const
{
 //std::cout << std::endl;
 std::cout << "TensorDenseAdpt{" << std::endl;
 std::cout << " Rank = " << Rank << std::endl;
 std::cout << " Dim extents:";
 for(unsigned int i=0; i<Rank; ++i) std::cout << " " << DimExtent[i];
 std::cout << std::endl;
 std::cout << " Data pointer: " << Body.get() << std::endl;
 std::cout << "}" << std::endl;
 return;
}

//Mutators:

/** Associates the tensor with an externally provided tensor body.
    Will fail if the tensor body is already present (defined).
    The new body may be null. **/
template <typename T>
void TensorDenseAdpt<T>::setBody(const std::shared_ptr<T> body)
{
 assert(!Body);
 Body=body;
 return;
}

/** Reassociates the tensor with another body. The new body may be null. **/
template <typename T>
void TensorDenseAdpt<T>::resetBody(const std::shared_ptr<T> body)
{
 if(Body) Body.reset();
 Body=body;
 return;
}

/** Allocates tensor body. **/
template <typename T>
void TensorDenseAdpt<T>::allocateBody()
{
 assert(!Body);
 auto vol = this->getVolume();
 assert(vol > 0);
 Body.reset(new T[vol], [](T * ptr){delete[] ptr;});
 assert(Body);
 return;
}

/** Reshapes the tensor to a different shape. If the tensor has a body,
    it will be nullified until the new body is supplied.  **/
template <typename T>
void TensorDenseAdpt<T>::reshape(const unsigned int rank,       //in: new tensor rank
                                 const std::size_t dimExtent[]) //in: new tensor dimension extents
{
 Body.reset(); //the old body to be gone
 DimExtent.reset(new std::size_t[rank]);
 for(unsigned int i=0; i<rank; ++i) DimExtent[i]=dimExtent[i];
 Rank=rank;
 return;
}
