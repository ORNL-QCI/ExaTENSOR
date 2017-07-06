/** C++ adapters for ExaTENSOR: Header

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

#ifndef _TENSOR_EXPRESSION_H
#define _TENSOR_EXPRESSION_H

#include <memory>

namespace exatensor {

template <typename T>
class TensorAdpt{

private:

 unsigned int Rank;                        //VAL: tensor rank (number of dimensions)
 std::unique_ptr<std::size_t[]> DimExtent; //VAL: tensor dimension extents
 std::shared_ptr<T> Body;                  //REF: pointer to the tensor body (tensor elements)

public:

 //Life cycle:
 TensorAdpt(unsigned int rank, std::size_t dimExtent[], std::shared_ptr<T> data):
 Rank(rank), DimExtent(new std::size_t[rank]), Body(data){
  for(unsigned int i=0; i<rank; ++i) DimExtent[i]=dimExtent[i];
 }

 TensorAdpt(const TensorAdpt & tensor):
 Rank(tensor.Rank), DimExtent(new std::size_t[tensor.rank]), Body(tensor.Body){
  for(unsigned int i=0; i<tensor.rank; ++i) DimExtent[i]=tensor.DimExtent[i];
 }

 TensorAdpt & operator=(const TensorAdpt & tensor){
  if(&tensor == this) return *this;
  if(tensor.Rank != Rank){
   DimExtent.reset(new std::size_t[tensor.Rank]);
   Rank=tensor.Rank;
  }
  std::copy(&tensor.DimExtent[0],&tensor.DimExtent[0]+Rank,&DimExtent[0]);
  Body=tensor.Body;
  return *this;
 }

 virtual ~TensorAdpt(){}

 //Accessors:
 unsigned int getRank(){return Rank;}

 const std::size_t * getDimExtents(){return DimExtent.get();}

 std::shared_ptr<T> & getBodyAccess(){return Body;}

 std::size_t getVolume(){
  std::size_t vol=1;
  for(unsigned int i=0; i<Rank; ++i) vol*=DimExtent[i];
  return vol;
 }

 std::size_t getSize(){return (this->getVolume())*sizeof(T);}

};

} //end namespace exatensor

#endif //_TENSOR_EXPRESSION_H
