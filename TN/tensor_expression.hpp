/** C++ adapters for ExaTENSOR: Header

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/07/17

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
#include <vector>
#include <assert.h>
#include <iostream>

#define _DEBUG_DIL

namespace exatensor {

/** Simple dense tensor wrapper with imported body **/
template <typename T>
class TensorDenseAdpt{

private:

 unsigned int Rank;                        //VAL: tensor rank (number of dimensions)
 std::unique_ptr<std::size_t[]> DimExtent; //VAL: tensor dimension extents
 std::shared_ptr<T> Body;                  //REF: pointer to the imported tensor body (tensor elements)

public:

 //Life cycle:
 TensorDenseAdpt(unsigned int rank, std::size_t dimExtent[], std::shared_ptr<T> data):
 Rank(rank), DimExtent(new std::size_t[rank]), Body(data){
  for(unsigned int i=0; i<rank; ++i) DimExtent[i]=dimExtent[i];
 }

 TensorDenseAdpt(const TensorDenseAdpt & tensor):
 Rank(tensor.Rank), DimExtent(new std::size_t[tensor.Rank]), Body(tensor.Body){
  for(unsigned int i=0; i<tensor.Rank; ++i) DimExtent[i]=tensor.DimExtent[i];
 }

 TensorDenseAdpt & operator=(const TensorDenseAdpt & tensor){
  if(&tensor == this) return *this;
  if(tensor.Rank != Rank){
   DimExtent.reset(new std::size_t[tensor.Rank]);
   Rank=tensor.Rank;
  }
  std::copy(&tensor.DimExtent[0],&tensor.DimExtent[0]+tensor.Rank,&DimExtent[0]);
  Body=tensor.Body;
  return *this;
 }

 virtual ~TensorDenseAdpt(){}

 //Accessors:
 unsigned int getRank() const {return Rank;}

 std::size_t getDimExtent(unsigned int dimension) const{
#ifdef _DEBUG_DIL
  assert(dimension < Rank);
#endif
  return DimExtent[dimension];
 }

 const std::size_t * getDimExtents() const {return DimExtent.get();}

 std::shared_ptr<T> & getBodyAccess() const {return Body;}

 std::size_t getVolume() const{
  std::size_t vol=1;
  for(unsigned int i=0; i<Rank; ++i) vol*=DimExtent[i];
  return vol;
 }

 std::size_t getSize() const {return (this->getVolume())*sizeof(T);}

 //Print:
 void printIt() const{
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

};


/** Tensor leg: Connection to another tensor **/
class TensorLeg{

private:

 unsigned int TensorId; //connected tensor id: 0 is output tensor (lhs), >0 is input tensor (rhs)
 unsigned int DimesnId; //connected tensor dimension: [0..rank-1], where "rank" is the rank of the connected tensor

public:

 //Life cycle:
 TensorLeg():TensorId(0),DimesnId(0){}

 TensorLeg(unsigned int tensorId, unsigned int dimesnId):
 TensorId(tensorId), DimesnId(dimesnId){}

 //Accesors:
 unsigned int getTensorId() const {return TensorId;}
 unsigned int getDimensionId() const {return DimesnId;}

 //Print:
 void printIt() const{
  std::cout << "{" << TensorId << ":" << DimesnId << "}";
  return;
 }

};


/** Tensor connected to other tensors via tensor legs **/
template<typename T>
class TensorConn{

private:

 TensorDenseAdpt<T> Tensor;   //tensor
 std::vector<TensorLeg> Legs; //tensor legs (connections to other tensors): [1..rank]

public:

 //Life cycle:
 TensorConn(const TensorDenseAdpt<T> & tensor, const std::vector<TensorLeg> & connections):
 Tensor(tensor), Legs(connections){
#ifdef _DEBUG_DIL
  assert(tensor.getRank() == connections.size());
#endif
 }

 virtual ~TensorConn(){}

 //Accessors:
 std::size_t getDimExtent(unsigned int dimension) const{
  return Tensor.getDimExtent(dimension);
 }

 const TensorLeg & getTensorLeg(unsigned int leg) const{
#ifdef _DEBUG_DIL
  assert(leg < Tensor.getRank());
#endif
  return Legs.at(leg);
 }

 unsigned int getNumLegs() const {return Tensor.getRank();}

 //Print:
 void printIt() const{
  //std::cout << std::endl;
  std::cout << "TensorConn{" << std::endl;
  Tensor.printIt();
  std::cout << "Legs:";
  for(unsigned int i=0; i<Tensor.getRank(); ++i){std::cout << " "; Legs.at(i).printIt();}
  std::cout << std::endl << "}" << std::endl;
  return;
 }

};


/** Tensor network (contraction of multiple tensors):
 A tensor network consists of tensors numerated from 0.
 Tensor 0 is always the output (lhs) tensor consisting of
 uncontracted legs. Tensors [1..max] are input (rhs) tensors.
 Legs of the input tensors that are left uncontracted define
 the output tensor (tensor 0) by definition. **/
template<typename T>
class TensorNetwork{

private:

 std::vector<TensorConn<T>> Tensors; //tensors: [0;1..num_rhs_tensors]

public:

 //Life cycle:

 TensorNetwork(){}

 virtual ~TensorNetwork(){}

 //Accessors:

 /** Returns the number of r.h.s. tensors in the tensor network.
     Note that the output (l.h.s.) tensor 0 is not counted here. **/
 unsigned int getNumTensors() const
 {
  return (unsigned int)(Tensors.size()-1);
 }

 //Mutators:

 /** Appends a tensor to the tensor network **/
 void appendTensor(const TensorDenseAdpt<T> & tensor,          //in: new tensor
                   const std::vector<TensorLeg> & connections) //in: connections of the new tensor to other tensors
 {
  auto num_tens = Tensors.size(); //current total number of tensors in the tensor network
  //Check the consistency of the new tensor candidate:
#ifdef _DEBUG_DIL
  assert(tensor.getRank() == connections.size());
  unsigned int i=0;
  for(auto it=connections.cbegin(); it != connections.cend(); ++it){
   const TensorLeg & leg = *it; //new tensor leg
   auto tens_id = leg.getTensorId(); //tensor to which the new leg is connected
   if(tens_id < num_tens){ //that tensor has already been appended into the tensor network
    TensorConn<T> & tensconn = Tensors[tens_id]; //reference to that tensor
    auto dimsn = leg.getDimensionId(); //specific dimension of that tensor
    const TensorLeg & other_leg = tensconn.getTensorLeg(dimsn); //leg on the other side
    assert(other_leg.getTensorId() == num_tens && other_leg.getDimensionId() == i); //legs connectivity must match
    assert(tensor.getDimExtent(i) == tensconn.getDimExtent(dimsn)); //dimension extents must match as well
   }else if(tens_id == num_tens){ //self-contraction
    auto dimsn = leg.getDimensionId(); //specific dimension of the same tensor
    assert(dimsn != i); //dimension of a tensor cannot be contracted with itself
    const TensorLeg & other_leg = connections.at(dimsn); //other leg of the same tensor (loop)
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

 //Transforms:

 /** Contracts two tensors in a given tensor network and returns the result as a new tensor network. **/
 TensorNetwork<T> * contractTensors(unsigned int tens_id1, //in: id of the 1st tensor: [1..max]
                                    unsigned int tens_id2) //in: id of the 2nd tensor: [1..max]
 {
  if(tens_id1 > tens_id2) std::swap(tens_id1,tens_id2);
  auto tensornet = new TensorNetwork<T>();

  return tensornet;
 }

 /** Appends another tensor network into the current tensor network. **/
 void appendNetwork(const TensorNetwork<T> & tensornet, //in: tensor network
                    const unsigned int dims0[],         //in: list of output dimensions of the curent tensor network
                    const unsigned int dims1[])         //in: list of matching output dimensions of the appended tensor network
 {

  return;
 }

 //Print:
 void printIt() const
 {
  std::cout << "TensorNetwork{" << std::endl;
  unsigned int NumInputTensors = this->getNumTensors();
  std::cout << "Number of input tensors = " << NumInputTensors << std::endl;
  for(unsigned int i = 0; i <= NumInputTensors; ++i) Tensors[i].printIt();
  std::cout << "}" << std::endl;
  return;
 }

};

} //end namespace exatensor

#endif //_TENSOR_EXPRESSION_H
