/** C++ adapters for ExaTENSOR: Header

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/07/18

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

/** Simple dense tensor wrapper with imported body. **/
template <typename T>
class TensorDenseAdpt{

private:

 unsigned int Rank;                        //VAL: tensor rank (number of dimensions)
 std::unique_ptr<std::size_t[]> DimExtent; //VAL: tensor dimension extents
 std::shared_ptr<T> Body;                  //REF: pointer to the imported tensor body (tensor elements)

public:

 //Life cycle:

 /** Constructs TensorDenseAdpt without a body (shape only). **/
 TensorDenseAdpt(unsigned int rank, std::size_t dimExtent[]):
 Rank(rank), DimExtent(new std::size_t[rank]), Body(nullptr)
 {
  for(unsigned int i=0; i<rank; ++i) DimExtent[i]=dimExtent[i];
 }

 /** Constructs TensorDensAdpt with an externally provided body. **/
 TensorDenseAdpt(unsigned int rank, std::size_t dimExtent[], std::shared_ptr<T> data):
 Rank(rank), DimExtent(new std::size_t[rank]), Body(data)
 {
  for(unsigned int i=0; i<rank; ++i) DimExtent[i]=dimExtent[i];
 }

 /** Copy constructor. **/
 TensorDenseAdpt(const TensorDenseAdpt & tensor):
 Rank(tensor.Rank), DimExtent(new std::size_t[tensor.Rank]), Body(tensor.Body)
 {
  for(unsigned int i=0; i<tensor.Rank; ++i) DimExtent[i]=tensor.DimExtent[i];
 }

 /** Copy assignment. **/
 TensorDenseAdpt & operator=(const TensorDenseAdpt & tensor)
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
 virtual ~TensorDenseAdpt(){}

 //Accessors:

 /** Returns tensor rank. **/
 unsigned int getRank() const {return Rank;}

 /** Returns the extent of the specific tensor dimension. **/
 std::size_t getDimExtent(unsigned int dimension) const
 {
#ifdef _DEBUG_DIL
  assert(dimension < Rank);
#endif
  return DimExtent[dimension];
 }

 /** Returns a pointer to the tensor dimension extents. **/
 const std::size_t * getDimExtents() const {return DimExtent.get();}

 /** Returns a shared pointer to the tensor body, NULL if there is no body. **/
 std::shared_ptr<T> & getBodyAccess() const {return Body;}

 /** Returns tensor volume (total number of elements). **/
 std::size_t getVolume() const
 {
  std::size_t vol=1;
  for(unsigned int i=0; i<Rank; ++i) vol*=DimExtent[i];
  return vol;
 }

 /** Returns tensor size in bytes. **/
 std::size_t getSize() const {return (this->getVolume())*sizeof(T);}

 /** Prints. **/
 void printIt() const
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
     Will fail if the tensor body is already present. **/
 void setBody(std::shared_ptr<T> body)
 {
  assert(!Body);
  Body=body;
  return;
 }

 /** Reassociates the tensor with another body. **/
 void resetBody(std::shared_ptr<T> body)
 {
  if(Body) Body.reset();
  Body=body;
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

 /** Default constructor (empty connection). **/
 TensorLeg():TensorId(0),DimesnId(0){}

 /** Leg (connection) constructor. **/
 TensorLeg(unsigned int tensorId,  //connected tensor id in the tensor network
           unsigned int dimesnId): //connected tensor dimension
 TensorId(tensorId), DimesnId(dimesnId){}

 //Accesors:

 /** Returns the connected tensor id: [0..max] **/
 unsigned int getTensorId() const {return TensorId;}

 /** Returns the connected tensor dimension: [0..rank-1] **/
 unsigned int getDimensionId() const {return DimesnId;}

 /** Print. **/
 void printIt() const
 {
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

 /** Constructor of a connected tensor. **/
 TensorConn(const TensorDenseAdpt<T> & tensor,           //tensor
            const std::vector<TensorLeg> & connections): //tensor connections (legs) to other tensors in a tensor network
 Tensor(tensor), Legs(connections)
 {
#ifdef _DEBUG_DIL
  assert(tensor.getRank() == connections.size());
#endif
 }

 /** Destructor. **/
 virtual ~TensorConn(){}

 //Accessors:

 /** Returns the extent of a specific tensor dimension. **/
 std::size_t getDimExtent(unsigned int dimension) const
 {
  return Tensor.getDimExtent(dimension);
 }

 /** Returns a specific tensor leg (connection to other tensors). **/
 const TensorLeg & getTensorLeg(unsigned int leg) const
 {
#ifdef _DEBUG_DIL
  assert(leg < Tensor.getRank());
#endif
  return Legs.at(leg);
 }

 /** Returns the total number of tensor legs (connections). **/
 unsigned int getNumLegs() const {return Tensor.getRank();}

 /** Prints. **/
 void printIt() const
 {
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

 /** Constructs an empty tensor network **/
 TensorNetwork(){}

 /** Destructor. **/
 virtual ~TensorNetwork(){}

 //Accessors:

 /** Returns the number of r.h.s. tensors in the tensor network.
     Note that the output (l.h.s.) tensor 0 is not counted here. **/
 unsigned int getNumTensors() const
 {
  return (unsigned int)(Tensors.size()-1);
 }

 /** Prins. **/
 void printIt() const
 {
  std::cout << "TensorNetwork{" << std::endl;
  unsigned int NumInputTensors = this->getNumTensors();
  std::cout << "Number of input tensors = " << NumInputTensors << std::endl;
  for(unsigned int i = 0; i <= NumInputTensors; ++i) Tensors[i].printIt();
  std::cout << "}" << std::endl;
  return;
 }

 //Mutators:

 /** Appends a tensor to the tensor network, either input or output.
     The output (lhs) tensor must be appended first (tensor 0). Each
     next appended tensor will be considered an input (rhs) tensor. **/
 void appendTensor(const TensorDenseAdpt<T> & tensor,          //in: new tensor, either input or output
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

 /** Appends a tensor to the tensor network by matching some or all of its
     dimensions, specified in dimsTens[], with the uncontracted dimensions
     of the tensor network. The input argument dimsNetw[] allows selecting
     specific uncontracted dimensions of the tensor network to connect to. **/
 void appendTensor(const TensorDenseAdpt<T> & tensor,       //in: tensor being appended
                   const unsigned int numConnections,       //in: number of connections between the tensor and tensor network
                   const unsigned int dimsTens[] = nullptr, //in: list of connected dimensions of the appended tensor
                   const unsigned int dimsNetw[] = nullptr) //in: list of matching output dimensions of the tensor network
 {
  auto tensRank = tensor.getRank();
  if(Tensors.size() > 0){ //non-empty tensor network
   auto & outTensor = Tensors[0]; //output (l.h.s.) tensor of the tensor network
   
  }else{ //empty tensor network
   
  }
  return;
 }

 /** Appends another tensor network into the current tensor network. **/
 void appendNetwork(const TensorNetwork<T> & tensornet,     //in: another tensor network
                    const unsigned int numConnections,      //in: number of connections between tensor networks
                    const unsigned int dimsNew[] = nullptr, //in: list of the output dimensions of the appended tensor network to connect with
                    const unsigned int dimsOld[] = nullptr) //in: list of the matching output dimensions of the current tensor network
 {

  return;
 }

 /** Associates the output (lhs) tensor with its externally provided body. **/
 void setOutputBody(std::shared_ptr<T> body)
 {
#ifdef _DIL_DEBUG
  assert(Tensors.size() > 0 && body);
#endif
  Tensors[0].setBody(body);
  return;
 }

 //Transforms:

 /** Contracts two tensors in a given tensor network and returns the result as a new tensor network. **/
 TensorNetwork<T> * contractTensors(unsigned int tens_id1, //in: id of the 1st tensor in the tensor network: [1..max]
                                    unsigned int tens_id2) //in: id of the 2nd tensor in the tensor network: [1..max]
 {
#ifdef _DEBUG_DIL
  assert((tens_id1 > 0 && tens_id1 <= this->getNumTensors()) && (tens_id2 > 0 && tens_id2 <= this->getNumTensors()));
#endif
  if(tens_id1 > tens_id2){unsigned int tmp = tens_id1; tens_id1 = tens_id2; tens_id2 = tmp;}
  auto tensornet = new TensorNetwork<T>();
  auto output = (this->Tensors)[0];
  for(unsigned int i = 0; i < tens_id1; ++i){
   
  }
  return tensornet;
 }

};

} //end namespace exatensor

#endif //_TENSOR_EXPRESSION_H
