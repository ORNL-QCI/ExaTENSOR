/** C++ adapters for ExaTENSOR: Header

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/07/21

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
 TensorDenseAdpt(unsigned int rank, const std::size_t dimExtent[]):
 Rank(rank), DimExtent(new std::size_t[rank]), Body(nullptr)
 {
  for(unsigned int i=0; i<rank; ++i) DimExtent[i]=dimExtent[i];
 }

 /** Constructs TensorDensAdpt with an externally provided body. **/
 TensorDenseAdpt(unsigned int rank, const std::size_t dimExtent[], std::shared_ptr<T> data):
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

 /** Returns a reference to a shared pointer to the tensor body, NULL if there is no body. **/
 const std::shared_ptr<T> & getBodyAccess() const {return Body;}

 /** Returns the tensor volume (total number of elements). **/
 std::size_t getVolume() const
 {
  std::size_t vol=1;
  for(unsigned int i=0; i<Rank; ++i) vol*=DimExtent[i];
  return vol;
 }

 /** Returns the tensor size in bytes. **/
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

 /** Reshapes the tensor to a different shape.
     The tensor is not allowed to have a body. **/
 void reshape(unsigned int rank,       //in: new tensor rank
              std::size_t dimExtent[]) //in: new tensor dimension extents
 {
#ifdef _DEBUG_DIL
  assert(!Body);
#endif
  DimExtent.reset(new std::size_t[rank]);
  for(unsigned int i=0; i<rank; ++i) DimExtent[i]=dimExtent[i];
  Rank=rank;
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

 //Mutators:

 /** Resets the tensor leg to another connection. **/
 void resetConnection(unsigned int tensorId, unsigned int dimesnId)
 {
  TensorId=tensorId; DimesnId=dimesnId;
  return;
 }

 /** Resets the tensor id in the tensor leg. **/
 void resetTensorId(unsigned int tensorId)
 {
  TensorId=tensorId;
  return;
 }

 /** Resets the tensor dimension id in the tensor leg. **/
 void resetDimensionId(unsigned int dimesnId)
 {
  DimesnId=dimesnId;
  return;
 }

};


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

 /** Returns the tensor rank. **/
 unsigned int getTensorRank() const {return Tensor.getRank();}

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
 unsigned int getNumLegs() const {return Legs.size();}

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

 //Mutators:

 /** Resets connection (leg). **/
 void resetConnection(unsigned int legId, const TensorLeg & tensorLeg)
 {
#ifdef _DEBUG_DIL
  assert(legId < Legs.size());
#endif
  Legs[legId].resetConnection(tensorLeg.getTensorId(),tensorLeg.getDimensionId());
  return;
 }

 /** Deletes the specified tensor dimension. **/
 void deleteDimension(unsigned int dimesn)
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
 void appendDimension(std::size_t dimExtent, const TensorLeg & leg)
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
  return (unsigned int)(Tensors.size()-1); //not counting the output tensor
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

 /** Explicitly appends a tensor to the tensor network, either input or
     output. The output (lhs) tensor must be appended first (tensor 0).
     Each next appended tensor will be considered an input (rhs) tensor. **/
 void appendTensor(const TensorDenseAdpt<T> & tensor,          //in: new tensor, either input (rhs) or output (lhs)
                   const std::vector<TensorLeg> & connections) //in: connections of the new tensor to other tensors via legs
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

 /** Appends a tensor to the tensor network by pairing some or all of its
     dimensions with the uncontracted dimensions of the tensor network.
     It is also fine to have none of the tensor legs be contracted with
     the tensor network, in which case they will simply be appended to
     the output tensor of the tensor network. **/
 void appendTensor(const TensorDenseAdpt<T> & tensor, //in: tensor being appended to the tensor network
                   const std::vector<std::pair<unsigned int, unsigned int>> & legPairs) //in: leg pairing: pair(tensor network leg id, tensor leg id)
 {
  auto tensRank = tensor.getRank(); //rank of the new tensor
  auto numLegs = legPairs.size(); //number of newly formed connections, can be zero
  auto nextTensorId = Tensors.size(); //id of the new tensor in the tensor network
#ifdef _DEBUG_DIL
  assert(numLegs <= tensRank);
#endif
  if(nextTensorId > 0){ //non-empty tensor network (at least the output and one input tensor)
#ifdef _DEBUG_DIL
   assert(nextTensorId >= 2);
#endif
   //get the current output tensor:
   auto & outTensor = Tensors[0]; //output (l.h.s.) tensor of the tensor network (TensorConn)
   auto outTensorRank = outTensor.getTensorRank(); //current rank of the output tensor
#ifdef _DEBUG_DIL
   assert(numLegs <= outTensorRank);
#endif
   std::vector<TensorLeg> legs(tensRank,TensorLeg(0,0)); //legs of the appended tensor
   //process all specified leg pairs:
   for(auto & legPair : legPairs){
    const auto ouLegId = legPair.first; //leg id in the tensor network
    const auto inLegId = legPair.second; //leg id in the appended tensor
    const auto & ouLeg = outTensor.getTensorLeg(ouLegId); //leg of the output tensor
    const auto tnTensorId = ouLeg.getTensorId(); //input tensor id (in the tensor network) that is connected to the output tensor
    const auto tnTensorLegId = ouLeg.getDimensionId(); //specific tensor leg of that input tensor
#ifdef _DEBUG_DIL
    assert(tnTensorId > 0);
#endif
    Tensors[tnTensorId].resetConnection(tnTensorLegId,TensorLeg(nextTensorId,inLegId)); //leg of the corresponding input tensor in the tensor network
    legs[inLegId].resetConnection(tnTensorId,tnTensorLegId); //connect the appended tensor to the tensor in the tensor network
    outTensor.resetConnection(ouLegId,TensorLeg(0,0)); //mark the disappeared dimension of the tensor network output tensor for subsequent deletion
   }
   //delete and shift other output tensor dimensions:
   if(numLegs > 0){
    unsigned int numDeleted=0;
    for(unsigned int i=0; i < outTensorRank; ++i){
     const auto j = i - numDeleted;
     const auto & ouLeg = outTensor.getTensorLeg(j); //leg of the output tensor
     const auto tnTensorId = ouLeg.getTensorId();
     if(tnTensorId == 0){ //delete
      outTensor.deleteDimension(j); //delete the disappeared dimension of the output tensor
      ++numDeleted;
     }else{ //shift
      const auto tnTensorLegId = ouLeg.getDimensionId();
      Tensors[tnTensorId].resetConnection(tnTensorLegId,TensorLeg(0,j));
     }
    }
#ifdef _DEBUG_DIL
    assert(numDeleted == numLegs);
#endif
   }
   //join the uncontracted dimensions of the appended tensor to the tensor network output tensor:
   outTensorRank = outTensor.getTensorRank(); //updated (decreased) rank of the output tensor
   for(unsigned int i = 0; i < tensRank; ++i){ //dimensions of the appended tensor
    auto & leg = legs[i];
    if(leg.getTensorId() == 0){ //uncontracted dimension
     leg.resetDimensionId(outTensorRank);
     outTensor.appendDimension(tensor.getDimExtent(i),TensorLeg(nextTensorId,i));
     outTensorRank = outTensor.getTensorRank();
    }
   }
   //append the input tensor to the tensor network:
   this->appendTensor(tensor,legs);
  }else{ //empty tensor network: two tensors will be appended (input and output)
#ifdef _DEBUG_DIL
   assert(numLegs == 0);
#endif
   //construct the output tensor without body (deferred):
   TensorDenseAdpt<T> outTensor(tensRank,tensor.getDimExtents());
   //construct the vector of legs for the output tensor:
   std::vector<TensorLeg> legs;
   for(unsigned int i = 0; i < tensRank; ++i) legs.push_back(TensorLeg(1,i));
   //append the output tensor to the tensor network:
   this->appendTensor(outTensor,legs);
   //construct the vector of legs for the 1st input tensor:
   for(auto & leg : legs) leg.resetTensorId(0);
   //append the input tensor to the tensor network:
   this->appendTensor(tensor,legs);
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
#ifdef _DEBUG_DIL
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
