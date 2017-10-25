/** C++ adapters for ExaTENSOR: Tensor network

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/10/25

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

#ifndef _EXA_TENSOR_NETWORK_H
#define _EXA_TENSOR_NETWORK_H

#include <memory>
#include <tuple>
#include <vector>
#include <queue>
#include <map>
#include <string>
#include <cstdint>
#include <assert.h>
#include <iostream>
#include <ctime>
#include <chrono>

#include <complex>

#include "type_deduct.hpp"

#include "tensor_conn.hpp"

#include "tensor_define.hpp"

#include "talsh.h"

#define _DEBUG_DIL

namespace exatensor {

//Types:
using ContractionSequence = std::vector<std::pair<unsigned int, unsigned int>>;

//Traits:
template <typename DTK>
struct TensorDataKind{
 static const int Type = 0;
};

template <>
struct TensorDataKind<float>{
 static const int Type = R4;
};

template <>
struct TensorDataKind<double>{
 static const int Type = R8;
};

template <>
struct TensorDataKind<std::complex<float>>{
 static const int Type = C4;
};

template <>
struct TensorDataKind<std::complex<double>>{
 static const int Type = C8;
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

//Data members:
 std::vector<TensorConn<T>> Tensors; //interconnected tensors: [0;1..num_rhs_tensors]

//Constants:
 static const unsigned int NumWalkersDefault = 1024; //number of walkers for tensor contraction sequence optimization

public:

//Life cycle:
 /** Constructs an empty tensor network **/
 TensorNetwork();
 /** Copy constructor. **/
 TensorNetwork(const TensorNetwork<T> & tensNetwork);
 /** Destructor. **/
 virtual ~TensorNetwork();

//Accessors:
 /** Returns TRUE if the tensor network is empty, FALSE otherwise. **/
 bool isEmpty() const;
 /** Returns the number of r.h.s. tensors in the tensor network.
     Note that the output (l.h.s.) tensor 0 is not counted here. **/
 unsigned int getNumTensors() const;
 /** Returns a const reference to a specific tensor from the tensor network. **/
 const TensorDenseAdpt<T> & getTensor(const unsigned int id) const;
 /** Returns a const reference to a specific tensor from the tensor network together with its connections. **/
 const TensorConn<T> & getTensorConn(const unsigned int id) const;
 /** Prints. **/
 void printIt() const;

//Mutators:
 /** Explicitly appends a tensor to the tensor network, either input or
     output. The output (lhs) tensor must be appended first (tensor 0).
     Each next appended tensor will be considered an input (rhs) tensor.
     This method should only be used when the tensor network is fully specified. **/
 void appendTensor(const TensorDenseAdpt<T> & tensor,           //in: new tensor, either input (rhs) or output (lhs)
                   const std::vector<TensorLeg> & connections); //in: connections of the new tensor to other tensors via legs
 /** Appends a tensor to the tensor network by pairing some or all of its
     dimensions with the uncontracted (output) dimensions of the tensor network.
     It is also fine to have none of the tensor legs be contracted with the tensor
     network, in which case they will simply be appended to the output tensor of
     the tensor network. In general, each newly appended tensor removes a number
     of output legs from the tensor network, shifts the numeration of the rest,
     and appends new output legs from itself to the tensor network. **/
 void appendTensor(const TensorDenseAdpt<T> & tensor, //in: tensor being appended to the tensor network
                   const std::vector<std::pair<unsigned int, unsigned int>> & legPairs); //in: leg pairing: pair<tensor network output leg id, tensor leg id>
 /** Appends a rank-2N tensor to a non-empty tensor network by pairing the first
     N legs of the tensor with the specific N output legs of the tensor network,
     subsequently replacing them with the other N legs of the tensor (in order).
     As a result, the number of the output legs of the tensor network won't change. **/
 void appendTensor(const TensorDenseAdpt<T> & Tensor, //in: rank-2N tensor being appended to the tensor network
                   const std::vector<unsigned int> & outLegs); //in: N output legs of the tensor network with which the first N legs of the tensor will be paired
 /** Appends another tensor network into the current tensor network by pairing
     the output legs of both. The remaining output legs of the two tensor networks
     will be placed in order, first tensor network preceding the second one. **/
 void appendNetwork(const TensorNetwork<T> & tensornet, //in: another tensor network
                    const std::vector<std::pair<unsigned int, unsigned int>> & legPairs); //in: leg pairing: pair<output leg id, output leg id>, may be empty
 /** Associates the output (lhs) tensor with its externally provided body (cannot be null). **/
 void setOutputBody(const std::shared_ptr<T> body);
 /** Allocates the output (lhs) tensor body and sets it to zero. **/
 void allocateOutputBody();
 /** Resets the body of an arbitrary tensor. The new body may be null. **/
 void resetTensorBody(const unsigned int tensId, const std::shared_ptr<T> body);
//Transforms:
 /** Contracts two tensors in a tensor network. Always the tensor with a smaller id will be replaced
     by a contracted product while the tensor with a larger id will be deleted from the tensor network,
     causing a shift in the tensor numeration that will affect all tensors with id > "tensId2". **/
 void contractTensors(unsigned int tensId1, //in: id of the 1st tensor in the tensor network: [1..max]
                      unsigned int tensId2, //in: id of the 2nd tensor in the tensor network: [1..max]
                      int * contrPattern = nullptr); //out: digital tensor contraction pattern
 /** Contracts two tensors in a tensor network and returns the result as a raw pointer to the new tensor network.
     Always the tensor with a smaller id will be replaced by a contracted product while the tensor
     with a larger id will be deleted from the tensor network, causing a shift in the tensor numeration
     that will affect all tensors with id > "tensId2". **/
 void contractTensors(const unsigned int tensId1, //in: id of the 1st tensor in the tensor network: [1..max]
                      const unsigned int tensId2, //in: id of the 2nd tensor in the tensor network: [1..max]
                      TensorNetwork<T> ** resultNetwork, //out: tensor network result (returns a pointer to it)
                      int * contrPattern = nullptr) const; //out: digital tensor contraction pattern
 /** Contracts two tensors in a tensor network and returns the result as a smart pointer to the new tensor network.
     Always the tensor with a smaller id will be replaced by a contracted product while the tensor
     with a larger id will be deleted from the tensor network, causing a shift in the tensor numeration
     that will affect all tensors with id > "tensId2". **/
 std::unique_ptr<TensorNetwork<T>> contractTensorsOut(const unsigned int tensId1, //in: id of the 1st tensor in the tensor network: [1..max]
                                                      const unsigned int tensId2, //in: id of the 2nd tensor in the tensor network: [1..max]
                                                      int * contrPattern = nullptr) const; //out: digital tensor contraction pattern
 /** Returns the computational cost of the specified contraction of two tensors. **/
 double getContractionCost(const unsigned int tensId1,         //in: id of the 1st r.h.s. tensor (>0)
                           const unsigned int tensId2,         //in: id of the 2nd r.h.s. tensor (>0)
                           double * arithmIntensity = nullptr, //out: arithmetic intensity
                           bool rescale = false) const;        //in: rescale the Flop cost due to arithmetic intensity
 /** Determines a pseudo-optimal sequence of tensor contractions
     for the given tensor network and numerically evaluates these
     tensor contractions to produce the value of the output tensor.
     If "contrSeq" already contains the previously determined
     contraction sequence, it will be used immediately. **/
 int evaluate(ContractionSequence & contrSeq,                     //inout: tensor contraction sequence (either empty or previously determined)
              const unsigned int numWalkers = NumWalkersDefault); //in: optimization depth
 /** Determines a pseudo-optimal sequence of tensor contractions
     for the given tensor network and numerically evaluates these
     tensor contractions to produce the value of the output tensor
     for which an externally provided body is specified.
     If "contrSeq" already contains the previously determined
     contraction sequence, it will be used immediately. **/
 int evaluate(ContractionSequence & contrSeq,                     //inout: tensor contraction sequence (either empty or previously determined)
              const std::shared_ptr<T> body,                      //in: externally provided pointer to the output tensor body
              const unsigned int numWalkers = NumWalkersDefault); //in: optimization depth

private:
 /** Determines the pseudo-optimal tensor contraction sequence and returns
     it as a vector of pairs of tensor id's to contract. Note that each
     subsequent pair will have its tensor id's refer to the corresponding
     reduced tensor network (tensor numeration changes after each contraction). **/
 void getContractionSequence(ContractionSequence & contrSeq,       //out: contraction sequence
                             const unsigned int numWalkers) const; //in: optimization depth
 /** Performs all tensor contractions, thus evaluating the value of the output tensor.
     Single-node version based on TAL-SH. **/
 int computeOutputLocal(const ContractionSequence & contrSeq); //in: contraction sequence

};

//Template definition:
#include "tensor_network.cpp"

} //end namespace exatensor

#endif //_EXA_TENSOR_NETWORK_H
