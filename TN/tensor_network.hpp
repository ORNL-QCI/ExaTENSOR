/** C++ adapters for ExaTENSOR: Tensor network

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/09/08

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
#include <assert.h>
#include <iostream>
#include <ctime>
#include <chrono>

#include "type_deduct.hpp"

#include "tensor_conn.hpp"

#include "talsh.h"

#define _DEBUG_DIL

namespace exatensor {

//Types:
using ContractionSequence = std::vector<std::pair<unsigned int, unsigned int>>;

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

public:

//Constants:
 static const unsigned int NumWalkersDefault = 1024; //number of walkers for tensor contraction sequence optimization

//Life cycle:
 /** Constructs an empty tensor network **/
 TensorNetwork();
 /** Copy constructor. **/
 TensorNetwork(const TensorNetwork<T> & tensNetwork);
 /** Destructor. **/
 virtual ~TensorNetwork();

//Accessors:
 /** Returns the number of r.h.s. tensors in the tensor network.
     Note that the output (l.h.s.) tensor 0 is not counted here. **/
 unsigned int getNumTensors() const;
 /** Returns a const reference to a specific tensor from the tensor network. **/
 const TensorDenseAdpt<T> & getTensor(const unsigned int id) const;
 /** Prints. **/
 void printIt() const;

//Mutators:
 /** Explicitly appends a tensor to the tensor network, either input or
     output. The output (lhs) tensor must be appended first (tensor 0).
     Each next appended tensor will be considered an input (rhs) tensor. **/
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
 /** Appends another tensor network into the current tensor network
     by pairing the dimensions of the output tensors of both. **/
 void appendNetwork(const TensorNetwork<T> & tensornet, //in: another tensor network
                    const std::vector<std::pair<unsigned int, unsigned int>> & legPairs); //in: leg pairing: pair<output leg id, output leg id>, may be empty
 /** Associates the output (lhs) tensor with its externally provided body. **/
 void setOutputBody(const std::shared_ptr<T> body);
 /** Resets the body of an arbitrary tensor. The new body may be null. **/
 void resetTensorBody(const unsigned int tensId, const std::shared_ptr<T> body);
//Transforms:
 /** Contracts two tensors in a tensor network. Always the tensor with a smaller id will be replaced
     by a contracted product while the tensor with a larger id will be deleted from the tensor network,
     causing a shift in the tensor numeration that will affect all tensors with id > "tensId2". **/
 void contractTensors(unsigned int tensId1,  //in: id of the 1st tensor in the tensor network: [1..max]
                      unsigned int tensId2); //in: id of the 2nd tensor in the tensor network: [1..max]
 /** Contracts two tensors in a tensor network and returns the result as a raw pointer to the new tensor network.
     Always the tensor with a smaller id will be replaced by a contracted product while the tensor
     with a larger id will be deleted from the tensor network, causing a shift in the tensor numeration
     that will affect all tensors with id > "tensId2". **/
 void contractTensors(const unsigned int tensId1, //in: id of the 1st tensor in the tensor network: [1..max]
                      const unsigned int tensId2, //in: id of the 2nd tensor in the tensor network: [1..max]
                      TensorNetwork<T> ** resultNetwork) const; //out: tensor network result (returns a pointer to it)
 /** Contracts two tensors in a tensor network and returns the result as a smart pointer to the new tensor network.
     Always the tensor with a smaller id will be replaced by a contracted product while the tensor
     with a larger id will be deleted from the tensor network, causing a shift in the tensor numeration
     that will affect all tensors with id > "tensId2". **/
 std::unique_ptr<TensorNetwork<T>> contractTensorsOut(const unsigned int tensId1,        //in: id of the 1st tensor in the tensor network: [1..max]
                                                      const unsigned int tensId2) const; //in: id of the 2nd tensor in the tensor network: [1..max]
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
 /** Performs all tensor contractions, thus evaluating the value of the output tensor. **/
 int computeOutput(const ContractionSequence & contrSeq);

};

//Template definition:
#include "tensor_network.cpp"

} //end namespace exatensor

#endif //_EXA_TENSOR_NETWORK_H
