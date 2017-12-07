/** C++ adapters for ExaTENSOR: Tensor network

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

//Life cycle:

/** Constructs an empty tensor network **/
template <typename T>
TensorNetwork<T>::TensorNetwork()
{
}

/** Copy constructor. **/
template <typename T>
TensorNetwork<T>::TensorNetwork(const TensorNetwork<T> & tensNetwork):
 Tensors(tensNetwork.Tensors)
{
}

/** Destructor. **/
template <typename T>
TensorNetwork<T>::~TensorNetwork()
{
}

//Accessors:

/** Returns TRUE if the tensor network is empty, FALSE otherwise. **/
template <typename T>
bool TensorNetwork<T>::isEmpty() const
{
 return (Tensors.size() == 0);
}

/** Returns TRUE if the tensor network is closed, FALSE otherwise. **/
template <typename T>
bool TensorNetwork<T>::isClosed() const
{
 return (Tensors[0].getTensorRank() == 0);
}

/** Returns the number of r.h.s. tensors in the tensor network.
    Note that the output (l.h.s.) tensor 0 is not counted here. **/
template <typename T>
unsigned int TensorNetwork<T>::getNumTensors() const
{
 return (unsigned int)(Tensors.size()-1); //not counting the output tensor
}

/** Returns a const reference to a specific tensor from the tensor network. **/
template <typename T>
const TensorDenseAdpt<T> & TensorNetwork<T>::getTensor(const unsigned int id) const
{
 return Tensors.at(id).getTensor();
}

/** Returns a const reference to a specific tensor from the tensor network together with its connections. **/
template <typename T>
const TensorConn<T> & TensorNetwork<T>::getTensorConn(const unsigned int id) const
{
 assert(id <= this->getNumTensors()); //tensor id: [0;1..num_rhs_tensors]
 return Tensors[id];
}

/** Prints. **/
template <typename T>
void TensorNetwork<T>::printIt() const
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
    Each next appended tensor will be considered an input (rhs) tensor.
    This method should only be used when the tensor network is fully specified. **/
template <typename T>
void TensorNetwork<T>::appendTensor(const TensorDenseAdpt<T> & tensor,          //in: new tensor, either input (rhs) or output (lhs)
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
    dimensions with the uncontracted (output) dimensions of the tensor network.
    It is also fine to have none of the tensor legs be contracted with the tensor
    network, in which case they will simply be appended to the output tensor of
    the tensor network. In general, each newly appended tensor removes a number
    of output legs from the tensor network, shifts the numeration of the rest,
    and appends new output legs from itself to the tensor network. **/
template <typename T>
void TensorNetwork<T>::appendTensor(const TensorDenseAdpt<T> & tensor, //in: tensor being appended to the tensor network
                                    const std::vector<std::pair<unsigned int, unsigned int>> & legPairs) //in: leg pairing: pair<tensor network output leg id, tensor leg id>
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
   for(unsigned int i = 0; i < outTensorRank; ++i){
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

/** Appends a rank-2N tensor to a non-empty tensor network by pairing the first
    N legs of the tensor with the specific N output legs of the tensor network,
    subsequently replacing them with the other N legs of the tensor (in order).
    As a result, the number of the output legs of the tensor network won't change. **/
template <typename T>
void TensorNetwork<T>::appendTensor(const TensorDenseAdpt<T> & tensor, //in: rank-2N tensor being appended to the tensor network
                                    const std::vector<unsigned int> & outLegs) //in: N output legs of the tensor network with which the first N legs of the tensor will be paired
{
 assert(this->getNumTensors() > 0); //tensor network must be non-emtpty (at least one input tensor)
 const auto tensRank = tensor.getRank(); //tensor rank
 const auto numLegPairs = outLegs.size(); //number of legs to pair
 auto & outTensor = Tensors[0]; //output tensor of the tensor network
 const auto outTensorRank = outTensor.getTensorRank(); //output tensor rank
 assert(numLegPairs <= outTensorRank && numLegPairs*2 == tensRank);
 unsigned int braLegId = 0; //will cover bra tensor legs (first N tensor legs)
 unsigned int ketLegId = numLegPairs; //will cover ket tensor legs (second N tensor legs)
 unsigned int newTensorId = Tensors.size(); //id of the the newly added tensor in the tensor network
 std::vector<TensorLeg> tensorLegs; //will be the legs of the newly appended tensor
 for(auto outLegId: outLegs){ //outLeg: output tensor leg id to pair with
  const auto oldOutLeg = outTensor.getTensorLeg(outLegId); //output tensor leg to pair with
  const auto inpTensId = oldOutLeg.getTensorId(); //input tensor id the output leg was paired with
  const auto inpTensLegId = oldOutLeg.getDimensionId(); //input tensor dimension id the output leg was paired with
  Tensors[inpTensId].resetConnection(inpTensLegId,TensorLeg(newTensorId,braLegId++)); //re-pair the input tensor leg with the newly appended tensor
  tensorLegs.emplace_back(TensorLeg(inpTensId,inpTensLegId)); //pair the newly appended tensor with the input tensor leg which was previously paired with the output tensor
  outTensor.resetConnection(outLegId,TensorLeg(newTensorId,ketLegId++)); //re-pair output tensor leg with an uncontracted leg of the newly appended tensor
 }
 for(auto outLegId: outLegs) tensorLegs.emplace_back(TensorLeg(0,outLegId)); //append ket (output connected) legs to the newly appended tensor
 Tensors.emplace_back(TensorConn<T>(tensor,tensorLegs)); //append the new tensor into the tensor network
 return;
}

/** Appends another tensor network into the current tensor network by pairing
    the output legs of both. The remaining output legs of the two tensor networks
    will be placed in order, first tensor network preceding the second one. **/
template <typename T>
void TensorNetwork<T>::appendNetwork(const TensorNetwork<T> & tensornet, //in: another tensor network
                                     const std::vector<std::pair<unsigned int, unsigned int>> & legPairs) //in: leg pairing: pair<output leg id, output leg id>, may be empty
{
 assert(!(this->isEmpty()) && !(tensornet.isEmpty()));
 const unsigned int numConnections = legPairs.size(); //number of connections between the two tensor networks
 const auto numTensors1 = this->getNumTensors(); //number of the r.h.s. tensors in the 1st tensor network
 const auto numTensors2 = tensornet.getNumTensors(); //number the r.h.s. tensors in the 2nd tensor network
 auto & outTensorConn1 = Tensors[0]; //connected output tensor of the 1st tensor network (mutable reference)
 const unsigned int numOutLegs1 = outTensorConn1.getNumLegs(); //original number of legs in the 1st output tensor
 const auto & outTensorConn2 = tensornet.getTensorConn(0); //connected output tensor of the 2nd tensor network (immutable reference)
 const unsigned int numOutLegs2 = outTensorConn2.getNumLegs(); //original number of legs in the 2st output tensor
 assert(numConnections <= numOutLegs1 && numConnections <= numOutLegs2);
 //Append new r.h.s. tensors into the 1st tensor network and shift tensor numeration (numeration starts from numTensors1+1):
 assert(Tensors.size() == numTensors1 + 1);
 for(unsigned int i = 1; i <= numTensors2; ++i){
  Tensors.push_back(tensornet.getTensorConn(i)); //copy r.h.s. tensors from the 2nd tensor network into the 1st one
  auto & tensConn = Tensors[numTensors1+i]; //newly appended tensor (mutable reference)
  const auto numLegs = tensConn.getNumLegs(); //number of legs in the newly appended tensor
  for(unsigned int j = 0; j < numLegs; ++j){ //loop over the legs of the newly appended tensor
   auto leg = tensConn.getTensorLeg(j); //leg of the newly appended tensor (mutable copy)
   const auto oldTensId = leg.getTensorId(); //id of the tensor the leg is connected with
   if(oldTensId > 0){ //shift numeration for r.h.s. tensors only
    leg.resetTensorId(numTensors1+oldTensId); //numeration is shifted by the number of r.h.s. tensors in the 1st tensor network
    tensConn.resetConnection(j,leg); //update the leg of the newly appended tensor
   }
  }
 }
 //Contract output tensor legs and remove them from the 1st output tensor, update input tensor legs:
 unsigned int outLegActive[numOutLegs2+1] = {1}; //mark legs of the 2nd output tensor as active (one additional element for avoiding zero length)
 for(auto & legPair: legPairs){ //match contracted output legs
  const auto legId1 = legPair.first; //corresponding leg id in the 1st output tensor
  const auto legId2 = legPair.second; //corresponding leg id in the 2nd output tensor
  assert(legId1 < numOutLegs1 && legId2 < numOutLegs2);
  outLegActive[legId2] = 0; //deactivate the contracted leg from the 2nd output tensor
  auto outLeg1 = outTensorConn1.getTensorLeg(legId1); //corresponding leg of the 1st output tensor (mutable copy)
  const auto & outLeg2 = outTensorConn2.getTensorLeg(legId2); //corresponding leg of the 2nd output tensor
  const auto rhsTensId1 = outLeg1.getTensorId(); //id of the rhs tensor connected to outLeg1
  const auto rhsTensLegId1 = outLeg1.getDimensionId(); //dimension of the rhs tensor connected to outLeg1
  const auto rhsTensId2 = outLeg2.getTensorId(); //id of the rhs tensor connected to outLeg2
  const auto rhsTensLegId2 = outLeg2.getDimensionId(); //dimension of the rhs tensor connected to outLeg2
  Tensors[rhsTensId1].resetConnection(rhsTensLegId1,TensorLeg(numTensors1+rhsTensId2,rhsTensLegId2));
  Tensors[numTensors1+rhsTensId2].resetConnection(rhsTensLegId2,TensorLeg(rhsTensId1,rhsTensLegId1));
  outLeg1.resetTensorId(0); //special tag to delete this leg later
  outTensorConn1.resetConnection(legId1,outLeg1); //replace the leg of the output tensor with a dead leg to later delete it
 }
 //Amputate dead legs from the 1st output tensor:
 unsigned int numDeletedLegs = 0;
 for(unsigned int i = 0; i < numOutLegs1; ++i){
  const unsigned int j = i - numDeletedLegs; //new (shifted) id of the output tensor leg
  const auto & outLeg = outTensorConn1.getTensorLeg(j);
  const auto tensId = outLeg.getTensorId();
  if(tensId == 0){ //dead leg (delete it)
   outTensorConn1.deleteDimension(j);
   numDeletedLegs++;
  }else{ //live leg (update the connected r.h.s. tensor leg)
   const auto tensLegId = outLeg.getDimensionId();
   Tensors[tensId].resetConnection(tensLegId,TensorLeg(0,j)); //update connection of the r.h.s. tensor to the output tensor
  }
 }
 //Append the remaining tensor legs from the 2nd output tensor to the 1st output tensor and update r.h.s. connections:
 auto numOutLegs = outTensorConn1.getNumLegs(); //number of output legs left in the 1st output tensor
 for(unsigned int i = 0; i < numOutLegs2; ++i){ //loop over the output legs of the 2nd output tensor
  if(outLegActive[i] != 0){ //leg is still active (uncontracted)
   auto newLeg = outTensorConn2.getTensorLeg(i); //uncontracted leg from the 2nd output tensor (mutable copy)
   const auto rhsTensId2 = newLeg.getTensorId();
   const auto rhsTensLegId2 = newLeg.getDimensionId();
   newLeg.resetTensorId(numTensors1+rhsTensId2);
   outTensorConn1.appendDimension(outTensorConn2.getDimExtent(i),newLeg);
   Tensors[numTensors1+rhsTensId2].resetConnection(rhsTensLegId2,TensorLeg(0,numOutLegs++));
  }
 }
 assert(numOutLegs == ((numOutLegs1 + numOutLegs2) - (numConnections * 2)));
 return;
}

/** Associates the output (lhs) tensor with its externally provided body (cannot be null). **/
template <typename T>
void TensorNetwork<T>::setOutputBody(const std::shared_ptr<T> body)
{
#ifdef _DEBUG_DIL
 assert(Tensors.size() > 0 && body);
#endif
 Tensors[0].setBody(body);
 return;
}

/** Allocates the output (lhs) tensor body and sets it to zero. **/
template <typename T>
void TensorNetwork<T>::allocateOutputBody()
{
#ifdef _DEBUG_DIL
 assert(Tensors.size() > 0);
#endif
 Tensors[0].allocateBody();
 Tensors[0].nullifyBody();
 return;
}

/** Resets the body of an arbitrary tensor. The new body may be null. **/
template <typename T>
void TensorNetwork<T>::resetTensorBody(const unsigned int tensId, const std::shared_ptr<T> body)
{
#ifdef _DEBUG_DIL
 assert(tensId < Tensors.size());
#endif
 Tensors[tensId].resetBody(body);
 return;
}

//Transforms:

/** Contracts two tensors in a given tensor network. Always the tensor with a smaller id will be replaced
    by a contracted product while the tensor with a larger id will be deleted from the tensor network,
    causing a shift in the tensor numeration that will affect all tensors with id > "tensId2". **/
template <typename T>
void TensorNetwork<T>::contractTensors(unsigned int tensId1, //in: id of the 1st tensor in the tensor network: [1..max]
                                       unsigned int tensId2, //in: id of the 2nd tensor in the tensor network: [1..max]
                                       int * contrPattern)   //out: digital tensor contraction pattern
{
#ifdef _DEBUG_DIL
 assert((tensId1 >= 1 && tensId1 <= this->getNumTensors()) && (tensId2 >= 1 && tensId2 <= this->getNumTensors()));
 assert(tensId1 != tensId2);
#endif
 if(tensId1 > tensId2){auto tmp = tensId1; tensId1 = tensId2; tensId2 = tmp;}
 unsigned int cp = 0; unsigned int un = 0;
 //Remove contracted legs from the 1st tensor:
 auto & tensor1 = Tensors[tensId1]; //1st contracted tensor
 auto rank1 = tensor1.getTensorRank(); //rank of the 1st tensor
 unsigned int j = 0; //number of removed legs
 for(unsigned int i = 0; i < rank1; ++i){
  auto legId = i - j;
  const auto & leg = tensor1.getTensorLeg(legId);
  const auto connTensId = leg.getTensorId();
  const auto connTensLegId = leg.getDimensionId();
  if(connTensId == tensId2){ //contracted dimension (remove)
   tensor1.deleteDimension(legId); ++j;
   if(contrPattern != nullptr) contrPattern[cp++] = -(connTensLegId + 1);
  }else{
   Tensors[connTensId].resetConnection(connTensLegId,TensorLeg(tensId1,legId));
   if(contrPattern != nullptr) contrPattern[cp++] = ++un;
  }
 }
 rank1 = tensor1.getTensorRank(); //updated rank of the 1st tensor
 //Remove contracted legs from the 2nd tensor:
 auto & tensor2 = Tensors[tensId2]; //2nd contracted tensor
 auto rank2 = tensor2.getTensorRank(); //rank of the 2nd tensor
 j = 0; //number of removed legs
 for(unsigned int i = 0; i < rank2; ++i){
  auto legId = i - j;
  const auto & leg = tensor2.getTensorLeg(legId);
  const auto connTensId = leg.getTensorId();
  const auto connTensLegId = leg.getDimensionId();
  if(connTensId == tensId1){ //contracted dimension (remove)
   tensor2.deleteDimension(legId); ++j;
   if(contrPattern != nullptr) contrPattern[cp++] = -(connTensLegId + 1);
  }else{
   Tensors[connTensId].resetConnection(connTensLegId,TensorLeg(tensId2,legId));
   if(contrPattern != nullptr) contrPattern[cp++] = ++un;
  }
 }
 rank2 = tensor2.getTensorRank(); //updated rank of the 2nd tensor
 //Append uncontracted legs of the 2nd tensor to the 1st tensor:
 for(unsigned int i = 0; i < rank2; ++i){
  const auto & leg = tensor2.getTensorLeg(i);
  const auto connTensId = leg.getTensorId();
  const auto connTensLegId = leg.getDimensionId();
  Tensors[connTensId].resetConnection(connTensLegId,TensorLeg(tensId1,rank1));
  tensor1.appendDimension(tensor2.getDimExtent(i),leg);
  rank1 = tensor1.getTensorRank(); //updated rank of the 1st tensor
 }
 //Delete the 2nd tensor and adjust tensor numeration in the tensor network:
 Tensors.erase(Tensors.begin()+tensId2);
 for(j = 0; j < Tensors.size(); ++j){
  auto & tensor = Tensors[j];
  for(unsigned int i = 0; i < tensor.getNumLegs(); ++i){
   const auto & leg = tensor.getTensorLeg(i);
   auto connTensId = leg.getTensorId();
#ifdef _DEBUG_DIL
   assert(connTensId != tensId2);
#endif
   if(connTensId > tensId2) tensor.resetConnection(i,TensorLeg(connTensId-1,leg.getDimensionId()));
  }
 }
 return;
}

/** Contracts two tensors in a tensor network and returns the result as a raw pointer to a new tensor network.
    Always the tensor with a smaller id will be replaced by a contracted product while the tensor
    with a larger id will be deleted from the tensor network, causing a shift in the tensor numeration
    that will affect all tensors with id > "tensId2". **/
template <typename T>
void TensorNetwork<T>::contractTensors(const unsigned int tensId1, //in: id of the 1st tensor in the tensor network: [1..max]
                                       const unsigned int tensId2, //in: id of the 2nd tensor in the tensor network: [1..max]
                                       TensorNetwork<T> ** resultNetwork, //out: tensor network result (returns a pointer to it)
                                       int * contrPattern) const          //out: digital tensor contraction pattern
{
 *resultNetwork = new TensorNetwork<T>(*this);
 (*resultNetwork)->contractTensors(tensId1,tensId2,contrPattern);
 return;
}

/** Contracts two tensors in a tensor network and returns the result as a smart pointer to a new tensor network.
    Always the tensor with a smaller id will be replaced by a contracted product while the tensor
    with a larger id will be deleted from the tensor network, causing a shift in the tensor numeration
    that will affect all tensors with id > "tensId2". **/
template <typename T>
std::unique_ptr<TensorNetwork<T>> TensorNetwork<T>::contractTensorsOut(const unsigned int tensId1, //in: id of the 1st tensor in the tensor network: [1..max]
                                                                       const unsigned int tensId2, //in: id of the 2nd tensor in the tensor network: [1..max]
                                                                       int * contrPattern) const   //out: digital tensor contraction pattern
{
 std::unique_ptr<TensorNetwork<T>> resultNetwork(new TensorNetwork<T>(*this));
 resultNetwork->contractTensors(tensId1,tensId2,contrPattern);
 return std::move(resultNetwork);
}

/** Returns the computational cost of the specified contraction of two tensors in a tensor network. **/
template <typename T>
double TensorNetwork<T>::getContractionCost(const unsigned int tensId1, //in: id of the 1st rhs tensor (>0)
                                            const unsigned int tensId2, //in: id of the 2nd rhs tensor (>0)
                                            double * arithmIntensity,   //out: arithmetic intensity
                                            bool rescale) const         //in: rescale the Flop cost due to arithmetic intensity
{
#ifdef _DEBUG_DIL
 assert((tensId1 >= 1 && tensId1 <= this->getNumTensors()) && (tensId2 >= 1 && tensId2 <= this->getNumTensors()));
 assert(tensId1 != tensId2);
#endif
 //Compute Flop count:
 const auto & tensor1 = Tensors[tensId1];
 double lTensVol = static_cast<double>(tensor1.getVolume());
 double cost = lTensVol;
 double cVol = 1.0;
 double lVol = 1.0;
 double rVol = 1.0;
 const auto & tensor2 = Tensors[tensId2];
 const auto rank2 = tensor2.getTensorRank();
 for(unsigned int legId = 0; legId < rank2; ++legId){
  const auto & leg = tensor2.getTensorLeg(legId);
  double legVol = (static_cast<double>(tensor2.getDimExtent(legId)));
  if(leg.getTensorId() != tensId1){ //right leg
   rVol*=legVol; cost*=legVol;
  }else{ //contracted leg
   cVol*=legVol;
  }
 }
 lVol = lTensVol / cVol;
 double rTensVol = cVol * rVol;
 double dTensVol = lVol * rVol;
 double arithInt = cost / (dTensVol + lTensVol + rTensVol);
 if(arithmIntensity != nullptr) *arithmIntensity = arithInt;
 //Rescale the cost due to arithmetic intensity:
 if(rescale){
  //`Finish
 }
 return cost;
}

/** Determines a pseudo-optimal sequence of tensor contractions
    for the given tensor network and numerically evaluates these
    tensor contractions to produce the value of the output tensor.
    If "contrSeq" already contains the previously determined
    contraction sequence, it will be used immediately. **/
template <typename T>
int TensorNetwork<T>::evaluate(ContractionSequence & contrSeq,
                               const unsigned int numWalkers)
{
 int error_code = 0; //success
 auto numTensors = this->getNumTensors(); //number of r.h.s. tensors in the tensor network
 auto numContr = contrSeq.size(); //number of tensor contractions in the contraction sequence
 if(numContr == 0){ //contraction sequence has not been determined yet
  this->getContractionSequence(contrSeq,numWalkers);
 }else{
  if(numContr != (numTensors - 1)) error_code=-1; //invalid number of tensor contractions in the contraction sequence
 }
 if(error_code == 0) error_code = this->computeOutputLocal(contrSeq);
 return error_code;
}

/** Determines a pseudo-optimal sequence of tensor contractions
    for the given tensor network and numerically evaluates these
    tensor contractions to produce the value of the output tensor
    for which an externally provided body is specified.
    If "contrSeq" already contains the previously determined
    contraction sequence, it will be used immediately. **/
template <typename T>
int TensorNetwork<T>::evaluate(ContractionSequence & contrSeq,
                               const std::shared_ptr<T> body,
                               const unsigned int numWalkers)
{
 int error_code = 0; //success
 this->setOutputBody(body);
 error_code = this->evaluate(contrSeq,numWalkers);
 return error_code;
}

/** Determines the pseudo-optimal tensor contraction sequence and returns
    it as a vector of pairs of the r.h.s. tensor id's to contract. Note that
    each subsequent pair will have its tensor id's refer to the corresponding
    reduced tensor network. **/
template <typename T>
void TensorNetwork<T>::getContractionSequence(ContractionSequence & contrSeq,
                                              const unsigned int numWalkers) const
{
 using ContrPath = std::tuple<TensorNetwork<T>,ContractionSequence,double>;

 std::cout << "#MSG(TensorNetwork<T>::getContractionSequence): Determining a pseudo-optimal tensor contraction sequence ... "; //debug

 auto timeBeg = std::chrono::high_resolution_clock::now();

 auto numContractions = this->getNumTensors() - 1; //number of contractions is one less than the number of r.h.s. tensors
 assert(numContractions > 0); //at least one tensor contraction is expected (two or more r.h.s. tensors)
 assert(contrSeq.size() == 0); //the contraction sequence must be empty on entrance

 ContractionSequence contrSeqEmpty;
 std::vector<ContrPath> inputPaths; //input: vector
 inputPaths.emplace_back(std::make_tuple(*this,contrSeqEmpty,0.0)); //initial configuration

 auto cmpPaths = [](const ContrPath & left, const ContrPath & right){return (std::get<2>(left) < std::get<2>(right));};
 std::priority_queue<ContrPath,std::vector<ContrPath>,decltype(cmpPaths)> priq(cmpPaths); //output: priority queue

 for(decltype(numContractions) pass = 0; pass < numContractions; ++pass){
  unsigned int numPassCands = 0;
  for(const auto & contrPath : inputPaths){
   const auto & parentTensNet = std::get<0>(contrPath); //parental tensor network
   const auto numTensors = parentTensNet.getNumTensors(); //number of r.h.s. tensors in the parental tensor network
   const auto & parentContrSeq = std::get<1>(contrPath); //parental contraction sequence
   for(unsigned int i = 1; i < numTensors; ++i){ //r.h.s. tensors are numbered from 1
    for(unsigned int j=i+1; j <= numTensors; ++j){ //r.h.s. tensors are numbered from 1
     double contrCost = parentTensNet.getContractionCost(i,j); //tensor contraction cost
     //std::cout << std::endl << "New candidate contracted pair of tensors is {" << i << "," << j << "} with cost " << contrCost; //debug
     auto tensNet = parentTensNet.contractTensorsOut(i,j); //contract tensors i and j and return a pointer to a new tensor network
     auto cSeq = parentContrSeq; cSeq.emplace_back(std::pair<unsigned int, unsigned int>(i,j)); //append a new pair of contracted tensor id's
     priq.emplace(std::make_tuple(*tensNet,cSeq,contrCost+std::get<2>(contrPath))); //cloning tensor network and contraction sequence
     if(priq.size() > numWalkers) priq.pop();
     numPassCands++;
    }
   }
  }
  std::cout << std::endl << "Pass " << pass << ": Total number of candidates considered = " << numPassCands; //debug
  inputPaths.clear();
  if(pass == numContractions - 1){ //last pass
   while(priq.size() > 1) priq.pop();
   contrSeq = std::get<1>(priq.top());
   std::cout << std::endl << "Best tensor contraction sequence cost found = " << std::get<2>(priq.top()); //debug
   priq.pop();
  }else{
   while(priq.size() > 0){
    inputPaths.emplace_back(priq.top());
    priq.pop();
   }
  }
 }

 auto timeEnd = std::chrono::high_resolution_clock::now();
 auto timeTot = std::chrono::duration_cast<std::chrono::duration<double>>(timeEnd-timeBeg);
 std::cout << std::endl << "Done (" << timeTot.count() << " sec):"; //debug
 for(const auto & cPair : contrSeq) std::cout << " {" << std::get<0>(cPair) << "," << std::get<1>(cPair) << "}"; //debug
 std::cout << std::endl; //debug
 return;
}

/** Performs all tensor contractions, thus evaluating the value of the output tensor.
    Single-node version based on TAL-SH. **/
template <typename T>
int TensorNetwork<T>::computeOutputLocal(const ContractionSequence & contrSeq)
{
 talsh_task_t tsk;
 talsh_tens_t *tens, *dtens, *ltens, *rtens;
 int errc,cpl,lRank,rRank,tDims[MAX_TENSOR_RANK],contrPtrnDig[MAX_TENSOR_RANK*2];
 char contrPtrnSym[512]; //should be large enough to contain an arbitrary binary tensor contraction specification

 int error_code = 0; //success
 std::cout << "#MSG(TensorNetwork<T>::computeOutputLocal): Computing ... "; //debug
 auto timeBeg = std::chrono::high_resolution_clock::now();

 std::map<std::uintptr_t,std::pair<talsh_tens_t*,std::shared_ptr<T>>> tensMap; //body address --> {TAL-SH tensor pointer, tensor body pointer}

 auto numContractions = contrSeq.size();
 assert(numContractions == (this->getNumTensors() - 1));
 std::unique_ptr<TensorNetwork<T>> tensNet(this);
 for(decltype(numContractions) contrNum = 0; contrNum < numContractions; ++contrNum){

  //Extract the pair of contracted input tensors:
  auto lid = std::get<0>(contrSeq[contrNum]); //left input tensor id
  auto rid = std::get<1>(contrSeq[contrNum]); //right input tensor id
  assert(lid > 0 && lid < rid); //r.h.s. tensor id is always > 0 (0 is the output tensor)

  //Construct the left tensor:
  auto & leftTensor = tensNet->getTensor(lid); //left tensor
  int tRank = leftTensor.getRank(); lRank=tRank;
  auto pDims = leftTensor.getDimExtents();
  for(int i = 0; i < tRank; ++i) tDims[i] = static_cast<int>(pDims[i]);
  auto pBody = leftTensor.getBodyAccess();
  void * tBody = static_cast<void*>(pBody.get());
  assert(tBody != nullptr); //input tensor must have been defined
  std::uintptr_t lbodyAddr = reinterpret_cast<std::uintptr_t>(tBody);
  auto tensPos = tensMap.find(lbodyAddr);
  if(tensPos != tensMap.end()){ //intermediate tensor
   tens = (tensPos->second).first;
  }else{ //external input tensor
   errc = talshTensorCreate(&tens); if(errc != TALSH_SUCCESS) return -1;
   errc = talshTensorConstruct(tens,TensorDataKind<T>::Type,tRank,tDims,talshFlatDevId(DEV_HOST,0),tBody);
   if(errc != TALSH_SUCCESS) return -1;
   tensMap.emplace(lbodyAddr,std::pair<talsh_tens_t*,std::shared_ptr<T>>(tens,std::shared_ptr<T>(nullptr,[](T * ptr){})));
  }
  ltens = tens;

  //Construct the right tensor:
  auto & rightTensor = tensNet->getTensor(rid); //right tensor
  tRank = rightTensor.getRank(); rRank = tRank;
  pDims = rightTensor.getDimExtents();
  for(int i = 0; i < tRank; ++i) tDims[i] = static_cast<int>(pDims[i]);
  pBody = rightTensor.getBodyAccess();
  tBody = static_cast<void*>(pBody.get());
  assert(tBody != nullptr); //input tensor must have been defined
  std::uintptr_t rbodyAddr = reinterpret_cast<std::uintptr_t>(tBody);
  tensPos = tensMap.find(rbodyAddr);
  if(tensPos != tensMap.end()){ //intermediate tensor
   tens = (tensPos->second).first;
  }else{ //external input tensor
   errc = talshTensorCreate(&tens); if(errc != TALSH_SUCCESS) return -1;
   errc = talshTensorConstruct(tens,TensorDataKind<T>::Type,tRank,tDims,talshFlatDevId(DEV_HOST,0),tBody);
   if(errc != TALSH_SUCCESS) return -1;
   tensMap.emplace(rbodyAddr,std::pair<talsh_tens_t*,std::shared_ptr<T>>(tens,std::shared_ptr<T>(nullptr,[](T * ptr){})));
  }
  rtens = tens;

  std::cout << std::endl << "Left body  = " << lbodyAddr << ": Ref count = " << tensMap[lbodyAddr].second.use_count(); //debug
  std::cout << std::endl << "Right body = " << rbodyAddr << ": Ref count = " << tensMap[rbodyAddr].second.use_count(); //debug
  //std::string tmp; std::cout << std::endl << "Press any key and enter to continue ..."; std::cin >> tmp; //debug

  //Construct the destination tensor:
  decltype(lid) did; //destination tensor id
  if(contrNum == numContractions - 1){ //last contraction
   assert(lid == 1 && rid == 2); // last contraction: 0+=1*2
   auto tmpTensNet = tensNet->contractTensorsOut(lid,rid,contrPtrnDig);
   did = 0; //output tensor
  }else{
   auto tmpTensNet = tensNet->contractTensorsOut(lid,rid,contrPtrnDig);
   if(contrNum == 0) auto pThis = tensNet.release(); //release "this"
   tensNet = std::move(tmpTensNet);
   did = lid;
  }
  auto & resultTensor = tensNet->getTensor(did); //destination tensor
  tRank = resultTensor.getRank();
  pDims = resultTensor.getDimExtents();
  for(int i = 0; i < tRank; ++i) tDims[i] = static_cast<int>(pDims[i]);
  errc = talshTensorCreate(&tens); if(errc != TALSH_SUCCESS) return -1;
  std::shared_ptr<T> sBody(nullptr);
  pBody = resultTensor.getBodyAccess();
  if(pBody){ //output tensor (already has an external body)
   tBody = static_cast<void*>(pBody.get());
   errc = talshTensorConstruct(tens,TensorDataKind<T>::Type,tRank,tDims,talshFlatDevId(DEV_HOST,0),tBody);
   if(errc != TALSH_SUCCESS) return -1;
  }else{ //intermediate tensor (body is owned by TAL-SH)
   errc = talshTensorConstruct(tens,TensorDataKind<T>::Type,tRank,tDims,talshFlatDevId(DEV_HOST,0));
   if(errc != TALSH_SUCCESS) return -1;
   errc = talshTensorGetBodyAccess(tens,&tBody,TensorDataKind<T>::Type,0,DEV_HOST);
   if(errc != TALSH_SUCCESS) return -1;
   sBody = std::shared_ptr<T>(static_cast<T*>(tBody),[](T * ptr){});
   tensNet->resetTensorBody(did,sBody);
  }
  std::uintptr_t dbodyAddr = reinterpret_cast<std::uintptr_t>(tBody);
  tensMap.emplace(dbodyAddr,std::pair<talsh_tens_t*,std::shared_ptr<T>>(tens,sBody));
  dtens = tens;
  std::cout << std::endl << "Intermediate body = " << dbodyAddr << ": Ref count = " << tensMap[dbodyAddr].second.use_count(); //debug

  //Get the symbolic tensor contraction pattern:
  get_contr_pattern_sym(&lRank,&rRank,contrPtrnDig,contrPtrnSym,&cpl,&errc); if(errc != TALSH_SUCCESS) return -1;

  //Perform the tensor contraction:
  errc = talshTensorContract(contrPtrnSym,dtens,ltens,rtens); if(errc != TALSH_SUCCESS) return -1;

  //Destruct processed TAL-SH tensor aliases:
  // Right input tensor:
  tens=nullptr; tens=(tensMap[rbodyAddr]).first;
  assert(tens != nullptr);
  errc = talshTensorDestroy(tens); if(errc != TALSH_SUCCESS) return -1;
  tensMap.erase(rbodyAddr);
  // Left input tensor:
  tens=nullptr; tens=(tensMap[lbodyAddr]).first;
  assert(tens != nullptr);
  errc = talshTensorDestroy(tens); if(errc != TALSH_SUCCESS) return -1;
  tensMap.erase(lbodyAddr);
  // Output tensor:
  if(contrNum == numContractions - 1){ //last tensor contraction
   tens=nullptr; tens=(tensMap[dbodyAddr]).first;
   assert(tens != nullptr);
   errc = talshTensorDestroy(tens); if(errc != TALSH_SUCCESS) return -1;
   tensMap.erase(dbodyAddr);
  }

  //std::cout << std::endl << "End of iteration " << contrNum << ". Press any key and enter to continue ..."; std::cin >> tmp; //debug

 }

 auto timeEnd = std::chrono::high_resolution_clock::now();
 auto timeTot = std::chrono::duration_cast<std::chrono::duration<double>>(timeEnd-timeBeg);
 std::cout << std::endl << "Done (" << timeTot.count() << " sec)" << std::endl; //debug
 return error_code;
}
