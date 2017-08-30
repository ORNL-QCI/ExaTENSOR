/** C++ adapters for ExaTENSOR: Tensor network

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/08/30

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
    Each next appended tensor will be considered an input (rhs) tensor. **/
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
    dimensions with the uncontracted dimensions of the tensor network.
    It is also fine to have none of the tensor legs be contracted with
    the tensor network, in which case they will simply be appended to
    the output tensor of the tensor network. **/
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

/** Appends another tensor network into the current tensor network
    by pairing the dimensions of the output tensors of both. **/
template <typename T>
void TensorNetwork<T>::appendNetwork(const TensorNetwork<T> & tensornet, //in: another tensor network
                                     const std::vector<std::pair<unsigned int, unsigned int>> & legPairs) //in: leg pairing: pair<output leg id, output leg id>, may be empty
{
 //`Finish
 return;
}

/** Associates the output (lhs) tensor with its externally provided body. **/
template <typename T>
void TensorNetwork<T>::setOutputBody(const std::shared_ptr<T> body)
{
#ifdef _DEBUG_DIL
 assert(Tensors.size() > 0 && body);
#endif
 Tensors[0].setBody(body);
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
                                       unsigned int tensId2) //in: id of the 2nd tensor in the tensor network: [1..max]
{
#ifdef _DEBUG_DIL
 assert((tensId1 >= 1 && tensId1 <= this->getNumTensors()) && (tensId2 >= 1 && tensId2 <= this->getNumTensors()));
 assert(tensId1 != tensId2);
#endif
 if(tensId1 > tensId2){auto tmp = tensId1; tensId1 = tensId2; tensId2 = tmp;}
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
   tensor1.deleteDimension(legId);
   ++j;
  }else{
   Tensors[connTensId].resetConnection(connTensLegId,TensorLeg(tensId1,legId));
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
   tensor2.deleteDimension(legId);
   ++j;
  }else{
   Tensors[connTensId].resetConnection(connTensLegId,TensorLeg(tensId2,legId));
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
                                       TensorNetwork<T> ** resultNetwork) const //out: tensor network result (returns a pointer to it)
{
 *resultNetwork = new TensorNetwork<T>(*this);
 (*resultNetwork)->contractTensors(tensId1,tensId2);
 return;
}

/** Contracts two tensors in a tensor network and returns the result as a smart pointer to a new tensor network.
    Always the tensor with a smaller id will be replaced by a contracted product while the tensor
    with a larger id will be deleted from the tensor network, causing a shift in the tensor numeration
    that will affect all tensors with id > "tensId2". **/
template <typename T>
std::unique_ptr<TensorNetwork<T>> TensorNetwork<T>::contractTensorsOut(const unsigned int tensId1,       //in: id of the 1st tensor in the tensor network: [1..max]
                                                                       const unsigned int tensId2) const //in: id of the 2nd tensor in the tensor network: [1..max]
{
 std::unique_ptr<TensorNetwork<T>> resultNetwork(new TensorNetwork<T>(*this));
 resultNetwork->contractTensors(tensId1,tensId2);
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
 if(error_code == 0) error_code = this->computeOutput(contrSeq);
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

 auto numContractions = this->getNumTensors() - 1; //number of contractions is one less than the number of r.h.s. tensors
 assert(numContractions > 0); //at least one tensor contraction is expected (two or more r.h.s. tensors)
 assert(contrSeq.size() == 0); //the contraction sequence must be empty on entrance

 ContractionSequence contrSeqEmpty;
 std::vector<ContrPath> inputPaths;
 inputPaths.emplace_back(make_tuple(*this,contrSeqEmpty,0.0)); //initial configuration

 auto cmpPaths = [](const ContrPath & left, const ContrPath & right){return (std::get<2>(left) < std::get<2>(right));};
 std::priority_queue<ContrPath,std::vector<ContrPath>,decltype(cmpPaths)> priq(cmpPaths); //output: priority queue

 for(decltype(numContractions) pass = 0; pass < numContractions; ++pass){
  for(const auto & contrPath : inputPaths){
   const auto & parentTensNet = std::get<0>(contrPath); //parental tensor network
   const auto numTensors = parentTensNet.getNumTensors(); //number of r.h.s. tensors in the parental tensor network
   const auto & parentContrSeq = std::get<1>(contrPath); //parental contraction sequence
   for(unsigned int i = 1; i < numTensors; ++i){ //r.h.s. tensors are numbered from 1
    for(unsigned int j=i+1; j <= numTensors; ++j){
     double contrCost = parentTensNet.getContractionCost(i,j); //tensor contraction cost
     auto tensNet = parentTensNet.contractTensorsOut(i,j); //contract tensors i and j
     auto cSeq = parentContrSeq; cSeq.emplace_back(std::pair<unsigned int, unsigned int>(i,j));
     priq.emplace(std::make_tuple(*tensNet,cSeq,contrCost+std::get<2>(contrPath)));
     if(priq.size() > numWalkers) priq.pop();
    }
   }
  }
  inputPaths.clear();
  if(pass == numContractions - 1){
   while(priq.size() > 1) priq.pop();
   contrSeq = std::get<1>(priq.top());
   priq.pop();
  }else{
   while(priq.size() > 0){
    inputPaths.emplace_back(priq.top());
    priq.pop();
   }
  }
 }
 std::cout << "Done: "; //debug
 for(const auto & cPair : contrSeq) std::cout << " {" << std::get<0>(cPair) << "," << std::get<1>(cPair) << "}"; //debug
 std::cout << std::endl; //debug
 return;
}

/** Performs all tensor contractions, thus evaluating the value of the output tensor. **/
template <typename T>
int TensorNetwork<T>::computeOutput(const ContractionSequence & contrSeq)
{
 int error_code = 0; //success
 std::cout << "#MSG(TensorNetwork<T>::computeOutput): Computing ... "; //debug
 //`Finish
 std::cout << "Done" << std::endl; //debug
 return error_code;
}
