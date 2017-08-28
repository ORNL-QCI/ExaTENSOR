/** C++ adapters for ExaTENSOR: Testing

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/08/03

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

#include <memory>
#include <complex>
#include <iostream>

#include "tensor_network.hpp"

int test_tensor_expression(){

 //Parameters:
 const std::size_t TENS_DIM_EXT=8;

 //Type aliases:
 using TensDataType = std::complex<double>;
 using Tensor = exatensor::TensorDenseAdpt<TensDataType>;
 using TensorLeg = exatensor::TensorLeg;
 using TensorConn = exatensor::TensorConn<TensDataType>;
 using TensorNetwork = exatensor::TensorNetwork<TensDataType>;
 using ContractionSequence = exatensor::ContractionSequence;
 using PairUnsignedInt = std::pair<unsigned int, unsigned int>;

 std::size_t dims[32];

 //Tensor dimension extents:
 unsigned int rank = 4; for(unsigned int i = 0; i < rank; ++i) dims[i]=TENS_DIM_EXT;
 //Tensor volume:
 std::size_t vol = 1; for(unsigned int i = 0; i < rank; ++i) vol*=dims[i];
 //Persistent tensor body:
 std::shared_ptr<TensDataType> body(new TensDataType[vol], [](TensDataType * p){delete[] p;});
 //Tensor 0:
 Tensor tensor0(rank,dims);

 //Tensor 0 accessors:
 std::cout << "Rank = " << tensor0.getRank() << std::endl;
 const std::size_t * dim = tensor0.getDimExtents();
 std::cout << "Dim extents: ";
 for(unsigned int i=0; i<tensor0.getRank(); ++i) std::cout << dim[i] << " ";
 std::cout << std::endl;
 std::cout << "Volume = " << tensor0.getVolume() << ": Size = " << tensor0.getSize() << std::endl;

 //Print tensor 0:
 tensor0.printIt();

 //Copy assignment:
 Tensor tensor1 = tensor0;

 //Print tensor 1:
 tensor1.printIt();

 //Copy constructor:
 Tensor tensor2(tensor1);

 //Print tensor 2:
 tensor2.printIt();

 //Create an empty tensor network:
 TensorNetwork tensnet0;
 //Append tensor 0:
 std::vector<TensorLeg> legs;
 legs.push_back(TensorLeg(1,3));
 legs.push_back(TensorLeg(1,0));
 legs.push_back(TensorLeg(2,1));
 legs.push_back(TensorLeg(2,2));
 tensnet0.appendTensor(tensor0,legs);
 //Append tensor 1:
 legs.clear();
 legs.push_back(TensorLeg(0,1));
 legs.push_back(TensorLeg(2,3));
 legs.push_back(TensorLeg(2,0));
 legs.push_back(TensorLeg(0,0));
 tensnet0.appendTensor(tensor1,legs);
 //Append tensor 2:
 legs.clear();
 legs.push_back(TensorLeg(1,2));
 legs.push_back(TensorLeg(0,2));
 legs.push_back(TensorLeg(0,3));
 legs.push_back(TensorLeg(1,1));
 tensnet0.appendTensor(tensor2,legs);
 legs.clear();
 //Print the tensor network:
 std::cout << std::endl;
 tensnet0.printIt();
 //Create another tensor and append it to the tensor network:
 Tensor tensor3(tensor2);
 std::vector<std::pair<unsigned int, unsigned int>> legPairs;
 legPairs.push_back(std::pair<unsigned int, unsigned int>(0,3));
 legPairs.push_back(std::pair<unsigned int, unsigned int>(3,0));
 tensnet0.appendTensor(tensor3,legPairs);
 legPairs.clear();
 //Print the tensor network:
 std::cout << std::endl;
 tensnet0.printIt();

/*
 //Contract tensors 1 and 2:
 double ai;
 std::cout << std::endl << "Contraction cost = " << tensnet0.getContractionCost(1,2,&ai);
 std::cout << "; Arithmetic intensity = " << ai;
 tensnet0.contractTensors(1,2);
 //Print the new tensor network:
 std::cout << std::endl;
 tensnet0.printIt();
*/

 TensorNetwork tensnet1;

 rank = 4; for(unsigned int i = 0; i < rank; ++i) dims[i]=TENS_DIM_EXT;
 Tensor tens1(rank,dims); vol = tens1.getVolume();
 std::shared_ptr<TensDataType> body1(new TensDataType[vol], [](TensDataType * p){delete[] p;});
 tens1.setBody(body1);

 rank = 4; for(unsigned int i = 0; i < rank; ++i) dims[i]=TENS_DIM_EXT;
 Tensor tens2(rank,dims); vol = tens2.getVolume();
 std::shared_ptr<TensDataType> body2(new TensDataType[vol], [](TensDataType * p){delete[] p;});
 tens2.setBody(body2);

 rank = 4; for(unsigned int i = 0; i < rank; ++i) dims[i]=TENS_DIM_EXT;
 Tensor tens3(rank,dims); vol = tens3.getVolume();
 std::shared_ptr<TensDataType> body3(new TensDataType[vol], [](TensDataType * p){delete[] p;});
 tens3.setBody(body3);

 rank = 5; for(unsigned int i = 0; i < rank; ++i) dims[i]=TENS_DIM_EXT;
 Tensor tens4(rank,dims); vol = tens4.getVolume();
 std::shared_ptr<TensDataType> body4(new TensDataType[vol], [](TensDataType * p){delete[] p;});
 tens4.setBody(body4);

 rank = 4; for(unsigned int i = 0; i < rank; ++i) dims[i]=TENS_DIM_EXT;
 Tensor tens5(rank,dims); vol = tens5.getVolume();
 std::shared_ptr<TensDataType> body5(new TensDataType[vol], [](TensDataType * p){delete[] p;});
 tens5.setBody(body5);

 rank = 5; for(unsigned int i = 0; i < rank; ++i) dims[i]=TENS_DIM_EXT;
 Tensor tens6(rank,dims); vol = tens6.getVolume();
 std::shared_ptr<TensDataType> body6(new TensDataType[vol], [](TensDataType * p){delete[] p;});
 tens6.setBody(body6);

 std::vector<PairUnsignedInt> connections;
 tensnet1.appendTensor(tens1,connections);
 connections.push_back(PairUnsignedInt(2,0)); connections.push_back(PairUnsignedInt(3,1));
 tensnet1.appendTensor(tens2,connections);
 connections.clear();
 connections.push_back(PairUnsignedInt(1,0)); connections.push_back(PairUnsignedInt(2,1));
 tensnet1.appendTensor(tens3,connections);
 connections.clear();
 connections.push_back(PairUnsignedInt(0,0)); connections.push_back(PairUnsignedInt(2,1)); connections.push_back(PairUnsignedInt(3,2));
 tensnet1.appendTensor(tens4,connections);
 connections.clear();
 connections.push_back(PairUnsignedInt(2,0)); connections.push_back(PairUnsignedInt(0,1));
 tensnet1.appendTensor(tens5,connections);
 connections.clear();
 connections.push_back(PairUnsignedInt(0,0)); connections.push_back(PairUnsignedInt(1,1));
 tensnet1.appendTensor(tens6,connections);
 connections.clear();

 const Tensor & tens0 = tensnet1.getTensor(0);
 vol = tens0.getVolume();
 std::shared_ptr<TensDataType> body0(new TensDataType[vol], [](TensDataType * p){delete[] p;});
 tensnet1.setOutputBody(body0);

 ContractionSequence contrSeq;
 int error_code = tensnet1.evaluate(contrSeq);

 //Done:
 return 0;
}

int main(int argc, char ** argv){
 return test_tensor_expression();
}
