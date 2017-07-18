/** C++ adapters for ExaTENSOR: Testing

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

#include <memory>
#include <complex>
#include <iostream>

#include "tensor_expression.hpp"

int test_tensor_expression(){

 //Type aliases:
 using TensData = std::complex<double>;
 using Tensor = exatensor::TensorDenseAdpt<TensData>;
 using TensorLeg = exatensor::TensorLeg;
 using TensorConn = exatensor::TensorConn<TensData>;
 using TensorNetwork = exatensor::TensorNetwork<TensData>;

 //Tensor dimension extents:
 unsigned int rank = 4;
 auto dims = new std::size_t[rank];
 for(unsigned int i=0; i<rank; ++i) dims[i]=8;

 //Tensor volume:
 std::size_t vol = 1;
 for(unsigned int i=0; i<rank; ++i) vol*=dims[i];

 //Tensor body:
 std::shared_ptr<TensData> body(new TensData[vol],[](TensData * p){delete[] p;});

 //Tensor 0:
 Tensor tensor0(rank,dims,body);

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
 TensorNetwork tensnet;
 //Legs of tensor 0:
 std::vector<TensorLeg> legs;
 legs.push_back(TensorLeg(1,3));
 legs.push_back(TensorLeg(1,0));
 legs.push_back(TensorLeg(2,1));
 legs.push_back(TensorLeg(2,2));
 tensnet.appendTensor(tensor0,legs);
 //Legs of tensor 1:
 legs.clear();
 legs.push_back(TensorLeg(0,1));
 legs.push_back(TensorLeg(2,3));
 legs.push_back(TensorLeg(2,0));
 legs.push_back(TensorLeg(0,0));
 tensnet.appendTensor(tensor1,legs);
 //Legs of tensor 2:
 legs.clear();
 legs.push_back(TensorLeg(1,2));
 legs.push_back(TensorLeg(0,2));
 legs.push_back(TensorLeg(0,3));
 legs.push_back(TensorLeg(1,1));
 tensnet.appendTensor(tensor2,legs);
 legs.clear();
 //Print the tensor network:
 tensnet.printIt();

 //Free dimension extents:
 delete[] dims;

 //Done:
 return 0;
}

int main(int argc, char ** argv){
 return test_tensor_expression();
}
