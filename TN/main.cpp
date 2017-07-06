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

 using TensData = std::complex<double>;
 using Tensor = exatensor::TensorAdpt<TensData>;

 unsigned int rank = 4;
 auto dims = new std::size_t[rank];
 for(unsigned int i=0; i<rank; ++i) dims[i]=8;

 std::size_t vol = 1;
 for(unsigned int i=0; i<rank; ++i) vol*=dims[i];

 std::shared_ptr<TensData> body(new TensData[vol],[](TensData * p){delete[] p;});

 Tensor tensor(rank,dims,body);

 std::cout << "Rank = " << tensor.getRank() << std::endl;

 const std::size_t * dim = tensor.getDimExtents();

 std::cout << "Dim extents: ";
 for(unsigned int i=0; i<tensor.getRank(); ++i) std::cout << dim[i] << " ";
 std::cout << std::endl;

 std::cout << "Volume = " << tensor.getVolume() << ": Size = " << tensor.getSize() << std::endl;

 delete[] dims;

 return 0;
}

int main(int argc, char ** argv){
 return test_tensor_expression();
}
