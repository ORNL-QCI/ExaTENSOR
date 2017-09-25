/** C++ adapters for ExaTENSOR: External tensor definition mechanism

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/09/25

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

#ifndef _EXA_TENSOR_DEFINE_H
#define _EXA_TENSOR_DEFINE_H

#include "talsh.h"

namespace exatensor {

using TensorShape = talsh_tens_shape_t;

class TensorDenseDefiner{

public:

 /** Destructor **/
 virtual ~TensorDenseDefiner();

 /** User-provided tensor definition class (eturns an error code). **/
 virtual int defineTensBody(void * bodyPtr,                                     //pointer to the tensor body
                            const int dataKind,                                 //tensor element data kind: {R4,R8,C4,C8}
                            const TensorShape * tensShape,                      //tensor shape
                            const long long int * tensSignature = nullptr) = 0; //tensor signature (optional)

}; //end class TensorDenseDefiner

} //end namespace exatensor

#endif //_EXA_TENSOR_DEFINE_H
