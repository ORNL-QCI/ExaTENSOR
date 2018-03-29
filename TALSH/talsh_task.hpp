/** ExaTensor::TAL-SH: C++ TAL-SH task
REVISION: 2018/03/29

Copyright (C) 2014-2017 Dmitry I. Lyakh (Liakh)
Copyright (C) 2014-2017 Oak Ridge National Laboratory (UT-Battelle)

This file is part of ExaTensor.

ExaTensor is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ExaTensor is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with ExaTensor. If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
**/

#ifndef _TALSH_TASK_HPP
#define _TALSH_TASK_HPP

#include "talsh.h" //TAL-SH C header

namespace talsh{

/** Tensor task. **/
class TensorTask{

public:

 TensorTask();

 TensorTask(const TensorTask & task_handle) = delete;

 TensorTask & operator=(const TensorTask & task_handle) = delete;

 ~TensorTask();

private:

 talsh_task_t * get_talsh_task_ptr();

 talsh_task_t talsh_task_;

 friend class Tensor;

};

} //namespace talsh

#endif //_TALSH_TASK_HPP
