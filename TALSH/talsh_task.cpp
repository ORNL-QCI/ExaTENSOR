/** ExaTensor::TAL-SH: Device-unified user-level C++ API implementation.
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

#include <assert.h>

#include "talsh_task.hpp"

namespace talsh{


TensorTask::TensorTask()
{
 int errc = talshTaskClean(&talsh_task_);
 assert(errc == TALSH_SUCCESS);
}


TensorTask::~TensorTask()
{
 if(talshTaskIsEmpty(&talsh_task_) != YEP){
  int stats;
  int errc = talshTaskWait(&talsh_task_,&stats);
  assert(errc == TALSH_SUCCESS);
 }
 int errc = talshTaskDestruct(&talsh_task_);
 assert(errc == TALSH_SUCCESS);
}


talsh_task_t * TensorTask::get_talsh_task_ptr()
{
 return &talsh_task_;
}


} //namespace talsh
