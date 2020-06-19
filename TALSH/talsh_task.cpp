/** ExaTensor::TAL-SH: Device-unified user-level C++ API implementation.
REVISION: 2020/06/19

Copyright (C) 2014-2020 Dmitry I. Lyakh (Liakh)
Copyright (C) 2014-2020 Oak Ridge National Laboratory (UT-Battelle)

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

#include "talsh_task.hpp"
#include "talshxx.hpp"

#include <iostream>
#include <cassert>

namespace talsh{

TensorTask::TensorTask():
 num_tensors_(0)
{
 int errc = talshTaskClean(&talsh_task_);
 assert(errc == TALSH_SUCCESS);
}


TensorTask::~TensorTask()
{
 this->wait();
 int errc = talshTaskDestruct(&talsh_task_);
 assert(errc == TALSH_SUCCESS);
}


bool TensorTask::isEmpty()
{
 return (talshTaskIsEmpty(&talsh_task_) == YEP);
}


void TensorTask::clean()
{
 int status = talshTaskStatus(&talsh_task_);
 if(status != TALSH_TASK_ERROR && status != TALSH_TASK_EMPTY && status != TALSH_TASK_COMPLETED) this->wait();
 for(int i = 0; i < num_tensors_; ++i){
  Tensor * tensor = used_tensors_[i];
  if(tensor->getWriteTask() == this) tensor->resetWriteTask();
 }
 num_tensors_ = 0;
 int errc = talshTaskDestruct(&talsh_task_);
 assert(errc == TALSH_SUCCESS);
 return;
}


bool TensorTask::wait()
{
 int stats = TALSH_TASK_COMPLETED;
 if(talshTaskIsEmpty(&talsh_task_) != YEP){
  int errc;
  int completed = talshTaskComplete(&talsh_task_,&stats,&errc);
  if(errc != TALSH_SUCCESS) std::cout << "#ERROR(TAL-SH:TensorTask.wait): Task completion check failed: Error " << errc << std::endl; //debug
  assert(errc == TALSH_SUCCESS);
  if(completed != YEP){
   errc = talshTaskWait(&talsh_task_,&stats);
   if(errc != TALSH_SUCCESS) std::cout << "#ERROR(TAL-SH:TensorTask.wait): Task completion wait failed: Error " << errc << std::endl; //debug
   assert(errc == TALSH_SUCCESS);
  }
  if(stats != TALSH_TASK_COMPLETED){ //debug
   std::cout << "#ERROR(TAL-SH:TensorTask.wait): Task completed with error: Status " << stats << std::endl;
   talshTaskPrint(&talsh_task_);
   assert(stats == TALSH_TASK_COMPLETED);
  }
 }
 this->clean();
 return (stats == TALSH_TASK_COMPLETED);
}


bool TensorTask::test(int * status)
{
 bool res = true;
 if(talshTaskIsEmpty(&talsh_task_) != YEP){
  int errc;
  res = (talshTaskComplete(&talsh_task_,status,&errc) == YEP);
  if(errc != TALSH_SUCCESS) std::cout << "#ERROR(TAL-SH:TensorTask.test): Task completion check failed: Error " << errc << std::endl; //debug
  assert(errc == TALSH_SUCCESS);
 }else{ //empty task: Completed = TRUE
  *status = TALSH_TASK_EMPTY;
 }
 if(res) this->clean();
 return res;
}


int TensorTask::getExecutionDevice(int * device_kind)
{
 return talshTaskDevId(&talsh_task_,device_kind);
}


unsigned int TensorTask::getNumTensorArguments() const
{
 return num_tensors_;
}


const Tensor * TensorTask::getTensorArgument(unsigned int arg_num) const
{
 if(arg_num < num_tensors_) return used_tensors_[arg_num];
 return nullptr;
}


const talshTensArg_t * TensorTask::getTensorArgumentImages(int * num_arguments)
{
 return talshTaskTensArgs(&talsh_task_,num_arguments);
}


int TensorTask::getTensorArgumentCoherence()
{
 return talshTaskArgCoherence(&talsh_task_);
}


talsh_task_t * TensorTask::getTalshTaskPtr()
{
 return &talsh_task_;
}

} //namespace talsh
