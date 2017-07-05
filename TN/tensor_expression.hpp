/** C++ adapters for ExaTENSOR **/

#ifndef _TENSOR_EXPRESSION_H
#define _TENSOR_EXPRESSION_H

#include <memory>

namespace exatensor {

template <typename T>
class TensorAdapter{

private:

 unsigned int rank;              //tensor rank (number of dimensions)
 unsigned long long  *dimExtent; //tensor dimension extents
 std::shared_ptr<T> pBody;       //pointer to the tensor body (tensor elements)

public:

 TensorAdapter(int Rank, unsigned long long  DimExtent[], std::shared_ptr<T> DataPtr):
 rank(Rank){}

 ~TensorAdapter(){}

};

} //end namespace talsh

#endif //_TENSOR_EXPRESSION_H
