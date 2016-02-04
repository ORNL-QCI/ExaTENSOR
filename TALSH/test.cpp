/** TALSH::C/C++ API testing. **/

#include <stdio.h>
#include <stdlib.h>

#include "talsh.h"

#ifdef __cplusplus
extern "C"{
#endif
void test_talsh_c(int * ierr);
#ifdef __cplusplus
}
#endif

void test_talsh_c(int * ierr)
{
 int host_arg_max,errc;
 size_t host_buf_size;
 tensBlck_t *t0,*t1,*t2; //tensor blocks (pointers)
 int r0=4,r1=4,r2=4; //tensor block ranks
 int dims0[]={40,40,40,40}; //tensor block dimensions
 int dims1[]={40,40,40,40};
 int dims2[]={40,40,40,40};
 int cptrn0[MAX_TENSOR_RANK*2]; //tensor contraction pattern

 *ierr=0;
//Initialize Host/GPU argument buffers and NV-TAL:
 host_buf_size=1000000000;
 printf(" Initializing NV-TAL ... ");
 errc=arg_buf_allocate(&host_buf_size,&host_arg_max,0,0);
 printf(" Status %d: Host argument buffer size = %lu; Host arg max = %d\n",errc,host_buf_size,host_arg_max);
 if(errc){*ierr=1; return;}
//Create tensor blocks:
 //Tensor block 0:
 printf(" Creating tensor block 0 ... ");
 errc=tensBlck_create(&t0);
 printf(" Status %d\n",errc);
 if(errc){*ierr=1; return;}
 printf(" Constructing shape of tensor block 0 ... ");
 errc=tensBlck_construct(t0,YEP,r0,dims0);
 printf(" Status %d: Tensor block volume = %lu\n",errc,tensBlck_volume(t0));
 if(errc){*ierr=1; return;}
 printf(" Attaching body to tensor block 0 ... ");
 errc=tensBlck_attach_body(t0,R8);
 printf(" Status %d\n",errc);
 if(errc){*ierr=1; return;}
 //Tensor block 1:
 printf(" Creating tensor block 1 ... ");
 errc=tensBlck_create(&t1);
 printf(" Status %d\n",errc);
 if(errc){*ierr=1; return;}
 printf(" Constructing shape of tensor block 1 ... ");
 errc=tensBlck_construct(t1,YEP,r1,dims1);
 printf(" Status %d: Tensor block volume = %lu\n",errc,tensBlck_volume(t1));
 if(errc){*ierr=1; return;}
 printf(" Attaching body to tensor block 1 ... ");
 errc=tensBlck_attach_body(t1,R8);
 printf(" Status %d\n",errc);
 if(errc){*ierr=1; return;}
 //Tensor block 2:
 printf(" Creating tensor block 2 ... ");
 errc=tensBlck_create(&t2);
 printf(" Status %d\n",errc);
 if(errc){*ierr=1; return;}
 printf(" Constructing shape of tensor block 2 ... ");
 errc=tensBlck_construct(t2,YEP,r2,dims2);
 printf(" Status %d: Tensor block volume = %lu\n",errc,tensBlck_volume(t2));
 if(errc){*ierr=1; return;}
 printf(" Attaching body to tensor block 2 ... ");
 errc=tensBlck_attach_body(t2,R8);
 printf(" Status %d\n",errc);
 if(errc){*ierr=1; return;}




//Destroy tensor blocks:
 //Tensor block 2:
 printf(" Destroying tensor block 2 ... ");
 errc=tensBlck_destroy(t2);
 printf(" Status %d\n",errc);
 if(errc){*ierr=1; return;}
 //Tensor block 1:
 printf(" Destroying tensor block 1 ... ");
 errc=tensBlck_destroy(t1);
 printf(" Status %d\n",errc);
 if(errc){*ierr=1; return;}
 //Tensor block 0:
 printf(" Destroying tensor block 0 ... ");
 errc=tensBlck_destroy(t0);
 printf(" Status %d\n",errc);
 if(errc){*ierr=1; return;}
//Free Host/GPU argument buffers and shutdown NV-TAL:
 printf(" Shutting down NV-TAL ... ");
 errc=arg_buf_deallocate(0,0);
 printf(" Status: %d\n",errc);
 if(errc){*ierr=1; return;}
 return;
}
