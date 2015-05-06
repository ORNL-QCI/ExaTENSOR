/** MPI Fortran wrappers for Fortran-unfriendly MPI functions. **/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

//Prototypes:
#ifdef __cplusplus
extern "C" {
#endif
void MPI_Get_Displacement(const void * location, MPI_Aint * disp, int * ierr);
#ifdef __cplusplus
}
#endif

//Function definitions:
void MPI_Get_Displacement(const void * location, MPI_Aint * disp, int * ierr){
*disp=0;
*ierr=MPI_Get_address(location,disp);
//`Check for signed overflow for Fortran
return;
}
