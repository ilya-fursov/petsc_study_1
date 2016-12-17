
static char help[] = "Examine DMDA and numberings\n";


#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>


//#undef __FUNCT__
//#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscInitialize(&argc,&argv,(char*)0,help);
  
  PetscErrorCode ierr;
  int rank, size, N;
  DM da;
  
  ierr  = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
  ierr  = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);
  ierr  = PetscOptionsGetInt(NULL,NULL,"-n",&N,NULL); CHKERRQ(ierr);

  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,N,1,1,NULL,&da); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"N = %d\n", N); CHKERRQ(ierr);

 
  ierr = DMDestroy(&da); CHKERRQ(ierr);
  ierr = PetscFinalize();
  
  return 0;
  //PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------- */







