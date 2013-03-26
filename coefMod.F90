module coef

    use param
    implicit none
#include <finclude/petscvecdef.h>
#include <finclude/petscmatdef.h>


    real(KIND=PREC) :: Q(NXYZA), AP(NXYZA), AS(NXYZA), AN(NXYZA),&
                       AW(NXYZA), AE(NXYZA), AB(NXYZA), AT(NXYZA),&
                       AF(NFACEAL)
                       
    Mat :: Amat
    Vec :: solvec,bvec
    PetscInt, Parameter :: i1=1,i2=2,i5=5,i7=7
    PetscInt :: row, col(7),col1
    PetscScalar :: val(7),val1,valq,valt,tol,minres

end module coef
