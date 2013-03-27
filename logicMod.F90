module logic

    implicit none
    logical, parameter :: LTIME=.false.
    logical :: CONVERGED=.false.
    logical :: STARTED,FOUND
    integer :: pcInd, kspInd
    integer :: reasonInt, itsInt
    integer,parameter :: PROCOFFSET=10,BLOCKOFFSET=20
    real*8 :: tges

end module logic
