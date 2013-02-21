!#########################################################
program main
!#########################################################

    use ch
    use coef
    use ind
    use logic
    use petsc_ksp_module
    use sc
    use var
    implicit none
#include <finclude/petscsys.h>

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

    open(unit=2,FILE='grid.out')
    rewind 2

    open(unit=10,FILE='ERR.out')
    rewind 10

    call init
    call setUpKSP

    ITIMS=1
    ITIME=1
!....START TIME LOOP
    do ITIM=ITIMS,ITIME
        if (LTIME) then
            TIME=TIME+DT
            do K=1,NK
            do I=1,NI
            do J=1,NJ
                IJK=LK(K)+LI(I)+J
                TO(IJK)=T(IJK)
            end do
            end do
            end do
       end if

       call updateBd

!....START OUTER ITERATIONS
        LSG=10000
        !LSG=1
        do LS=1,LSG
            print *, 'OUTER ITERATION: ', LS
            call calcsc
            if (CONVERGED) then
                call writeVtk
                exit
            end if
        end do
        !call setField
    end do

    ! Free work space

    call cleanUp(A,b,sol)
    call PetscFinalize(ierr)

end program main

!#########################################################
subroutine init
!#########################################################
    
    use bc
    use coef
    use geo
    use ind
    use logic
    use param
    use var
    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscksp.h>
#include <finclude/petscpc.h>

    PetscErrorCode :: ierr

    NI=NXA
    NJ=NYA
    NK=NZA
    NIJ=NI*NJ
    NIJK=NXYZA
    NIM=NI-1
    NJM=NJ-1
    NKM=NK-1
    NICV=NIM-1
    NJCV=NJM-1
    NKCV=NKM-1
    N=NICV*NJCV*NKCV
    NDIRA=NDIR
    NIJCV=NICV*NJCV

    read(2,*)  (ITB(1,I),IK=1,NI*NK)
    read(2,*)  (ITB(2,I),IK=1,NI*NK)
    read(2,*)  (JTB(1,J),JK=1,NJ*NK)
    read(2,*)  (JTB(2,J),JK=1,NJ*NK)
    read(2,*)  (KTB(1,K),IJ=1,NI*NJ)

    read(2,*)  (KTB(2,K),IJ=1,NI*NJ)
    read(2,*)  (LK(K),K=1,NK)
    read(2,*)  (LI(I),I=1,NI)
    read(2,*)  (CTD(I),I=0,NIJK-1)
    read(2,*)  (DTC(I),I=1,NIJK)

    read(2,*)  (IJKD(I),I=1,NDIRA)
    read(2,*)  (IJKPD(I),I=1,NDIRA)
    read(2,*)  (IJKD1(I),I=1,NDIRA)
    read(2,*)  (IJKD2(I),I=1,NDIRA)
    read(2,*)  (IJKD3(I),I=1,NDIRA)

    read(2,*)  (IJKD4(I),I=1,NDIRA)
    read(2,*)  (X(I),I=1,NIJK)
    read(2,*)  (Y(I),I=1,NIJK)
    read(2,*)  (Z(I),I=1,NIJK)
    read(2,*)  (XC(I),I=1,NIJK)

    read(2,*)  (YC(I),I=1,NIJK)
    read(2,*)  (ZC(I),I=1,NIJK)
    read(2,*)  (FX(I), I=1,NIJK)
    read(2,*)  (FY(I), I=1,NIJK)
    read(2,*)  (FZ(I), I=1,NIJK)

    read(2,*)  DX,DY,DZ, VOL
    read(2,*)  (SRDD(I),I=1,NDIRA)

    TIME=1.0d0

    if (LTIME) then
        call setField
        call writeVtk
    end if

    ! Create Matrix
    call MatCreate(PETSC_COMM_WORLD,A,ierr)
    call MatSetType(A,MATSEQAIJ,ierr)
    call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N,ierr)
    call MatSetFromOptions(A,ierr) 
    
    ! Increase performance during matrix assembly due to preallocation
    call MatSeqAIJSetPreallocation(A,i7,PETSC_NULL_INTEGER,ierr)

    ! Create Vector
    call VecCreate(PETSC_COMM_WORLD,sol,ierr)
    call VecSetSizes(sol,PETSC_DECIDE,N,ierr)
    call VecSetFromOptions(sol,ierr)
    call VecDuplicate(sol,b,ierr)

end subroutine init

!#########################################################
subroutine setField
!#########################################################

    use geo
    use ind
    use mms
    use var
    implicit none

    IP=0
    do I=1,NI
        do J=1,NJ
            IJ=LI(I)+J
            T(IJ)=phi(XC(IJ), YC(IJ), 0.0d0, TIME)
        end do
    end do

end subroutine setField

!#########################################################
subroutine writeVtk
!#########################################################

    use ch
    use geo
    use ind
    use var
    implicit none

    character(10) :: TIME_CH
    !write(TIME_CH,'(f6.4)') TIME
    write(TIME_CH,'(I1)') ITIM
    !write(TIME_CH,'(I1)') LS
    VTKFILE='grid_'//trim(TIME_CH)//'.vtk'
    print *, ' *** GENERATING .VTK *** '
    open (unit=4,FILE=VTKFILE)
    write(4,'(A)') '# vtk DataFile Version 3.0'
    write(4,'(A)') 'grid'
    write(4,'(A)') 'ASCII'
    write(4,'(A)') 'DATASET STRUCTURED_GRID'
    write(4,'(A I6 I6 I6)') 'DIMENSIONS', NIM, NJM,NKM
    write(4,'(A I9 A)') 'Points', NIM*NJM*NKM, ' float'
    do K=1,NKM
    do I=1,NIM
    do J=1,NJM
        IJK=LK(K)+LI(I)+J
        write(4,'(E20.10,1X,E20.10,1X,E20.10)'), X(IJK), Y(IJK),Z(IJK)
    end do
    end do
    end do
    write(4,'(A10, I9)') 'CELL_DATA ',(NICV*NJCV*NKCV)
    write(4,'(A15)') 'SCALARS T float'
    write(4,'(A20)') 'LOOKUP_TABLE default'
    do K=2,NKM
    do I=2,NIM
    do J=2,NJM
        IJK=LK(K)+LI(I)+J
        write(4,'(F12.8)') T(IJK)
    end do
    end do
    end do

end subroutine writeVtk

!########################################################
subroutine updateBd
!#########################################################

    use bc
    use geo
    use ind
    use mms
    use sc
    use param
    use var
    implicit none
    real(KIND=PREC) :: XE,YE,ZE,XN,YN,ZN,XT,YT,ZT

!...Dirichlet BC
    do ID=1,NDIRA
        IJKB=IJKD(ID)
        T(IJKB)=phi(XC(IJKB),YC(IJKB),ZC(IJKB),TIME)
    end do

!...Calculate MassFluxes
    do K=2,NKM
    do I=2,NIM-1
    do J=2,NJM
        IJK=LK(K)+LI(I)+J
        F1(IJK)=RHO*DY*DZ*VX
    end do
    end do
    end do

    do K=2,NKM
    do I=2,NIM
    do J=2,NJM-1
        IJK=LK(K)+LI(I)+J
        F2(IJK)=RHO*DX*DZ*VY
    end do
    end do
    end do

    do K=2,NKM-1
    do I=2,NIM
    do J=2,NJM
        IJK=LK(K)+LI(I)+J
        F3(IJK)=RHO*DX*DY*VZ
    end do
    end do
    end do
    

end subroutine updateBd

!#########################################################
subroutine calcSc
!#########################################################

    use bc
    use coef
    use geo
    use grad
    use ind
    use mms
    use petsc_ksp_module
    use sc
    use logic
    use var
    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>

    real*8 :: APT, URF
    PetscLogDouble :: time1, time2

    URF=1.0d0

    call gradfi(T,DTX,DTY,DTZ)

    ! initialize Q and AP
    do K=2,NKM
    do I=2,NIM
    do J=2,NJM
        IJK=LK(K)+LI(I)+J
        Q(IJK)=src(XC(IJK), YC(IJK), ZC(IJK), TIME)*VOL
        AP(IJK)=0.0d0
    end do
    end do
    end do

    ! FLUXES THROUGH EAST FACE
    do K=2,NKM
    do I=2,NIM-1
    do J=2,NJM
        IJK=LK(K)+LI(I)+J
        call fluxsc(IJK,IJK+NJ,IJK-NIJ-1,IJK-1,IJK,F1(IJK),FX(IJK),AW(IJK+NJ),AE(IJK))
    end do
    end do
    end do


    ! FLUXES THROUGH NORTH FACE
    do K=2,NKM
    do I=2,NIM
    do J=2,NJM-1
        IJK=LK(K)+LI(I)+J
        call fluxsc(IJK,IJK+1,IJK-NIJ,IJK,IJK-NJ,F2(IJK),FY(IJK),AS(IJK+1),AN(IJK))
    end do
    end do
    end do

    ! FLUXES THROUGH TOP FACE
    do K=2,NKM-1
    do I=2,NIM
    do J=2,NJM
        IJK=LK(K)+LI(I)+J
        call fluxsc(IJK,IJK+NIJ,IJK-NJ-1,IJK-NJ,IJK,F3(IJK),FZ(IJK),AB(IJK+NIJ),AT(IJK))
    end do
    end do
    end do

!
!.....UNSTEADY TERM CONTRIBUTION
!
    if(LTIME) then
        do K=2,NKM
        do I=2,NIM
        do J=2,NJM
            IJK=LK(K)+LI(I)+J
            APT=RHO*VOL/DT
            Q(IJK)=Q(IJK)+APT*TO(IJK)
            AP(IJK)=AP(IJK)+APT
        end do
        end do
        end do
    end if

    call temp

!
!.....FINAL COEFFICIENT AND SOURCE MATRIX FOR FI-EQUATION
!
    do K=2,NKM
    do I=2,NIM
    do J=2,NJM
        IJK=LK(K)+LI(I)+J
        AP(IJK)=(AP(IJK)-AE(IJK)-AW(IJK)-AN(IJK)-AS(IJK)-AT(IJK)-AB(IJK))/URF
        Q(IJK)=Q(IJK)+(1.0d0-URF)*AP(IJK)*T(IJK)
    end do
    end do
    end do
    
    ! Assemble matrix and Vector

    ! Assembly inner CV matrix coefficients

    do K=3,NKM-1
    do I=3,NIM-1
    do J=3,NJM-1
        IJK=LK(K)+LI(I)+J
        IJKP=DTC(IJK)
        row=IJKP
        !
        col(1)=IJKP-NIJCV
        col(2)=IJKP-NJCV
        col(3)=IJKP-1
        col(4)=IJKP
        col(5)=IJKP+1
        col(6)=IJKP+NJCV
        col(7)=IJKP+NIJCV
        !
        val(1)=AB(IJK)
        val(2)=AW(IJK)
        val(3)=AS(IJK)
        val(4)=AP(IJK)
        val(5)=AN(IJK)
        val(6)=AE(IJK)
        val(7)=AT(IJK)
        valq=Q(IJK)
        !
        call MatSetValues(A,i1,row,i7,col,val,INSERT_VALUES,ierr)
        call VecSetValue(b,row,valq,INSERT_VALUES,ierr)
        !
    end do
    end do
    end do

    ! Assembly matrix coefficients of boundary cv

    do ID=1,NDIRA
        IJK=IJKPD(ID)
        IJKP=DTC(IJK)
        row=IJKP
        col=(/-1,-1,-1,IJKP,-1,-1,-1/)
        !
        val(4)=AP(IJK)
        valq=Q(IJK)
        !
        if (AS(IJK).ne.0) then
            val(3)=AS(IJK)
            col(3)=IJKP-1
        end if
        if (AN(IJK).ne.0) then
            val(5)=AN(IJK)
            col(5)=IJKP+1
        end if
        if (AW(IJK).ne.0) then
            val(2)=AW(IJK)
            col(2)=IJKP-NJCV
        end if
        if (AE(IJK).ne.0) then
            val(6)=AE(IJK)
            col(6)=IJKP+NJCV
        end if
        if (AB(IJK).ne.0) then
            val(1)=AB(IJK)
            col(1)=IJKP-NIJCV
        end if
        if (AT(IJK).ne.0) then
            val(7)=AT(IJK)
            col(7)=IJKP+NIJCV
        end if
        !
        call MatSetValues(A,i1,row,i7,col,val,INSERT_VALUES,ierr)
        call VecSetValue(b,row,valq,INSERT_VALUES,ierr)
        !
    end do


    ! assembly matrix and right hand vector

    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call VecAssemblyBegin(b,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call VecAssemblyEnd(b,ierr)

    call PetscGetCPUTime(time1,ierr)
    call solveSys(A,b,sol,N,LS,tol)
    call PetscGetCPUTime(time2,ierr)

    if (CONVERGED) then
        return
    end if

    tges=time2-time1

    do K=2,NIM
    do I=2,NIM
        do J=2,NJM
            IJK=LK(K)+LI(I)+J
            row=DTC(IJK)
            call VecGetValues(sol,i1,row,valt,ierr)
            T(IJK)=valt
        end do
    end do
    end do

    call calcErr

end subroutine calcSc

!################################################################
subroutine gradfi(FI,DFX,DFY,DFZ)
!################################################################

    use bc
    use geo
    use grad
    use ind
    use param
    use sc
    implicit none

    real(KIND=PREC), intent(in out) :: FI(NIJK), DFX(NIJK), DFY(NIJK), DFZ(NIJK)
    integer :: IJK1, IJK2, IJK3, IJK4

!
!...INITIALIZE FIELDS
!
    do IJK=1,NIJK
        DFX(IJK)=0.0d0
        DFY(IJK)=0.0d0
        DFZ(IJK)=0.0d0
    end do

!
!..CONTRRIBUTION FROM INNER EAST SIDES
!..WERTE AUF DEN CV-FLÄCHEN GEWICHTET MIT DER FLÄCHE
!
    do K=2,NKM
    do I=2,NIM-1
    do J=2,NJM
        IJK=LK(K)+LI(I)+J
        FIE=FI(IJ+NJ)*FX(IJ)+FI(IJ)*(1.0d0-FX(IJ))

        IJK4=IJK
        IJK3=IJK4-1
        IJK2=IJK3-NIJ
        IJK1=IJK4-NIJ
        !
        call normalArea(IJK,IJK,IJK2,IJK3,IJK4,AR,DN,NX,NY,NZ)
        !
        SX=AR*NX
        SY=AR*NY
        SZ=AR*NZ

        DFXE=FIE*SX
        DFYE=FIE*SY
        DFZE=FIE*SZ

        DFX(IJK)=DFX(IJK)+DFXE
        DFY(IJK)=DFY(IJK)+DFYE
        DFZ(IJK)=DFZ(IJK)+DFZE
        DFX(IJK+NJ)=DFX(IJK+NJ)-DFXE
        DFY(IJK+NJ)=DFY(IJK+NJ)-DFYE
        DFZ(IJK+NJ)=DFZ(IJK+NJ)-DFZE
    end do
    end do
    end do

!
!.....CONTRRIBUTION FROM INNER NORTH SIDES
!
    do K=2,NKM
    do I=2,NIM
    do J=2,NJM-1
        IJK=LK(K)+LI(I)+J
        FIN=FI(IJ+1)*FY(IJ)+FI(IJ)*(1.0d0-FY(IJ))

        IJK3=IJK
        IJK4=IJK3-NJ
        IJK2=IJK3-NIJ
        IJK1=IJK4-NIJ
        !
        call normalArea(IJK,IJK,IJK2,IJK3,IJK4,AR,DN,NX,NY,NZ)
        !
        SX=AR*NX
        SY=AR*NY
        SZ=AR*NZ

        DFXN=FIN*SX
        DFYN=FIN*SY
        DFZN=FIN*SZ

        DFX(IJK)=DFX(IJK)+DFXN
        DFY(IJK)=DFY(IJK)+DFYN
        DFZ(IJK)=DFZ(IJK)+DFZN
        DFX(IJK+1)=DFX(IJK+1)-DFXN
        DFY(IJK+1)=DFY(IJK+1)-DFYN
        DFZ(IJK+1)=DFZ(IJK+1)-DFZN
    end do
    end do
    end do
!
!.....CONTRRIBUTION FROM INNER TOP SIDES
!
    do K=2,NKM-1
    do I=2,NIM
    do J=2,NJM
        IJK=LK(K)+LI(I)+J
        FIN=FI(IJK+NIJ)*FZ(IJK)+FI(IJK)*(1.0d0-FZ(IJK))

        IJK4=IJK
        IJK3=IJK4-NJ
        IJK1=IJK4-1
        IJK2=IJK3-1
        !
        call normalArea(IJK,IJK,IJK2,IJK3,IJK4,AR,DN,NX,NY,NZ)
        !
        SX=AR*NX
        SY=AR*NY
        SZ=AR*NZ

        DFXN=FIN*SX
        DFYN=FIN*SY
        DFZN=FIN*SZ

        DFX(IJK)=DFX(IJK)+DFXN
        DFY(IJK)=DFY(IJK)+DFYN
        DFZ(IJK)=DFZ(IJK)+DFZN
        DFX(IJK+1)=DFX(IJK+1)-DFXN
        DFY(IJK+1)=DFY(IJK+1)-DFYN
        DFZ(IJK+1)=DFZ(IJK+1)-DFZN
    end do
    end do
    end do
!
!.....CONTRIBUTION FROM WALL BOUNDARIES
!
    do I=1,NDIRA
        IJKB=IJKD(I)
        IJKP=IJKPD(I)
        !IJK1=IJKD1(I)
        IJK2=IJKD2(I)
        IJK3=IJKD3(I)
        IJK4=IJKD4(I)
        !
        call normalArea(IJKP,IJKB,IJK2,IJK3,IJK4,AR,DN,NX,NY,NZ)
        !
        SX=AR*NX
        SY=AR*NY
        SZ=AR*NZ

        DFX(IJKPD(I))=DFX(IJKPD(I))+FI(IJKD(I))*SX
        DFY(IJKPD(I))=DFY(IJKPD(I))+FI(IJKD(I))*SY
        DFZ(IJKPD(I))=DFZ(IJKPD(I))+FI(IJKD(I))*SZ
    end do

!
!.....CALCULATE GRADIENT COMPONENTS AT CV-CENTERS
!
    do K=2,NKM
    do I=2,NIM
    do J=2,NJM
        IJK=LK(K)+LI(I)+J
        DFX(IJK)=DFX(IJK)/VOL
        DFY(IJK)=DFY(IJK)/VOL
        DFZ(IJK)=DFZ(IJK)/VOL
    end do
    end do
    end do


end subroutine gradfi

!################################################################
    subroutine fluxSc(IJKP,IJKN,IJK2,IJK3,IJK4,FM,FAC,CAP,CAN)
!################################################################
        
    use bc
    use coef
    use flux
    use geo
    use sc
    use param
    use var
    implicit none
    real(KIND=PREC), intent(in) :: FM,FAC
    integer, intent(in) :: IJKP,IJKN,IJK2,IJK3,IJK4
    real(KIND=PREC), intent(out) :: CAN, CAP
    real(KIND=PREC) :: ZERO,SMALL

    ZERO=0.0d0
    SMALL=1.0E-20

    G=1.0d0

    FACP=1.0d0-FAC
    FII=T(IJKN)*FAC+T(IJKP)*FACP
    DFXI=DTX(IJKN)*FAC+DTX(IJKP)*FACP
    DFYI=DTY(IJKN)*FAC+DTY(IJKP)*FACP
    DFZI=DTZ(IJKN)*FAC+DTZ(IJKP)*FACP
    
    !
    !.....SURFACE AND DISTANCE VECTOR COMPONENTS, DIFFUSION COEFF.
    !
    call normalArea(IJKP,IJKN,IJK2,IJK3,IJK4,AR,DN,NX,NY,NZ)
    !
    VSOL=ALPHA*(AR/DN+SMALL)
    !
    !.....EXPLICIT CONVECTIVE AND DIFFUSIVE FLUXES
    !
    FCFIE=FM*FII
    FDFIE=ALPHA*(DFXI*AR*NX+DFYI*AR*NY+DFZI*AR*NZ)
    !
    !.....IMPLICIT CONVECTIVE AND DIFFUSIVE FLUXES
    !
    FCFII=MIN(FM,ZERO)*T(IJKN)+MAX(FM,ZERO)*T(IJKP)
    FDFII=VSOL*(DFXI*DN*NX+DFYI*DN*NY+DFZI*DN*NZ)
    !
    !.....COEFFICIENTS, DEFERRED CORRECTION, SOURCE TERMS
    !
    CAN=-VSOL+MIN(FM,ZERO)
    CAP=-VSOL-MAX(FM,ZERO)
    FFIC=G*(FCFIE-FCFII)
    Q(IJKP)=Q(IJKP)-FFIC+FDFIE-FDFII
    Q(IJKN)=Q(IJKN)+FFIC-FDFIE+FDFII

end subroutine fluxsc

!################################################################
subroutine temp
!################################################################

    use bc
    use coef
    use geo
    use ind
    use sc
    use param
    use var
    implicit none

    real(KIND=PREC) :: COEFC,COEFD,ZERO
    integer :: IJK1, IJK2, IJK3, IJK4
    ZERO=0.0d0

!
!.....DIRICHLET BOUNDARY CONDITION (NO WALL!)
!
    do I=1,NDIRA
        IJKP=IJKPD(I)
        IJKB=IJKD(I)
        !IJK1=IJKD1(I)
        IJK2=IJKD2(I)
        IJK3=IJKD3(I)
        IJK4=IJKD4(I)
        !
        call normalArea(IJKB,IJKP,IJK2,IJK3,IJK4,AR,DN,NX,NY,NZ)
        !
        SX=AR*NX
        SY=AR*NY
        SZ=AR*NZ
        !
        COEFC=RHO*(SX*VX+SY*VY+SZ*VZ)
        COEFD=ALPHA*SRDD(I)
        AP(IJKP)=AP(IJKP)+COEFD-COEFC
        Q(IJKP)=Q(IJKP)+(COEFD-COEFC)*T(IJKB)
    end do

end subroutine temp

!################################################################
subroutine calcErr
!################################################################

    use geo
    use ind
    use logic
    use mms
    use param
    use var
    implicit none

    real(KIND=PREC) :: E,ER
    E=0.0d0
    ER=0.0d0

    do K=2,NKM
    do I=2,NIM
    do J=2,NJM
        IJK=LK(K)+LI(I)+J
        !ER=T(IJ)-phi(XC(IJ),YC(IJ),0.0d0,TIME)
        !E=abs(T(IJ)-phi(XC(IJ),YC(IJ),0.0d0,TIME))
        E=E+abs(T(IJK)-phi(XC(IJK),YC(IJK),ZC(IJK),TIME))
        !ER=max(E,ER)
    end do
    end do
    end do
    
    rewind 10
    write(10, *), E/dble(N), tges, reasonInt, itsInt, LS
    !write(10, *), ER 
    print *, 'ERROR ', E/dble(N), tges, reasonInt, itsInt, LS
    !print *,'ERROR: ', ER

end subroutine calcErr

!################################################################
subroutine linSys(N,A,B,X) 
!################################################################
    
    use param
    implicit none
    integer, intent(in) :: N
    real(KIND=PREC), dimension(N, N), intent(in) :: A
    real(KIND=PREC), dimension(N), intent(in) :: B
    real(KIND=PREC), dimension(N), intent(inout) :: X
    integer :: K, I, J
    real(KIND=PREC), dimension(N) :: X_0
    real(KIND=PREC) :: SIGMA, EPS, R
    
    EPS=1e-6
    R=EPS
    
    open(unit=11,FILE='coef.out')

    do I=1,N
        X(I)=0
        write(11, *), (A(I,J), J=1,N), B(I)
    end do
    K=0

    do 
        if (R.GE.EPS) then
            K=K+1
            do I=1,N
                X_0(I)=X(I)
            end do
            do I=1,N
                SIGMA=0
                do J=1,N
                    if(J.NE.I) then
                            SIGMA=SIGMA+A(I,J)*X_0(J)
                    end if
                end do
                X(I)=(B(I)-SIGMA)/A(I,I)
            end do
            R=0
            do I=1,N
                R=R+abs(X(I)-X_0(I))
            end do
            R=R/N
        else
            exit
        end if
    end do

end subroutine linSys
