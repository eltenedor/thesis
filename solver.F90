!#########################################################
program main
!#########################################################

    use ch
    use ind
    use logic
    use sc
    use var
    implicit none


    open(unit=2,FILE='grid.out')
    rewind 2

    open(unit=10,FILE='ERR.out')
    rewind 10

    call init

    ITIMS=1
    ITIME=1
    do ITIM=ITIMS,ITIME
        if (LTIME) then
            TIME=TIME+DT
           do I=1,NI
                do J=1,NJ
                   IJ=LI(I)+J
                   TO(IJ)=T(IJ)
               end do
           end do
       end if

        call updateBd
        call calcSc
        !call setField
        call writeVtk

    end do

end program main


!#########################################################
subroutine init
!#########################################################
    
    use bc
    use geo
    use ind
    use logic
    use param
    use var
    implicit none

    NI=NXA
    NJ=NYA
    NIJ=NXYA
    NIM=NI-1
    NJM=NJ-1
    NICV=NIM-1
    NJCV=NJM-1
    N=NICV*NJCV
    NWALI=NWA

    read(2,*) (ITB(1,I),I=1,NI), (ITB(2,I),I=1,NI),&
            (JTB(1,J),J=1,NJ), (JTB(2,J),J=1,NJ),& 
            (LI(I),I=1,NI),(CTD(IJP),IJP=0,NIJ-1),(DTC(I),I=1,NIJ),&
            (IJW(I),I=1,NWALI), (IJPW(I),I=1,NWALI),(IJW1(I),I=1,NWALI),&
            (IJW2(I),I=1,NWALI)
    read(2,*) (X(I),I=1,NIJ), (Y(I),I=1,NIJ), (XC(I),I=1,NIJ),&
            (YC(I),I=1,NIJ), (FX(I), I=1,NIJ), (FY(I), I=1,NIJ),DX,DY, VOL,&
            (SRDW(I),I=1,NWALI),(XTW(I),I=1,NWALI),(YTW(I),I=1,NWALI)

    TIME=1.0d0

    if (LTIME) then
        call setField
        call writeVtk
    end if
        
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
    VTKFILE='grid_'//trim(TIME_CH)//'.vtk'
    print *, ' *** GENERATING .VTK *** '
    open (unit=4,FILE=VTKFILE)
    write(4,'(A)') '# vtk DataFile Version 3.0'
    write(4,'(A)') 'grid'
    write(4,'(A)') 'ASCII'
    write(4,'(A)') 'DATASET STRUCTURED_GRID'
    write(4,'(A I6 I6 I6)') 'DIMENSIONS', NIM, NJM,1
    write(4,'(A I9 A)') 'Points', NIM*NJM, ' float'
    DO I=1,NIM
        DO J=1,NJM
            IJ=LI(I)+J
            write(4,'(E20.10,1X,E20.10,1X,E20.10)'), X(IJ), Y(IJ),0.0
        END DO
    END DO
    write(4,'(A10, I9)') 'CELL_DATA ', (NICV*NJCV)
    write(4,'(A15)') 'SCALARS T float'
    write(4,'(A20)') 'LOOKUP_TABLE default'
    do I=2,NIM
        do J=2,NJM
            IJ=LI(I)+J
            write(4,'(F12.8)') T(IJ)
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
    use var
    implicit none

!...Dirichlet BC
    do IW=1,NWALI
        IJB=IJW(IW)
        T(IJB)=phi(XC(IJB),YC(IJB),0.0d0,TIME)
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
    use sc
    use logic
    use var
    implicit none
    
#include <finclude/petscsys.h>
#include <finclude/petscksp.h>
#include <finclude/petscpc.h>

    real*8 :: APT, URF
    KSP :: ksp
    PC :: pc
    PetscReal :: tol
    PetscErrorCode :: ierr

    URF=1.0

    !call gradfi(T,DTX,DTY)

    do I=2,NIM
        do IJ=LI(I)+2,LI(I)+NIM
            Q(IJ)=src(XC(IJ), YC(IJ), 0.0d0, TIME)*VOL
            !print *, Q(IJ)
            AP(IJ)=0.0d0
        end do
    end do

    do I=2,NIM-1
        do IJ=LI(I)+2,LI(I)+NJM
        call fluxsc(IJ,IJ+NJ,IJ,IJ-1,RHO*VX*DY,FX(IJ),AW(IJ+NJ),AE(IJ))
        !print *, IJ, AW(IJ)
        end do
    end do

    do I=2,NIM
        do IJ=LI(I)+2,LI(I)+NJM-1
        call fluxsc(IJ,IJ+1,IJ-NJ,IJ,RHO*VY*DX,FY(IJ),AS(IJ+1),AN(IJ))
        end do
    end do

!
!.....UNSTEADY TERM CONTRIBUTION
!
    if(LTIME) then
        DO I=2,NIM
        DO IJ=LI(I)+2,LI(I)+NJM
          APT=RHO*VOL/DT
          Q(IJ)=Q(IJ)+APT*TO(IJ)
          AP(IJ)=AP(IJ)+APT
        END DO
        END DO
    end if

    call temp

!
!.....FINAL COEFFICIENT AND SOURCE MATRIX FOR FI-EQUATION
!
      DO I=2,NIM
      DO IJ=LI(I)+2,LI(I)+NJM
        AP(IJ)=(AP(IJ)-AE(IJ)-AW(IJ)-AN(IJ)-AS(IJ))/URF
        Q(IJ)=Q(IJ)+(1.-URF)*AP(IJ)*T(IJ)
      END DO
      END DO
    
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

    ! Create Matrix

    call MatCreate(PETSC_COMM_WORLD,A,ierr)
    call MatSetType(A,MATSEQAIJ,ierr)
    call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N,ierr)
    call MatSetFromOptions(A,ierr) 
    
    ! Increase performance during matrix assembly due to preallocation
    call MatSeqAIJSetPreallocation(A,i5,PETSC_NULL_INTEGER,ierr)

    ! Not necessary anymore, since preallocation is used
    !call MatSetUp(A,ierr)

    ! Create Vector

    call VecCreate(PETSC_COMM_WORLD,sol,ierr)
    call VecSetSizes(sol,PETSC_DECIDE,N,ierr)
    call VecSetFromOptions(sol,ierr)
    call VecDuplicate(sol,b,ierr)

    ! Assemble matrix and Vector

    ! Assembly inner CV matrix coefficients

    do I=3,NIM-1
        do IJ=LI(I)+3,LI(I)+NJM-1
            IJP=DTC(IJ)
            row=IJP
            !
            col(1)=IJP-NJCV
            col(2)=IJP-1
            col(3)=IJP
            col(4)=IJP+1
            col(5)=IJP+NJCV
            !
            val(1)=AW(IJ)
            val(2)=AS(IJ)
            val(3)=AP(IJ)
            val(4)=AN(IJ)
            val(5)=AE(IJ)
            valq=Q(IJ)
            !
            call MatSetValues(A,i1,row,i5,col,val,INSERT_VALUES,ierr)
            call VecSetValue(b,row,valq,INSERT_VALUES,ierr)
            !
        end do
    end do

    ! Assembly matrix coefficients of boundary cv

    do IW=1,NWALI
        IJ=IJPW(IW)
        IJP=DTC(IJ)
        row=IJP
        col=(/-1,-1,IJP,-1,-1/)
        !
        val(3)=AP(IJ)
        valq=Q(IJ)
        !
        if (AS(IJ).ne.0) then
            val(2)=AS(IJ)
            col(2)=IJP-1
        end if
        if (AN(IJ).ne.0) then
            val(4)=AN(IJ)
            col(4)=IJP+1
        end if
        if (AW(IJ).ne.0) then
            val(1)=AW(IJ)
            col(1)=IJP-NJCV
        end if
        if (AE(IJ).ne.0) then
            val(5)=AE(IJ)
            col(5)=IJP+NJCV
        end if
        !
        call MatSetValues(A,i1,row,i5,col,val,INSERT_VALUES,ierr)
        call VecSetValue(b,row,valq,INSERT_VALUES,ierr)
        !
    end do

    ! assembly matrix and right hand vector

    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call VecAssemblyBegin(b,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call VecAssemblyEnd(b,ierr)

    if(N.le.10) call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)

    ! create linear solver context

    call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)

    ! Set operators

    call KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN,ierr)

    ! Set linear solver defaults

    call KSPGetPC(ksp,pc,ierr)
    call PCSetType(pc,PCJACOBI,ierr)
    tol=1.0e-10
    call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_DOUBLE_PRECISION, &
         & PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr)

    ! Set runtime options

    call KSPSetFromOptions(ksp,ierr)

    ! Solve the linear system

    call KSPSolve(ksp,b,sol,ierr)

    ! View solver info

    call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
    if(N.le.10) then
        call VecView(sol,PETSC_VIEWER_STDOUT_WORLD,ierr)
        call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr)
    end if

    do I=2,NIM
        do J=2,NJM
            IJ=LI(I)+J
            col=DTC(IJ)
            call VecGetValues(sol,i1,col,valt,ierr)
            T(IJ)=valt
        end do
    end do

    ! Free work space

    call VecDestroy(sol,ierr)
    call VecDestroy(b,ierr)
    call MatDestroy(A,ierr)
    call KSPDestroy(ksp,ierr)
    call PetscFinalize(ierr)

    call calcErr

end subroutine calcSc

!################################################################
subroutine gradfi(FI,DFX,DFY)
!################################################################

    use bc
    use geo
    use grad
    use ind
    use sc
    implicit none

    real(KIND=PREC), intent(in out) :: FI(NIJ), DFX(NIJ), DFY(NIJ)

    do IJ=1,NIJ
        DFX(IJ)=0
        DFX(IJ)=0
    end do

!
!..CONTRRIBUTION FROM INNER EAST SIDES
!
    do I=2,NIM-1
        do IJ=LI(I)+2,LI(I)*NJM
            FIE=FI(IJ+NJ)*FX(IJ)+FI(IJ)*(1.0d0-FX(IJ))
            SX=Y(IJ)-Y(IJ-1)
            SY=X(IJ-1)-X(IJ)
            DFXE=FIE*SX
            DFYE=FIE*SY

            DFX(IJ)=DFX(IJ)+DFXE
            DFY(IJ)=DFY(IJ)+DFYE
            DFX(IJ+NJ)=DFX(IJ+NJ)-DFXE
            DFY(IJ+NJ)=DFY(IJ+NJ)-DFYE
        end do
    end do
!
!.....CONTRRIBUTION FROM INNER NORTH SIDES
!
    do I=2,NIM
        do IJ=LI(I)+2,LI(I)+NJM-1
            FIN=FI(IJ+1)*FY(IJ)+FI(IJ)*(1.0d0-FY(IJ))
            SX=Y(IJ-NJ)-Y(IJ)
            SY=X(IJ)-X(IJ-NJ)
            DFXN=FIN*SX
            DFYN=FIN*SY

            DFX(IJ)=DFX(IJ)+DFXN
            DFY(IJ)=DFY(IJ)+DFYN
            DFX(IJ+1)=DFX(IJ+1)-DFXN
            DFY(IJ+1)=DFY(IJ+1)-DFYN
        end do
    end do
!
!.....CONTRIBUTION FROM WALL BOUNDARIES
!
    do I=1,1+NWALI
        SX=Y(IJW1(I))-Y(IJW2(I))
        SY=X(IJW2(I))-X(IJW1(I))
        DFX(IJPW(I))=DFX(IJPW(I))+FI(IJW(I))*SX
        DFY(IJPW(I))=DFY(IJPW(I))+FI(IJW(I))*SY
    end do

!
!.....CALCULATE GRADIENT COMPONENTS AT CV-CENTERS
!
    do I=2,NIM
        do IJ=LI(I)+2,LI(I)+NJM
            DFX(IJ)=DFX(IJ)/VOL
            DFY(IJ)=DFY(IJ)/VOL
        end do
    end do

end subroutine gradfi

!################################################################
    subroutine fluxSc(IJP, IJN, IJ1, IJ2, FM, FAC, CAP, CAN)
    !################################################################
        
        use bc
        use coef
        use flux
        use geo
        use sc
        use var
        implicit none
        real*8, intent(in) :: FM, FAC
        real*8, intent(in out) :: CAP, CAN
        integer, intent(in) :: IJP, IJN, IJ1, IJ2

        G=0.5d0

        !FACP=1.0d0-FAC
        !FII=T(IJN)*FAC+T(IJP)*FACP
        !DFXI=DTX(IJN)*FAC+DTX(IJP)*FACP
        !DFYI=DTY(IJN)*FAC+DTY(IJP)*FACP
        
    !
    !.....SURFACE AND DISTANCE VECTOR COMPONENTS, DIFFUSION COEFF.
    !
        SX=(Y(IJ1)-Y(IJ2))
        SY=(X(IJ2)-X(IJ1))
        XPN=XC(IJN)-XC(IJP)
        YPN=YC(IJN)-YC(IJP)
        VSOL=ALPHA*SQRT((SX**2+SY**2)/(XPN**2+YPN**2))
    !
    !.....EXPLICIT CONVECTIVE AND DIFFUSIVE FLUXES
    !
        !FCFIE=FM*FII
        !FDFIE=ALPHA*(DFXI*SX+DFYI*SY)
    !
    !.....IMPLICIT CONVECTIVE AND DIFFUSIVE FLUXES
    !
        !FCFII=MIN(FM,ZERO)*T(IJN)+MAX(FM,ZERO)*T(IJP)
        !FDFII=VSOL*(DFXI*XPN+DFYI*YPN)
    !
    !.....COEFFICIENTS, DEFERRED CORRECTION, SOURCE TERMS
    !
        CAN=-VSOL+MIN(FM,ZERO)
        CAP=-VSOL-MAX(FM,ZERO)
        !FFIC=G*(FCFIE-FCFII)
        !Q(IJP)=Q(IJP)-FFIC+FDFIE-FDFII
        !Q(IJN)=Q(IJN)+FFIC-FDFIE+FDFII

end subroutine fluxsc

!################################################################
subroutine temp
!################################################################

    use bc
    use coef
    use geo
    use ind
    use sc
    use var
    implicit none

    real*8 :: COEFC,COEFD,ZERO
    integer :: IJ1, IJ2
    ZERO=0.0d0

!
!.....DIRICHLET BOUNDARY CONDITION (NO WALL!)
!
      DO IW=1,NWALI
        IJP=IJPW(IW)
        IJB=IJW(IW)
        IJ1=IJW1(IW)
        IJ2=IJW2(IW)
        SX=(Y(IJ1)-Y(IJ2))
        SY=(X(IJ2)-X(IJ1))
        COEFC=RHO*(SX*VX+SY*VY)
        COEFD=ALPHA*SRDW(IW)
        AP(IJP)=AP(IJP)+COEFD
        Q(IJP)=Q(IJP)+COEFD*T(IJB)-min(COEFC*T(IJB),ZERO)
        !Q(IJP)=Q(IJP)+COEFD*T(IJB)-COEFC*T(IJB)
      END DO
      !print *, COEFF

end subroutine temp

!################################################################
subroutine calcErr
!################################################################

    use geo
    use ind
    use mms
    use var
    implicit none

    real*8 :: E,ER
    E=0.0d0


    do I=2,NIM
        do IJ=LI(I)+2,LI(I)+NIM
            !ER=T(IJ)-phi(XC(IJ),YC(IJ),0.0d0,TIME)
            E=E+abs(T(IJ)-phi(XC(IJ),YC(IJ),0.0d0,TIME))
            !print *, T(IJ), XC(IJ), YC(IJ), phi(XC(IJ),YC(IJ),0.0d0,TIME),ER
        end do
    end do
    
    write(10, *), E/N
    print *, E/N

end subroutine calcErr

