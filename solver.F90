!############### BLOCK STRUCUTRED VERSION ################
!#########################################################
program main
!#########################################################

    use charModule
    use coefModule
    use indexModule
    use controlModule
    use petsc_ksp_module
    use scalarModule
    use varModule
    implicit none
#include <finclude/petscsys.h>

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

    open(unit=9,FILE='ERR.out')
    rewind 9

    print *, 'STARTING SOLVER'
    call init
    print *, 'SETTING UP KSP'
    call setUpKSP

    ITIMS=1
    ITIME=1
!....START TIME LOOP
    print *, 'STARTING TIMELOOP'
    do ITIM=ITIMS,ITIME
        if (LTIME) then
            TIME=TIME+DT
            do B=1,NB
                call setBlockInd(B)
                do K=1,NK
                do I=1,NI
                do J=1,NJ
                    !IJK=LK(K)+LI(I)+J
                    IJK=IJKST+(K-1)*NK*NJ+(I-1)*NJ+J
                    TO(IJK)=T(IJK)
                end do
                end do
                end do
            end do
       end if

       print *, 'UPDATING DIRICHLET BOUNDARIES'
       call updateDir

!....START OUTER ITERATIONS
        LSG=10000
        !LSG=2
        do LS=1,LSG
            print *, 'OUTER ITERATION: ', LS
            print *, '  UPDATING GHOST VALUES'
            call updateGhost
            !print *, '  CALCULATING VELOCITY FIELD'
            !call calcuvw
            print *, '  SOLVING TRANSPORT EQUATION'
            call calcsc
            if (CONVERGED) then
                call writeVtk
                exit
            end if
        end do
        !call setField
    end do

    ! Free work space

    call cleanUp(Amat,bvec,solvec)
    call PetscFinalize(ierr)

end program main

!#########################################################
subroutine init
!#########################################################
    
    use boundaryModule
    use coefModule
    use charModule
    use geoModule
    use indexModule
    !use preProcInd
    use controlModule
    use parameterModule
    use varModule
    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscksp.h>
#include <finclude/petscpc.h>

    PetscErrorCode :: ierr

    PROC=0
    write(PROC_CH,'(I1)') PROC
    PROCUNIT=PROC+PROCOFFSET
    PROCFILE='proc_'//trim(PROC_CH)//'.inp'

    open(UNIT=PROCUNIT,FILE=PROCFILE)
    rewind PROCUNIT 

    read(PROCUNIT,*) NB
    !print *, NB
    !
    !allocate(B_GLO(NB))
    !
    read(PROCUNIT,*) (B_GLO(B),B=1,NB)
    !print *, (B_GLO(B),B=1,NB)

    write(BLOCK_CH,'(I1)') B_GLO(1)
    BLOCKUNIT=BLOCKOFFSET+B_GLO(1)
    BLOCKFILE='grid_'//trim(BLOCK_CH)//'.out'
    open(UNIT=BLOCKUNIT,FILE=BLOCKFILE)
    rewind BLOCKUNIT
    read(BLOCKUNIT,*) NI,NJ,NK,NIJK,NBLOCK,NDIR,NFACE,N,IJKST
    
    IBL(1)=0
    JBL(1)=0
    KBL(1)=0
    IJKBL(1)=0
    IJKBLOCKBL(1)=0
    IJKDIRBL(1)=0
    NIBL(1)=NI
    NJBL(1)=NJ
    NKBL(1)=NK
    NIJKBL(1)=NIJK
    NBLOCKBL(1)=NBLOCK
    NDIRBL(1)=NDIR
    NFACEBL(1)=NFACE
    IJKPROC=IJKST
    NBL(1)=N

    do B=2,NB
        !read(PROCUNIT,*) B_GLO(B)
        BLOCKUNIT=BLOCKOFFSET+B_GLO(B)
        write(BLOCK_CH,'(I1)') B_GLO(B)
        BLOCKFILE='grid_'//trim(BLOCK_CH)//'.out'
        !print *, BLOCKFILE
        open(UNIT=BLOCKUNIT,FILE=BLOCKFILE)
        read(BLOCKUNIT,*) NI,NJ,NK,NIJK,NBLOCK,NDIR,NFACE,N,IJKST

        BB=B-1
        IBL(B)=IBL(BB)+NIBL(BB)
        JBL(B)=JBL(BB)+NJBL(BB)
        KBL(B)=KBL(BB)+NKBL(BB)
        IJKBL(B)=IJKBL(BB)+NIJKBL(BB)
        IJKBLOCKBL(B)=IJKBLOCKBL(BB)+NBLOCKBL(BB)
        IJKDIRBL(B)=IJKDIRBL(BB)+NDIRBL(BB)
        FACEBL(B)=FACEBL(BB)+NFACEBL(BB)
        NIBL(B)=NI
        NJBL(B)=NJ
        NKBL(B)=NK
        NIJKBL(B)=NIJK
        NBLOCKBL(B)=NBLOCK
        NDIRBL(B)=NDIR
        NFACEBL(B)=NFACEBL(BB)
        !IJKBL_GLO(B)=IJKST
        NBL(B)=N
    end do
    
    ! calculate processor load
    N=sum(NBL)

    print *, N

    do B=1,NB
        BLOCKUNIT=BLOCKOFFSET+B_GLO(B)
        call setBlockInd(B)
        !read(BLOCKUNIT,*) (LK(KST+K),K=1,NK)
        !read(BLOCKUNIT,*) (LI(IST+I),I=1,NI)
        read(BLOCKUNIT,*) (X(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT,*) (Y(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT,*) (Z(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT,*) (XC(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT,*) (YC(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT,*) (ZC(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT,*) (IJKBBL(IJKBLOCKST+IJK),IJK=1,NBLOCK)
        read(BLOCKUNIT,*) (IJKPBL(IJKBLOCKST+IJK),IJK=1,NBLOCK)
        read(BLOCKUNIT,*) (IJKBL1(IJKBLOCKST+IJK),IJK=1,NBLOCK)
        read(BLOCKUNIT,*) (IJKBL2(IJKBLOCKST+IJK),IJK=1,NBLOCK)
        read(BLOCKUNIT,*) (IJKBL3(IJKBLOCKST+IJK),IJK=1,NBLOCK)
        read(BLOCKUNIT,*) (IJKBL4(IJKBLOCKST+IJK),IJK=1,NBLOCK)
        read(BLOCKUNIT,*) (IJKBDI(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKPDI(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKDI1(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKDI2(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKDI3(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKDI4(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (FX(IJKST+IJK), IJK=1,NIJK)
        read(BLOCKUNIT,*) (FY(IJKST+IJK), IJK=1,NIJK)
        read(BLOCKUNIT,*) (FZ(IJKST+IJK), IJK=1,NIJK)
        read(BLOCKUNIT,*) DXBL(B),DYBL(B),DZBL(B),VOLBL(B)
        read(BLOCKUNIT,*) (SRDDI(IJKDIRST+IJK),IJK=1,NDIR)
        !19
        read(BLOCKUNIT,*) (L(FACEST+I),I=1,NFACE)
        read(BLOCKUNIT,*) (R(FACEST+I),I=1,NFACE)
        read(BLOCKUNIT,*) (XF(FACEST+I),I=1,NFACE)
        read(BLOCKUNIT,*) (YF(FACEST+I),I=1,NFACE)
        read(BLOCKUNIT,*) (ZF(FACEST+I),I=1,NFACE)
        read(BLOCKUNIT,*) (FF(FACEST+I),I=1,NFACE)
        read(BLOCKUNIT,*) (ARF(FACEST+I),I=1,NFACE)
        read(BLOCKUNIT,*) (DNF(FACEST+I),I=1,NFACE)
        read(BLOCKUNIT,*) (XPNF(FACEST+I),I=1,NFACE)
        read(BLOCKUNIT,*) (YPNF(FACEST+I),I=1,NFACE)
        read(BLOCKUNIT,*) (ZPNF(FACEST+I),I=1,NFACE)
        read(BLOCKUNIT,*) (NXF(FACEST+I),I=1,NFACE)
        read(BLOCKUNIT,*) (NYF(FACEST+I),I=1,NFACE)
        read(BLOCKUNIT,*) (NZF(FACEST+I),I=1,NFACE)
    end do
    
    BLOCKUNIT=BLOCKOFFSET+B_GLO(1)
    read(BLOCKUNIT,*) (MIJK(IJK),IJK=1,NXYZA)

    do B=1,NB
        BLOCKUNIT=BLOCKOFFSET+B_GLO(B)
        close(UNIT=BLOCKUNIT)
    end do

    TIME=0.0d0

    if (LTIME) then
        call setField
        !call writeVtk
    end if

    ! Create Matrix
    call MatCreate(PETSC_COMM_WORLD,Amat,ierr)
    call MatSetType(Amat,MATSEQAIJ,ierr)
    !call MatSetType(Amat,MATMPIAIJ,ierr)
    call MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,N,N,ierr)
    !call MatSetSizes(Amat,N,N,PETSC_DECIDE,PETSC_DECIDE,ierr)
    call MatSetFromOptions(Amat,ierr) 
    
    ! Increase performance during matrix assembly due to preallocation
    call MatSeqAIJSetPreallocation(Amat,i7,PETSC_NULL_INTEGER,ierr)
    !call MatMPIAIJSetPreallocation(Amat,i7,PETSC_NULL_INTEGER,ierr)

    ! Create Vector
    call VecCreate(PETSC_COMM_WORLD,solvec,ierr)
    call VecSetSizes(solvec,PETSC_DECIDE,N,ierr)
    !call VecSetSizes(solvec,N,PETSC_DECIDE,ierr)
    call VecSetFromOptions(solvec,ierr)
    call VecDuplicate(solvec,bvec,ierr)

    !stop

end subroutine init

!#########################################################
subroutine setField
!#########################################################

    use geoModule
    use indexModule
    use mmsModule
    use varModule
    implicit none

    do B=1,NB
        call setBlockInd(B)
        do K=1,NK
        do I=1,NI
        do J=1,NJ
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            T(IJK)=phi(XC(IJK),YC(IJK),ZC(IJK),TIME)
        end do
        end do
        end do
    end do

end subroutine setField

!#########################################################
subroutine writeVtk
!#########################################################

    use charModule
    use geoModule
    use indexModule
    use controlModule
    use varModule
    implicit none

    do B=1,NB
        call setBlockInd(B)
        BLOCKUNIT=B_GLO(B)+BLOCKOFFSET
        !write(TIME_CH,'(f6.4)') TIME
        write(TIME_CH,'(I1)') ITIM
        write(BLOCK_CH,'(I1)') B_GLO(B)
        !write(TIME_CH,'(I1)') LS
        VTKFILE='grid_'//trim(BLOCK_CH)//'_'//trim(TIME_CH)//'.vtk'
        print *, ' *** GENERATING .VTK *** '
        !
        open (UNIT=BLOCKUNIT,FILE=VTKFILE)
        rewind BLOCKUNIT
        write(BLOCKUNIT,'(A)') '# vtk DataFile Version 3.0'
        write(BLOCKUNIT,'(A)') 'grid'
        write(BLOCKUNIT,'(A)') 'ASCII'
        write(BLOCKUNIT,'(A)') 'DATASET STRUCTURED_GRID'
        write(BLOCKUNIT,'(A I6 I6 I6)') 'DIMENSIONS', NIM, NJM,NKM
        write(BLOCKUNIT,'(A I9 A)') 'Points', NIM*NJM*NKM, ' float'
        !
        do K=1,NKM
        do I=1,NIM
        do J=1,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            write(BLOCKUNIT,'(E20.10,1X,E20.10,1X,E20.10)'), X(IJK), Y(IJK),Z(IJK)
        end do
        end do
        end do
        !
        write(BLOCKUNIT,'(A10, I9)') 'CELL_DATA ',(NICV*NJCV*NKCV)
        write(BLOCKUNIT,'(A15)') 'SCALARS T float'
        write(BLOCKUNIT,'(A20)') 'LOOKUP_TABLE default'
        !
        do K=2,NKM
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            write(BLOCKUNIT,'(F12.8)') T(IJK)
        end do
        end do
        end do
        close(UNIT=BLOCKUNIT)
    end do

end subroutine writeVtk

!########################################################
subroutine updateDir
!#########################################################

    use boundaryModule
    use geoModule
    use indexModule
    use mmsModule
    use scalarModule
    use parameterModule
    use varModule
    implicit none

    do B=1,NB
        call setBlockInd(B)
    !...Dirichlet BC
        do IJKDIR=IJKDIRST+1,IJKDIRST+NDIR
            IJKB=IJKBDI(IJKDIR)
            T(IJKB)=phi(XC(IJKB),YC(IJKB),ZC(IJKB),TIME)
        end do
    end do

end subroutine updateDir

!########################################################
subroutine updateGhost
!#########################################################

    use boundaryModule
    use geoModule
    use indexModule
    use mmsModule
    use scalarModule
    use parameterModule
    use varModule
    implicit none

    do B=1,NB
        call setBlockInd(B)

        !...Block BC
        do F=FACEST+1,FACEST+NFACE
            TR(F)=T(R(F)) 
        end do

        !...Calculate Mass Fluxes at all(!) cell faces
        do K=2,NKM
        !do I=2,NIM-1
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            F1(IJK)=RHO*DY*DZ*VX
        end do
        end do
        end do

        do K=2,NKM
        do I=2,NIM
        !do J=2,NJM-1
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            F2(IJK)=RHO*DX*DZ*VY
        end do
        end do
        end do

        !do K=2,NKM-1
        do K=2,NKM
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            F3(IJK)=RHO*DX*DY*VZ
        end do
        end do
        end do
    end do
    
end subroutine updateGhost

!#########################################################
subroutine calcSc
!#########################################################

    use boundaryModule
    use coefModule
    use geoModule
    use gradModule
    use indexModule
    !use preProcInd
    use mmsModule
    use petsc_ksp_module
    use scalarModule
    use controlModule
    use varModule
    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>

    real*8 :: APT, URF
    PetscLogDouble :: time1, time2

    URF=1.0d0

    !print *, '  CALCULATE CV-CENTER GRADIENTS'
    call gradfi(T,TR,DTX,DTY,DTZ)
    call updateGrad

    ! START BLOCK LOOP
    do B=1,NB
        call setBlockInd(B)

        ! initialize Q and AP
        do K=2,NKM
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            Q(IJK)=src(XC(IJK), YC(IJK), ZC(IJK), TIME)*VOL
            AP(IJK)=0.0d0
        end do
        end do
        end do

        !do F=FACEST+1,FACEST+NFACE
        !    AF(F)=0.0d0
        !end do

        ! FLUXES THROUGH EAST FACE
        do K=2,NKM
        do I=2,NIM-1
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            call fluxsc(IJK,IJK+NJ,IJK-NIJ-1,IJK-1,IJK,F1(IJK),FX(IJK),AW(IJK+NJ),AE(IJK))
        end do
        end do
        end do


        ! FLUXES THROUGH NORTH FACE
        do K=2,NKM
        do I=2,NIM
        do J=2,NJM-1
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            call fluxsc(IJK,IJK+1,IJK-NIJ,IJK,IJK-NJ,F2(IJK),FY(IJK),AS(IJK+1),AN(IJK))
        end do
        end do
        end do

        ! FLUXES THROUGH TOP FACE
        do K=2,NKM-1
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
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
                IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
                APT=RHO*VOL/DT
                Q(IJK)=Q(IJK)+APT*TO(IJK)
                AP(IJK)=AP(IJK)+APT
            end do
            end do
            end do
        end if

    !
    !.....DIRICHLET BOUNDARY CONTRIBUTION
    !
        call dirBdFlux

    !
    !.....BLOCK BOUNDARY CONTRIBUTION
    !
        call blockBdFlux

    !
    !.....FINAL COEFFICIENT AND SOURCE MATRIX FOR FI-EQUATION
    !
        do K=2,NKM
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            AP(IJK)=(AP(IJK)-AE(IJK)-AW(IJK)-AN(IJK)-AS(IJK)-AT(IJK)-AB(IJK))/URF
            Q(IJK)=Q(IJK)+(1.0d0-URF)*AP(IJK)*T(IJK)
        end do
        end do
        end do
        
        ! Assemble matrix and Vector

        ! Assemble inner CV matrix coefficients
        print *, 'Assemble inner CV matrix coefficients'
        do K=3,NKM-1
        do I=3,NIM-1
        do J=3,NJM-1
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            IJKP=MIJK(IJK)
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
            !print *, 'TEST'
            !print *, row, col
            !
            call MatSetValues(Amat,i1,row,i7,col,val,INSERT_VALUES,ierr)
            call VecSetValue(bvec,row,valq,INSERT_VALUES,ierr)
            !
        end do
        end do
        end do

        ! Assembly matrix coefficients of dirichlet boundaries
        do IJKDIR=IJKDIRST+1,IJKDIRST+NDIR
            IJK=IJKPDI(IJKDIR)
            IJKP=MIJK(IJK)
            row=IJKP
            col=(/-1,-1,-1,IJKP,-1,-1,-1/)
            !print *, IJKDIR,IJKPDI(IJKDIR),row
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
            call MatSetValues(Amat,i1,row,i7,col,val,INSERT_VALUES,ierr)
            call VecSetValue(bvec,row,valq,INSERT_VALUES,ierr)
            !
        end do

        ! Assembly matrix coefficients of block boundary cells
        do IJKBLOCK=IJKBLOCKST+1,IJKBLOCKST+NBLOCK
            IJK=IJKPBL(IJKBLOCK)
            IJKP=MIJK(IJK)
            row=IJKP
            col=(/-1,-1,-1,IJKP,-1,-1,-1/)
            !print *, IJKBLOCK,IJKPBL(IJKBLOCK),row
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
            call MatSetValues(Amat,i1,row,i7,col,val,INSERT_VALUES,ierr)
            call VecSetValue(bvec,row,valq,INSERT_VALUES,ierr)
            !
        end do

        ! assembly matrix coefficients of block boundaries
        do F=FACEST+1,FACEST+NFACE
            row=MIJK(L(F))
            col1=MIJK(R(F))
            !print *, F,row,col1
            val1=AF(F)
            !call MatSetValues(Amat,i1,col1,val1,INSERT_VALUES,ierr)
            call MatSetValue(Amat,row,col1,val1,INSERT_VALUES,ierr)
        end do

    end do

    ! assembly matrix and right hand vector

    print *, '  STARTING MATRIX ASSEMBLY'
    call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr)
    call VecAssemblyBegin(bvec,ierr)
    call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr)
    call VecAssemblyEnd(bvec,ierr)

    !call MatView(Amat,PETSC_VIEWER_STDOUT_WORLD,ierr)
    !print *, ''
    !call VecView(bvec,PETSC_VIEWER_STDOUT_WORLD,ierr)
    !stop

    print *, '  SOLVING LINEAR SYSTEM'
    call PetscGetCPUTime(time1,ierr)
    call solveSys(Amat,bvec,solvec,N,LS,tol)
    call PetscGetCPUTime(time2,ierr)

    if (CONVERGED) then
        return
    end if

    tges=time2-time1

    do B=1,NB
        call setBlockInd(B)
        do K=2,NKM
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            row=MIJK(IJK)
            call VecGetValues(solvec,i1,row,valt,ierr)
            T(IJK)=valt
        end do
        end do
        end do
    end do

    call calcErr

end subroutine calcSc

!################################################################
subroutine gradfi(FI,FIR,DFX,DFY,DFZ)
!################################################################

    use boundaryModule
    use geoModule
    use gradModule
    use indexModule
    use parameterModule
    use scalarModule
    implicit none

    real(KIND=PREC), intent(in out) :: FI(NIJK),FIR(NFACE),DFX(NIJK),DFY(NIJK),DFZ(NIJK)
    integer :: IJK1, IJK2, IJK3, IJK4

    do B=1,NB
        call setBlockInd(B)
    !
    !...INITIALIZE FIELDS
    !
        !print *, '      INITIALIZING FIELDS'
        do IJK=IJKST+1,IJKST+NIJK
            DFX(IJK)=0.0d0
            DFY(IJK)=0.0d0
            DFZ(IJK)=0.0d0
        end do

    !
    !..CONTRRIBUTION FROM INNER EAST SIDES
    !
        !print *, '      CALC CONTRIBUTION FROM INNER EAST SIDES'
        do K=2,NKM
        do I=2,NIM-1
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            FIE=FI(IJK+NJ)*FX(IJK)+FI(IJK)*(1.0d0-FX(IJK))

            IJK4=IJK
            IJK3=IJK4-1
            IJK2=IJK3-NIJ
            IJK1=IJK4-NIJ
            !
            call normalArea(IJK,IJK,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
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
        !print *, '      CALC CONTRIBUTION FROM INNER NORTH SIDES'
        do K=2,NKM
        do I=2,NIM
        do J=2,NJM-1
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            FIN=FI(IJK+1)*FY(IJK)+FI(IJK)*(1.0d0-FY(IJK))

            IJK3=IJK
            IJK4=IJK3-NJ
            IJK2=IJK3-NIJ
            IJK1=IJK4-NIJ
            !
            call normalArea(IJK,IJK,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
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
        !print *, '      CALC CONTRIBUTION FROM INNER TOP SIDES'
        do K=2,NKM-1
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            FIN=FI(IJK+NIJ)*FZ(IJK)+FI(IJK)*(1.0d0-FZ(IJK))

            IJK4=IJK
            IJK3=IJK4-NJ
            IJK1=IJK4-1
            IJK2=IJK3-1
            !
            call normalArea(IJK,IJK,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
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
        !print *, '      CALC CONTRIBUTION FROM WALL BOUNDARIES'
        do IJKDIR=IJKDIRST+1,IJKDIRST+NDIR
            IJKB=IJKBDI(IJKDIR)
            IJKP=IJKPDI(IJKDIR)
            !IJK1=IJKDI1(IJKDIR)
            IJK2=IJKDI2(IJKDIR)
            IJK3=IJKDI3(IJKDIR)
            IJK4=IJKDI4(IJKDIR)
            !
            call normalArea(IJKP,IJKB,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
            !
            SX=AR*NX
            SY=AR*NY
            SZ=AR*NZ

            DFX(IJKPDI(IJKDIR))=DFX(IJKPDI(IJKDIR))+FI(IJKBDI(IJKDIR))*SX
            DFY(IJKPDI(IJKDIR))=DFY(IJKPDI(IJKDIR))+FI(IJKBDI(IJKDIR))*SY
            DFZ(IJKPDI(IJKDIR))=DFZ(IJKPDI(IJKDIR))+FI(IJKBDI(IJKDIR))*SZ
        end do
    !
    !.....CONTRIBUTION FROM BLOCK BOUNDARIES
    !
        !print *, '      CALC CONTRIBUTION FROM BLOCK BOUNDARIES'
        do F=FACEST+1,FACEST+NFACE
            FIN=FIR(F)*FF(F)+FI(L(F))*(1.0d0-FF(F))
            
            SX=ARF(F)*NXF(F)
            SY=ARF(F)*NYF(F)
            SZ=ARF(F)*NZF(F)
            
            DFXN=FIN*SX
            DFYN=FIN*SY
            DFZN=FIN*SZ
            
            DFX(L(F))=DFX(L(F))+DFXN
            DFY(L(F))=DFY(L(F))+DFYN
            DFZ(L(F))=DFZ(L(F))+DFZN
        end do
    !
    !.....CALCULATE GRADIENT COMPONENTS AT CV-CENTERS
    !
        do K=2,NKM
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            DFX(IJK)=DFX(IJK)/VOL
            DFY(IJK)=DFY(IJK)/VOL
            DFZ(IJK)=DFZ(IJK)/VOL
        end do
        end do
        end do
    end do

end subroutine gradfi

!################################################################
subroutine updateGrad
!################################################################

    use boundaryModule
    !use geoModule
    use indexModule
    !use mmsModule
    !use scalarModule
    use parameterModule
    use varModule
    implicit none

    do B=1,NB
        call setBlockInd(B)
        do F=FACEST+1,FACEST+NFACE
          DTXR(F)=DTX(R(F))  
          DTYR(F)=DTY(R(F))
          DTZR(F)=DTZ(R(F))
        end do
    end do

end subroutine updateGrad

!################################################################
subroutine fluxSc(IJKP,IJKN,IJK2,IJK3,IJK4,FM,FAC,CAP,CAN)
!################################################################
        
    use boundaryModule
    use coefModule
    use fluxModule
    use geoModule
    use indexModule, only: MIJK
    use scalarModule
    use parameterModule
    use varModule
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
    !
    call normalArea(IJKP,IJKN,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
    !
    VSOL=ALPHA*AR/(DN+SMALL)
    !
    !.....EXPLICIT CONVECTIVE AND DIFFUSIVE FLUXES
    !
    FCFIE=FM*FII
    FDFIE=ALPHA*(DFXI*AR*NX+DFYI*AR*NY+DFZI*AR*NZ)
    !
    !.....IMPLICIT CONVECTIVE AND DIFFUSIVE FLUXES
    !
    FCFII=MIN(FM,ZERO)*T(IJKN)+MAX(FM,ZERO)*T(IJKP)
    FDFII=VSOL*(DFXI*XPN+DFYI*YPN+DFZI*ZPN)
    !print *, MIJK(IJKP),T(IJKN),T(IJKP),FII,FCFIE
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
subroutine dirBdFlux
!################################################################

    use boundaryModule
    use coefModule
    use geoModule
    use indexModule
    use scalarModule
    use parameterModule
    use varModule
    implicit none

    real(KIND=PREC) :: COEFC,COEFD,ZERO
    integer :: IJK1, IJK2, IJK3, IJK4
    ZERO=0.0d0

!
!.....DIRICHLET BOUNDARY CONDITION (NO WALL!)
!
    do IJKDIR=IJKDIRST+1,IJKDIRST+NDIR
        IJKP=IJKPDI(IJKDIR)
        IJKB=IJKBDI(IJKDIR)
        !IJK1=IJKDI1(IJKDIR)
        IJK2=IJKDI2(IJKDIR)
        IJK3=IJKDI3(IJKDIR)
        IJK4=IJKDI4(IJKDIR)
        !
        call normalArea(IJKB,IJKP,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
        !
        SX=AR*NX
        SY=AR*NY
        SZ=AR*NZ
        !
        COEFC=RHO*(SX*VX+SY*VY+SZ*VZ)
        COEFD=ALPHA*SRDDI(IJKDIR)
        AP(IJKP)=AP(IJKP)+COEFD-COEFC
        Q(IJKP)=Q(IJKP)+(COEFD-COEFC)*T(IJKB)
        !print *, MIJK(IJKP), IJKB, T(IJKB), COEFD, COEFC
    end do

end subroutine dirBdFlux

!################################################################
subroutine blockBdFlux
!################################################################

    use boundaryModule
    use coefModule
    use fluxModule
    use geoModule
    use indexModule
    use scalarModule
    use parameterModule
    use varModule
    implicit none
    real(KIND=PREC) :: ZERO,SMALL
    real(KIND=PREC) :: FAC,FM

    ZERO=0.0d0
    SMALL=1.0E-20

    do F=FACEST+1,FACEST+NFACE
        FAC=FF(F)
        FACP=1.0d0-FAC
        FII=TR(F)*FAC+T(L(F))*FACP
        DFXI=DTXR(F)*FAC+DTX(L(F))*FACP
        DFYI=DTYR(F)*FAC+DTY(L(F))*FACP
        DFZI=DTZR(F)*FAC+DTZ(L(F))*FACP

        FM=F1(L(F))*NXF(F)+F2(L(F))*NYF(F)+F3(L(F))*NZF(F)

        VSOL=ALPHA*ARF(F)/(DNF(F)+SMALL)
        !
        !.....EXPLICIT CONVECTIVE AND DIFFUSIVE FLUXES
        !
        FCFIE=FM*FII
        FDFIE=ALPHA*(DFXI*ARF(F)*NXF(F)+DFYI*ARF(F)*NYF(F)+DFZI*ARF(F)*NZF(F))
        !
        !.....IMPLICIT CONVECTIVE AND DIFFUSIVE FLUXES
        !
        FCFII=MIN(FM,ZERO)*TR(F)+MAX(FM,ZERO)*T(L(F))
        FDFII=VSOL*(DFXI*XPNF(F)*NXF(F)+DFYI*YPNF(F)*NYF(F)+DFZI*ZPNF(F)*NZF(F))
        !print *, MIJK(L(F)),FDFII,FDFIE
        !
        !.....COEFFICIENTS, DEFERRED CORRECTION, SOURCE TERMS
        !
        !CAN=-VSOL+MIN(FM,ZERO)
        !CAP=-VSOL-MAX(FM,ZERO)
        !AF(F)=-VSOL-MAX(FM,ZERO)
        ! Boundary must be treated inoutflow like
        AF(F)=-VSOL+FM
        AP(L(F))=AP(L(F))-AF(F)
        FFIC=G*(FCFIE-FCFII)
        Q(L(F))=Q(L(F))-FFIC+FDFIE-FDFII
        !Q(L(F))=Q(L(F))-FFIC
        !Q(IJKN)=Q(IJKN)+FFIC-FDFIE+FDFII
    end do

end subroutine blockBdFlux

!################################################################
subroutine calcErr
!################################################################

    use geoModule
    use indexModule
    use controlModule
    use mmsModule
    use parameterModule
    use varModule
    implicit none

    real(KIND=PREC) :: E,ER
    E=0.0d0
    ER=0.0d0

    do B=1,NB
        call setBlockInd(B)
        do K=2,NKM
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            !ER=T(IJ)-phi(XC(IJ),YC(IJ),0.0d0,TIME)
            !E=abs(T(IJ)-phi(XC(IJ),YC(IJ),0.0d0,TIME))
            E=E+abs(T(IJK)-phi(XC(IJK),YC(IJK),ZC(IJK),TIME))
            !ER=max(E,ER)
        end do
        end do
        end do
    end do
    
    rewind 9
    write(9, *), E/dble(N), tges, reasonInt, itsInt, LS
    print *, 'ERROR ', E/dble(N), tges, reasonInt, itsInt, LS

end subroutine calcErr
