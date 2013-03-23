!############### BLOCK STRUCUTRED VERSION ################
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

    !open(unit=2,FILE='grid.out')
    !rewind 2

    !open(unit=10,FILE='ERR.out')
    !rewind 10

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
                !IJK=LK(K)+LI(I)+J
                IJK=(K-1)*NK*NI+(I-1)*NI+J
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

    call cleanUp(Amat,bvec,solvec)
    call PetscFinalize(ierr)

end program main

!#########################################################
subroutine init
!#########################################################
    
    use bc
    use coef
    use geo
    !use ind
    use preProcInd
    use logic
    use param
    use var
    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscksp.h>
#include <finclude/petscpc.h>

    PetscErrorCode :: ierr
    integer :: BLOCKUNIT,PROC,PROCUNIT
    !integer, dimension(:), allocatable :: B_GLO
    character(len=20) :: BLOCKFILE,BLOCK_CH,PROCFILE,PROC_CH

    PROC=0
    write(PROC_CH,'(I1)') PROC
    PROCUNIT=PROC+PROCOFFSET
    PROCFILE='proc_'//trim(PROC_CH)//'.inp'

    open(UNIT=PROCUNIT,FILE=PROCFILE)
    rewind PROCUNIT 

    read(PROCUNIT,*) NB
    print *, NB
    !
    !allocate(B_GLO(NB))
    !
    read(PROCUNIT,*) (B_GLO(B),B=1,NB)
    print *, (B_GLO(B),B=1,NB)

    write(BLOCK_CH,'(I1)') B_GLO(1)
    BLOCKUNIT=BLOCKOFFSET+B_GLO(1)
    BLOCKFILE='grid_'//trim(BLOCK_CH)//'.out'
    open(UNIT=BLOCKUNIT,FILE=BLOCKFILE)
    rewind BLOCKUNIT
    read(BLOCKUNIT,*)   NI,NJ,NK,NIJK,NDIR,NFACE,N,IJKST
    
    IBL(1)=0
    JBL(1)=0
    KBL(1)=0
    IJKBL(1)=0
    IJKDIRBL(1)=0
    NIBL(1)=NI
    NJBL(1)=NJ
    NKBL(1)=NK
    NIJKBL(1)=NIJK
    NDIRBL(1)=NDIR
    NFACEBL(1)=NFACE
    IJKBL_GLO(1)=IJKST
    NBL(1)=N

    do B=2,NB
        !read(PROCUNIT,*) B_GLO(B)
        BLOCKUNIT=BLOCKOFFSET+B_GLO(B)
        write(BLOCK_CH,'(I1)') B_GLO(B)
        BLOCKFILE='grid_'//trim(BLOCK_CH)//'.out'
        print *, BLOCKFILE
        open(UNIT=BLOCKUNIT,FILE=BLOCKFILE)
        read(BLOCKUNIT,*)   NI,NJ,NK,NIJK,NDIR,NFACE,N,IJKST

        BB=B-1
        IBL(B)=IBL(BB)+NIBL(BB)
        JBL(B)=JBL(BB)+NJBL(BB)
        KBL(B)=KBL(BB)+NKBL(BB)
        IJKBL(B)=IJKBL(BB)+NIJKBL(BB)
        IJKDIRBL(B)=IJKDIRBL(BB)+NDIRBL(BB)
        FACEBL(B)=FACEBL(BB)+NFACEBL(BB)
        NIBL(B)=NI
        NJBL(B)=NJ
        NKBL(B)=NK
        NIJKBL(B)=NIJK
        NDIRBL(B)=NDIR
        NFACEBL(B)=NFACEBL(BB)
        IJKBL_GLO(B)=IJKST
        NBL(B)=N
    end do
    
    ! calculate processor load
    N=sum(NBL)

    print *, N

    !call indAllocate

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
        read(BLOCKUNIT,*) (IJKBDI(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKPDI(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKDI1(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKDI2(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKDI3(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKDI4(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (FX(IJKST+IJK), IJK=1,NIJK)
        read(BLOCKUNIT,*) (FY(IJKST+IJK), IJK=1,NIJK)
        read(BLOCKUNIT,*) (FZ(IJKST+IJK), IJK=1,NIJK)
        read(BLOCKUNIT,*) !DX(B),DY(B),DZ(B),VOL(B)
        read(BLOCKUNIT,*) (SRDDI(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (L(FACEST+I),I=1,NFACE)
        read(BLOCKUNIT,*) !(R(FACEST+I),I=1,NFACE)
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

    TIME=1.0d0

    if (LTIME) then
        call setField
        call writeVtk
    end if

    ! Create Matrix
    call MatCreate(PETSC_COMM_WORLD,Amat,ierr)
    call MatSetType(Amat,MATSEQAIJ,ierr)
    call MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,N,N,ierr)
    call MatSetFromOptions(Amat,ierr) 
    
    ! Increase performance during matrix assembly due to preallocation
    call MatSeqAIJSetPreallocation(Amat,i7,PETSC_NULL_INTEGER,ierr)

    ! Create Vector
    call VecCreate(PETSC_COMM_WORLD,solvec,ierr)
    call VecSetSizes(solvec,PETSC_DECIDE,N,ierr)
    call VecSetFromOptions(solvec,ierr)
    call VecDuplicate(solvec,bvec,ierr)

    stop

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
        IJK=(K-1)*NI*NJ+(I-1)*NI+J
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
    do IDI=1,NDIRA
        IJKB=IJKBDI(IDI)
        T(IJKB)=phi(XC(IJKB),YC(IJKB),ZC(IJKB),TIME)
    end do

!...Calculate MassFluxes
    do K=2,NKM
    do I=2,NIM-1
    do J=2,NJM
        IJK=(K-1)*NI*NJ+(I-1)*NI+J
        F1(IJK)=RHO*DY*DZ*VX
    end do
    end do
    end do

    do K=2,NKM
    do I=2,NIM
    do J=2,NJM-1
        IJK=(K-1)*NI*NJ+(I-1)*NI+J
        F2(IJK)=RHO*DX*DZ*VY
    end do
    end do
    end do

    do K=2,NKM-1
    do I=2,NIM
    do J=2,NJM
        IJK=(K-1)*NI*NJ+(I-1)*NI+J
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
    use ind, only : NDIRA,NJCV,NIJCV,LS
    use preProcInd
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
        IJK=(K-1)*NI*NJ+(I-1)*NI+J
        Q(IJK)=src(XC(IJK), YC(IJK), ZC(IJK), TIME)*VOL
        AP(IJK)=0.0d0
    end do
    end do
    end do

    ! FLUXES THROUGH EAST FACE
    do K=2,NKM
    do I=2,NIM-1
    do J=2,NJM
        IJK=(K-1)*NI*NJ+(I-1)*NI+J
        call fluxsc(IJK,IJK+NJ,IJK-NIJ-1,IJK-1,IJK,F1(IJK),FX(IJK),AW(IJK+NJ),AE(IJK))
    end do
    end do
    end do


    ! FLUXES THROUGH NORTH FACE
    do K=2,NKM
    do I=2,NIM
    do J=2,NJM-1
        IJK=(K-1)*NI*NJ+(I-1)*NI+J
        call fluxsc(IJK,IJK+1,IJK-NIJ,IJK,IJK-NJ,F2(IJK),FY(IJK),AS(IJK+1),AN(IJK))
    end do
    end do
    end do

    ! FLUXES THROUGH TOP FACE
    do K=2,NKM-1
    do I=2,NIM
    do J=2,NJM
        IJK=(K-1)*NI*NJ+(I-1)*NI+J
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
            IJK=(K-1)*NI*NJ+(I-1)*NI+J
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
        IJK=(K-1)*NI*NJ+(I-1)*NI+J
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
        IJK=(K-1)*NI*NJ+(I-1)*NI+J
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
        !
        call MatSetValues(Amat,i1,row,i7,col,val,INSERT_VALUES,ierr)
        call VecSetValue(bvec,row,valq,INSERT_VALUES,ierr)
        !
    end do
    end do
    end do

    ! Assembly matrix coefficients of boundary cv

    do IDI=1,NDIRA
        IJK=IJKPDI(IDI)
        IJKP=MIJK(IJK)
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
        call MatSetValues(Amat,i1,row,i7,col,val,INSERT_VALUES,ierr)
        call VecSetValue(bvec,row,valq,INSERT_VALUES,ierr)
        !
    end do


    ! assembly matrix and right hand vector

    call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr)
    call VecAssemblyBegin(bvec,ierr)
    call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr)
    call VecAssemblyEnd(bvec,ierr)

    call PetscGetCPUTime(time1,ierr)
    call solveSys(Amat,bvec,solvec,N,LS,tol)
    call PetscGetCPUTime(time2,ierr)

    if (CONVERGED) then
        return
    end if

    tges=time2-time1

    do K=2,NIM
    do I=2,NIM
        do J=2,NJM
            IJK=(K-1)*NI*NJ+(I-1)*NI+J
            row=MIJK(IJK)
            call VecGetValues(solvec,i1,row,valt,ierr)
            T(IJK)=valt
        end do
    end do
    end do

    !call calcErr

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
        IJK=(K-1)*NI*NJ+(I-1)*NI+J
        FIE=FI(IJ+NJ)*FX(IJ)+FI(IJ)*(1.0d0-FX(IJ))

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
    do K=2,NKM
    do I=2,NIM
    do J=2,NJM-1
        IJK=(K-1)*NI*NJ+(I-1)*NI+J
        FIN=FI(IJ+1)*FY(IJ)+FI(IJ)*(1.0d0-FY(IJ))

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
    do K=2,NKM-1
    do I=2,NIM
    do J=2,NJM
        IJK=(K-1)*NI*NJ+(I-1)*NI+J
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
    do IDI=1,NDIRA
        IJKB=IJKBDI(IDI)
        IJKP=IJKPDI(IDI)
        !IJK1=IJKDI1(IDI)
        IJK2=IJKDI2(IDI)
        IJK3=IJKDI3(IDI)
        IJK4=IJKDI4(IDI)
        !
        call normalArea(IJKP,IJKB,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
        !
        SX=AR*NX
        SY=AR*NY
        SZ=AR*NZ

        DFX(IJKPDI(IDI))=DFX(IJKPDI(IDI))+FI(IJKBDI(IDI))*SX
        DFY(IJKPDI(IDI))=DFY(IJKPDI(IDI))+FI(IJKBDI(IDI))*SY
        DFZ(IJKPDI(IDI))=DFZ(IJKPDI(IDI))+FI(IJKBDI(IDI))*SZ
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
    !
    call normalArea(IJKP,IJKN,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
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
    FDFII=VSOL*(DFXI*XPN+DFYI*YPN+DFZI*ZPN)
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
    do IDI=1,NDIRA
        IJKP=IJKPDI(IDI)
        IJKB=IJKBDI(IDI)
        !IJK1=IJKDI1(IDI)
        IJK2=IJKDI2(IDI)
        IJK3=IJKDI3(IDI)
        IJK4=IJKDI4(IDI)
        !
        call normalArea(IJKB,IJKP,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
        !
        SX=AR*NX
        SY=AR*NY
        SZ=AR*NZ
        !
        COEFC=RHO*(SX*VX+SY*VY+SZ*VZ)
        COEFD=ALPHA*SRDDI(IDI)
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
