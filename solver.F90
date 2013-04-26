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
    
!
!==========================================================
!....INPUT DATA AND INITIALIZATION
!==========================================================
!
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    !call PetscGetTime(time1,ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

    open(unit=9,FILE='ERR.out')
    rewind 9

    if (rank .eq. 0) print *, 'STARTING SOLVER'
    call init
    print '(I2,A)', rank, ': SETTING UP KSP'
    call setUpKSP
!
!==========================================================
!....START TIME LOOP
!==========================================================
!
    ITIMS=1
    !ITIME=32
    ITIME=1
    print '(I2,A)',rank, ': STARTING TIMELOOP'
    do ITIM=ITIMS,ITIME
        print *, 'TIMESTEP: ', ITIM
!
!.....SHIFT SOLUTIONS IN TIME (OLD = CURRENT)
!
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

       print '(I2,A)', rank, ': UPDATING BOUNDARIES'
       call updateBd
       print '(I2,A)', rank, ': SETTING ANALYTICAL SOLUTION'
       call setSolution
!
!==========================================================
!....START OUTER ITERATIONS
!==========================================================
!
        !LSG=100
        LSG=1
        !LSG=2
        do LS=1,LSG
            if (CONVERGED) then
                !call writeVtk
                exit
            else
                if (rank .eq. 1) print *, 'OUTER ITERATION: ', LS
                print '(I2,A)', rank, ': UPDATING GHOST VALUES'
                call updateGhost
                !print *, '  CALCULATING VELOCITY FIELD'
                !call calcuvw
                print '(I2,A)', rank, ': SOLVING TRANSPORT EQUATION'
                call calcsc
            end if
        end do
        CONVERGED=.false.
        !call setField
    end do
!
!==========================================================
!....FREE WORK SPACE AND FINALIZE PROGRAM
!==========================================================
!
    call cleanUp(A_Mat,B_Vec,SOL_Vec)
    call PetscGetTime(time2,ierr)
    tges=time2-time1
    if (rank .eq. 0) print *, 'TGES ', tges
    if (rank .eq. 0) write(9,*), 'TGES ', tges
    call PetscFinalize(ierr)

end program main

!==========================================================
!>  reads in grid geometry and boundary conditions and
!>  initializes matrices and vectors
!#########################################################
subroutine init
!#########################################################
    
    use boundaryModule
    use coefModule
    use controlModule
    use charModule
    use geoModule
    use indexModule
    use parameterModule
    use scalarModule
    use varModule
    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscksp.h>
#include <finclude/petscpc.h>

    ! MUSS NOCH ANGEPASST WERDEN, ENTWEDER PROCFILE ODER NUR AUS SCALAR MODULE!
    DT=0.0d0
    !DT=1.0d0
    !DT=0.5d0
    !DT=0.25d0
    !DT=0.125d0
    !DT=0.0625d0
    !DT=0.0625d0/3.8258604184154104d0
    !DT=0.03125d0
    !DT=0.015625d0
    !DT=0.0078125d0
    !DT=0.00390625d0
    !if (rank.eq.0) then
        !print *, 'ENTER DT (0 - stationary)'
        !read *, DT
        if (DT.gt.0.0d0) then
            LTIME=.true.
            print *, 'DT: ', DT
            !print *, 'ENTER T_0'
            !read *, T_0
        else
            LTIME=.false.
        end if
    !end if

    PROC=rank
    write(PROC_CH,*) PROC
    PROCUNIT=PROC+PROCOFFSET
    PROCFILE='proc_'//trim(adjustl(PROC_CH))//'.inp'

    open(UNIT=PROCUNIT,FILE=PROCFILE)
    print *, 'opening ... ', PROCFILE
    rewind PROCUNIT 

    read(PROCUNIT,*) NB
    print *, NB
    read(PROCUNIT,*) (B_GLO(B),B=1,NB)
    print *, (B_GLO(B),B=1,NB)

    write(BLOCK_CH,*) B_GLO(1)
    BLOCKUNIT=BLOCKOFFSET+B_GLO(1)
    BLOCKFILE='grid_'//trim(adjustl(BLOCK_CH))//'.out'
    open(UNIT=BLOCKUNIT,FILE=BLOCKFILE)
    print *, 'opening ... ', BLOCKFILE
    rewind BLOCKUNIT
    read(BLOCKUNIT,*) NI,NJ,NK,NIJK,NDIR,NNEU,NWAL,NBLO,NFACE,N,IJKPROC,IJKPROC_GLO
    !print *, NDIR,NNEU
    
    IBL(1)=0
    JBL(1)=0
    KBL(1)=0
    IJKBL(1)=0 
    IJKDIRBL(1)=0
    IJKNEUBL(1)=0
    IJKWALBL(1)=0
    IJKBLOBL(1)=0
    NIBL(1)=NI
    NJBL(1)=NJ
    NKBL(1)=NK
    NIJKBL(1)=NIJK
    NDIRBL(1)=NDIR
    NNEUBL(1)=NNEU
    NWALBL(1)=NWAL
    NBLOBL(1)=NBLO
    NFACEBL(1)=NFACE
    !IJKPROC=IJKST
    !IJKPROC_GLO=IJKPROC
    NBL(1)=N
    !print *, 1,N,IJKPROC

    do B=2,NB
        !read(PROCUNIT,*) B_GLO(B)
        BLOCKUNIT=BLOCKOFFSET+B_GLO(B)
        write(BLOCK_CH,*) B_GLO(B)
        BLOCKFILE='grid_'//trim(adjustl(BLOCK_CH))//'.out'
        !print *, BLOCKFILE
        open(UNIT=BLOCKUNIT,FILE=BLOCKFILE)
        !print *, 'opening ... ', BLOCKFILE
        rewind BLOCKUNIT
        read(BLOCKUNIT,*) NI,NJ,NK,NIJK,NDIR,NNEU,NWAL,NBLO,NFACE,N

        BB=B-1
        IBL(B)=IBL(BB)+NIBL(BB)
        JBL(B)=JBL(BB)+NJBL(BB)
        KBL(B)=KBL(BB)+NKBL(BB)
        IJKBL(B)=IJKBL(BB)+NIJKBL(BB)
        IJKDIRBL(B)=IJKDIRBL(BB)+NDIRBL(BB)
        IJKNEUBL(B)=IJKNEUBL(BB)+NNEUBL(BB)
        IJKWALBL(B)=IJKWALBL(BB)+NWALBL(BB)
        IJKBLOBL(B)=IJKBLOBL(BB)+NBLOBL(BB)
        FACEBL(B)=FACEBL(BB)+NFACEBL(BB)
        NIBL(B)=NI
        NJBL(B)=NJ
        NKBL(B)=NK
        NIJKBL(B)=NIJK
        NDIRBL(B)=NDIR
        NNEUBL(B)=NNEU
        NWALBL(B)=NWAL
        NBLOBL(B)=NBLO
        NFACEBL(B)=NFACE
        !IJKBL_GLO(B)=IJKST
        NBL(B)=N
        !print *, B,N
    end do
    
    ! calculate processor load
    N=sum(NBL)
    print '(I2,A,I10)',rank, ': PROCESSOR LOAD: ', N
    NIJKPROC=sum(NIJKBL)
    NDIRPROC=sum(NDIRBL)
    NNEUPROC=sum(NNEUBL)
    NWALPROC=sum(NWALBL)
    NBLOPROC=sum(NBLOBL)
    NFACEPROC=sum(NFACEBL)
    !print *, NIJKPROC,NDIRPROC,NNEUPROC,NWALPROC,NBLOPROC,NFACEPROC

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

        read(BLOCKUNIT,*) (IJKBBLO(IJKBLOST+IJK),IJK=1,NBLO)
        read(BLOCKUNIT,*) (IJKPBLO(IJKBLOST+IJK),IJK=1,NBLO)
        read(BLOCKUNIT,*) (IJKBLO1(IJKBLOST+IJK),IJK=1,NBLO)
        read(BLOCKUNIT,*) (IJKBLO2(IJKBLOST+IJK),IJK=1,NBLO)
        read(BLOCKUNIT,*) (IJKBLO3(IJKBLOST+IJK),IJK=1,NBLO)
        read(BLOCKUNIT,*) (IJKBLO4(IJKBLOST+IJK),IJK=1,NBLO)

        read(BLOCKUNIT,*) (IJKBDIR(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKPDIR(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKDIR1(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKDIR2(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKDIR3(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*) (IJKDIR4(IJKDIRST+IJK),IJK=1,NDIR)

        read(BLOCKUNIT,*) (IJKBNEU(IJKNEUST+IJK),IJK=1,NNEU)
        read(BLOCKUNIT,*) (IJKPNEU(IJKNEUST+IJK),IJK=1,NNEU)
        read(BLOCKUNIT,*) (IJKNEU1(IJKNEUST+IJK),IJK=1,NNEU)
        read(BLOCKUNIT,*) (IJKNEU2(IJKNEUST+IJK),IJK=1,NNEU)
        read(BLOCKUNIT,*) (IJKNEU3(IJKNEUST+IJK),IJK=1,NNEU)
        read(BLOCKUNIT,*) (IJKNEU4(IJKNEUST+IJK),IJK=1,NNEU)

        read(BLOCKUNIT,*) (IJKBWAL(IJKWALST+IJK),IJK=1,NWAL)
        read(BLOCKUNIT,*) (IJKPWAL(IJKWALST+IJK),IJK=1,NWAL)
        read(BLOCKUNIT,*) (IJKWAL1(IJKWALST+IJK),IJK=1,NWAL)
        read(BLOCKUNIT,*) (IJKWAL2(IJKWALST+IJK),IJK=1,NWAL)
        read(BLOCKUNIT,*) (IJKWAL3(IJKWALST+IJK),IJK=1,NWAL)
        read(BLOCKUNIT,*) (IJKWAL4(IJKWALST+IJK),IJK=1,NWAL)

        read(BLOCKUNIT,*) (FX(IJKST+IJK), IJK=1,NIJK)
        read(BLOCKUNIT,*) (FY(IJKST+IJK), IJK=1,NIJK)
        read(BLOCKUNIT,*) (FZ(IJKST+IJK), IJK=1,NIJK)
        read(BLOCKUNIT,*) DXBL(B),DYBL(B),DZBL(B),VOLBL(B)
        !read(BLOCKUNIT,*) (SRDDIR(IJKDIRST+I),I=1,NDIR)
        !read(BLOCKUNIT,*) (SRDNEU(IJKNEUST+I),I=1,NNEU)
        read(BLOCKUNIT,*) (SRDWAL(IJKWALST+I),I=1,NWAL)
        !
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
    read(BLOCKUNIT,*) (MIJK(IJK),IJK=1,NXYZA*NPROCS)

    do B=1,NB
        BLOCKUNIT=BLOCKOFFSET+B_GLO(B)
        close(UNIT=BLOCKUNIT)
    end do

    call distributeLoad(N)

    TIME=0.0d0
    if (LTIME) then
        call setField
    end if

    print *, 'REMAPPING VALUES'
    do B=1,NB
        call setBlockInd(B)
        do IJKBLO=IJKBLOST+1,IJKBLOST+NBLO
            IJKBBLO(IJKBLO)=IJKBBLO(IJKBLO)-IJKPROC
            IJKPBLO(IJKBLO)=IJKPBLO(IJKBLO)-IJKPROC
            IJKBLO1(IJKBLO)=IJKBLO1(IJKBLO)-IJKPROC
            IJKBLO2(IJKBLO)=IJKBLO2(IJKBLO)-IJKPROC
            IJKBLO3(IJKBLO)=IJKBLO3(IJKBLO)-IJKPROC
            IJKBLO4(IJKBLO)=IJKBLO4(IJKBLO)-IJKPROC
        end do
        !
        do IJKDIR=IJKDIRST+1,IJKDIRST+NDIR
            IJKBDIR(IJKDIR)=IJKBDIR(IJKDIR)-IJKPROC
            IJKPDIR(IJKDIR)=IJKPDIR(IJKDIR)-IJKPROC
            IJKDIR1(IJKDIR)=IJKDIR1(IJKDIR)-IJKPROC
            IJKDIR2(IJKDIR)=IJKDIR2(IJKDIR)-IJKPROC
            IJKDIR3(IJKDIR)=IJKDIR3(IJKDIR)-IJKPROC
            IJKDIR4(IJKDIR)=IJKDIR4(IJKDIR)-IJKPROC
        end do
        !
        do IJKNEU=IJKNEUST+1,IJKNEUST+NNEU
            IJKBNEU(IJKNEU)=IJKBNEU(IJKNEU)-IJKPROC
            IJKPNEU(IJKNEU)=IJKPNEU(IJKNEU)-IJKPROC
            IJKNEU1(IJKNEU)=IJKNEU1(IJKNEU)-IJKPROC
            IJKNEU2(IJKNEU)=IJKNEU2(IJKNEU)-IJKPROC
            IJKNEU3(IJKNEU)=IJKNEU3(IJKNEU)-IJKPROC
            IJKNEU4(IJKNEU)=IJKNEU4(IJKNEU)-IJKPROC
        end do
        !
        do IJKWAL=IJKWALST+1,IJKWALST+NWAL
            IJKBWAL(IJKWAL)=IJKBWAL(IJKWAL)-IJKPROC
            IJKPWAL(IJKWAL)=IJKPWAL(IJKWAL)-IJKPROC
            IJKWAL1(IJKWAL)=IJKWAL1(IJKWAL)-IJKPROC
            IJKWAL2(IJKWAL)=IJKWAL2(IJKWAL)-IJKPROC
            IJKWAL3(IJKWAL)=IJKWAL3(IJKWAL)-IJKPROC
            IJKWAL4(IJKWAL)=IJKWAL4(IJKWAL)-IJKPROC
        end do
    end do

end subroutine init

!========================================================
!>  sets the initial field values if a instationary
!>  problem is solved
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

!========================================================
!>  creates a .vtk containing grid geometry and scalar
!>  field values of the current block
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
        write(TIME_CH,*) ITIM
        write(BLOCK_CH,*) B_GLO(B)
        !write(TIME_CH,'(I1)') LS
        VTKFILE='grid_'//trim(adjustl(BLOCK_CH))//'_'//trim(adjustl(TIME_CH))//'.vtk'
        print *, ' *** GENERATING .VTK *** '
        !
        open (UNIT=BLOCKUNIT,FILE=VTKFILE)
        rewind BLOCKUNIT
        write(BLOCKUNIT,'(A)') '# vtk DataFile Version 3.0'
        write(BLOCKUNIT,'(A)') 'grid'
        write(BLOCKUNIT,'(A)') 'ASCII'
        write(BLOCKUNIT,'(A)') 'DATASET STRUCTURED_GRID'
        write(BLOCKUNIT,'(A,I6,I6,I6)') 'DIMENSIONS', NIM,NJM,NKM
        write(BLOCKUNIT,'(A6,I9,A6)') 'Points', NIM*NJM*NKM, ' float'
        !
        do K=1,NKM
        do I=1,NIM
        do J=1,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            write(BLOCKUNIT,'(E20.10,1X,E20.10,1X,E20.10)'), X(IJK),Y(IJK),Z(IJK)
        end do
        end do
        end do
        !
        write(BLOCKUNIT,'(A10,1X,I9)') 'CELL_DATA ',(NICV*NJCV*NKCV)
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

!=======================================================
!>  updates the boundary field values if a instationary
!>  problem is being solved
!########################################################
subroutine updateBd
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
        do IJKDIR=IJKDIRST+1,IJKDIRST+NDIR
            IJKB=IJKBDIR(IJKDIR)
            T(IJKB)=phi(XC(IJKB),YC(IJKB),ZC(IJKB),TIME)
        end do
        do IJKNEU=IJKNEUST+1,IJKNEUST+NNEU
            IJKB=IJKBNEU(IJKNEU)
            T(IJKB)=phi(XC(IJKB),YC(IJKB),ZC(IJKB),TIME)
        end do
        do IJKWAL=IJKWALST+1,IJKWALST+NWAL
            IJKB=IJKBWAL(IJKWAL)
            T(IJKB)=phi(XC(IJKB),YC(IJKB),ZC(IJKB),TIME)
        end do
    end do

end subroutine updateBd

!=========================================================
!>   gathers all necessary neighbour values needed for
!>   the correct calculation of diffusive and convective
!>   fluxes through block boundary faces
!#########################################################
subroutine updateGhost
!#########################################################

    use boundaryModule
    use coefModule
    use controlModule
    use geoModule
    use indexModule
    use mmsModule
    use scalarModule
    use parameterModule
    use varModule
    implicit none
#include <finclude/petscsys.h>
    integer :: IJK1,IJK2,IJK3,IJK4

    call VecToArr(NFACEAL,MIJK(R),SOL_Vec,TR)

    do B=1,NB
        call setBlockInd(B)

        !...Calculate Mass Fluxes at all(!) cell faces
        !if (rank .eq. 1) print *, 'CALCULATE MASS FLUXES'
        do K=2,NKM
        !do I=2,NIM-1
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            F1(IJK)=RHO*DY*DZ*VX
            !if (rank .eq. 1) print *, F1(IJK)
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

        !print *, 'mass fluxes'
        do IJKDIR=IJKDIRST+1,IJKDIRST+NDIR
            IJKB=IJKBDIR(IJKDIR)
            IJKP=IJKPDIR(IJKDIR)
            !IJK1=IJKDIR1(IJKDIR)
            IJK2=IJKDIR2(IJKDIR)
            IJK3=IJKDIR3(IJKDIR)
            IJK4=IJKDIR4(IJKDIR)
            call normalArea(IJKP,IJKB,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
            FDIR(IJKDIR)=RHO*AR*(VX*NX+VY*NY+VZ*NZ)
        end do

        do IJKNEU=IJKNEUST+1,IJKNEUST+NNEU
            IJKB=IJKBNEU(IJKNEU)
            IJKP=IJKPNEU(IJKNEU)
            !IJK1=IJKNEU1(IJKNEU)
            IJK2=IJKNEU2(IJKNEU)
            IJK3=IJKNEU3(IJKNEU)
            IJK4=IJKNEU4(IJKNEU)
            call normalArea(IJKP,IJKB,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
            !print *,MIJK(IJKP),AR,NX,NY,NZ
            FNEU(IJKNEU)=RHO*AR*(VX*NX+VY*NY+VZ*NZ)
        end do
    end do
    
end subroutine updateGhost

!=============================================================
!>     discretizes (FVM) and solves the scalar transport
!>     equation.
!#############################################################
subroutine calcSc
!#############################################################

    use boundaryModule
    use coefModule
    use controlModule
    use geoModule
    use gradModule
    use indexModule
    use mmsModule
    use petsc_ksp_module
    use scalarModule
    use varModule
    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>

    real*8 :: APT,URF,CB,CP
    integer :: IJK1,IJK2,IJK3,IJK4
    character(len=8) :: NOTHING

    URF=1.0d0

    call gradfi(T,TR,DTX,DTY,DTZ,DTX_Vec,DTY_Vec,DTZ_Vec)
    call updateGrad
!
!.....START BLOCK LOOP
!
    do B=1,NB
        call setBlockInd(B)
!
!.....INITIALIZE Q AND AP
!

        do K=2,NKM
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            Q(IJK)=src(XC(IJK), YC(IJK), ZC(IJK), TIME)*VOL
            AP(IJK)=0.0d0
        end do
        end do
        end do
!
!.....FLUXES THROUGH EAST FACE
!
        do K=2,NKM
        do I=2,NIM-1
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            call fluxsc(IJK,IJK+NJ,IJK-NIJ-1,IJK-1,IJK,F1(IJK),AW(IJK+NJ),AE(IJK),FX(IJK),1.0d0)
        end do
        end do
        end do
!
!.....FLUXES THROUGH NORTH FACE
!
        do K=2,NKM
        do I=2,NIM
        do J=2,NJM-1
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            call fluxsc(IJK,IJK+1,IJK-NIJ,IJK,IJK-NJ,F2(IJK),AS(IJK+1),AN(IJK),FY(IJK),1.0d0)
        end do
        end do
        end do
!
!.....FLUXES THROUGH TOP FACE
!
        do K=2,NKM-1
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            call fluxsc(IJK,IJK+NIJ,IJK-NJ-1,IJK-NJ,IJK,F3(IJK),AB(IJK+NIJ),AT(IJK),FZ(IJK),1.0d0)
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
!.....DIRICHLET BOUNARIES
!
        !print *, 'DIRICHLET BOUNDARIES'
        do IJKDIR=IJKDIRST+1,IJKDIRST+NDIR
            IJKP=IJKPDIR(IJKDIR)
            IJKB=IJKBDIR(IJKDIR)
            !IJK1=IJKDIR1(IJKDIR)
            IJK2=IJKDIR2(IJKDIR)
            IJK3=IJKDIR3(IJKDIR)
            IJK4=IJKDIR4(IJKDIR)
            DTX(IJKB)=DTX(IJKP)
            DTY(IJKB)=DTY(IJKP)
            DTZ(IJKB)=DTZ(IJKP)

            call fluxsc(IJKP,IJKB,IJK2,IJK3,IJK4,FDIR(IJKDIR),CP,CB,ONE,ZERO)

            AP(IJKP)=AP(IJKP)-CB
            Q(IJKP)=Q(IJKP)-CB*T(IJKB)
        end do
!
!.....NEUMANN ZERO GRADIENT BOUNARIES
!
        !print *, 'NEUMANN BOUNDARIES'
        do IJKNEU=IJKNEUST+1,IJKNEUST+NNEU
            IJKP=IJKPNEU(IJKNEU)
            IJKB=IJKBNEU(IJKNEU)
            T(IJKB)=T(IJKP)
            !IJK1=IJKNEU1(IJKNEU)
            IJK2=IJKNEU2(IJKNEU)
            IJK3=IJKNEU3(IJKNEU)
            IJK4=IJKNEU4(IJKNEU)
            DTX(IJKB)=DTX(IJKP)
            DTY(IJKB)=DTY(IJKP)
            DTZ(IJKB)=DTZ(IJKP)

            call fluxsc(IJKP,IJKB,IJK2,IJK3,IJK4,FNEU(IJKNEU),CP,CB,ONE,ZERO)

        end do
!
!.....WALL BOUNDARY CONTRIBUTION
!
        call walBdFlux
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
!        
!.....ASSEMBLE MATRIX AND RHS-VECTOR
!
        do K=3,NKM-1
        do I=3,NIM-1
        do J=3,NJM-1
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            IJKP=MIJK(IJK+IJKPROC)
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

            !print *,IJKPROC, row, col
            !
            call MatSetValues(A_Mat,i1,row,i7,col,val,INSERT_VALUES,ierr)
            call VecSetValue(B_Vec,row,valq,INSERT_VALUES,ierr)
            !
        end do
        end do
        end do
        !stop

        ! Assembly matrix coefficients of inlet boundaries
        !print *, 'Assembling inlet boundary terms'
        do IJKDIR=IJKDIRST+1,IJKDIRST+NDIR
            IJK=IJKPDIR(IJKDIR)
            IJKP=MIJK(IJK+IJKPROC)
            !print *, IJKP
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
            !print *,IJKPROC,MIJK(IJK),col,val
            !
            call MatSetValues(A_Mat,i1,row,i7,col,val,INSERT_VALUES,ierr)
            call VecSetValue(B_Vec,row,valq,INSERT_VALUES,ierr)
            !
            !print *, row
            if (row .eq. 276) then
                !print *, 'Matrix coefficients to set - Iteration No.', LS,IJK,IJKP, val
            end if

        end do

        ! Assembly matrix coefficients of neumann boundaries
        do IJKNEU=IJKNEUST+1,IJKNEUST+NNEU
            IJK=IJKPNEU(IJKNEU)
            IJKP=MIJK(IJK+IJKPROC)
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
            call MatSetValues(A_Mat,i1,row,i7,col,val,INSERT_VALUES,ierr)
            call VecSetValue(B_Vec,row,valq,INSERT_VALUES,ierr)
            !
        end do

        ! Assembly matrix coefficients of wall boundaries
        do IJKWAL=IJKWALST+1,IJKWALST+NWAL
            IJK=IJKPWAL(IJKWAL)
            IJKP=MIJK(IJK+IJKPROC)
            row=IJKP
            col=(/-1,-1,-1,IJKP,-1,-1,-1/)
            !print *, IJKWAL,IJKPDI(IJKWAL),row
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
            call MatSetValues(A_Mat,i1,row,i7,col,val,INSERT_VALUES,ierr)
            call VecSetValue(B_Vec,row,valq,INSERT_VALUES,ierr)
            !
        end do

        ! Assembly matrix coefficients of block boundary cells
        !print *, 'Assembly block boundary coefficients'
        do IJKBLO=IJKBLOST+1,IJKBLOST+NBLO
            IJK=IJKPBLO(IJKBLO)
            IJKP=MIJK(IJK+IJKPROC)
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
            !print *,IJKPROC,row,col
            call MatSetValues(A_Mat,i1,row,i7,col,val,INSERT_VALUES,ierr)
            call VecSetValue(B_Vec,row,valq,INSERT_VALUES,ierr)
            !
        end do

        ! Assembly off diagonal matrix coefficients of block boundaries
        !print *, 'FACES'
        do F=FACEST+1,FACEST+NFACE
            row=MIJK(L(F))
            col1=MIJK(R(F))
            !print *, IJKPROC,row,col1
            val1=AF(F)
            !call MatSetValues(A_Mat,i1,col1,val1,INSERT_VALUES,ierr)
            call MatSetValue(A_Mat,row,col1,val1,INSERT_VALUES,ierr)
        end do
    end do

    ! Assembly matrix and right hand vector
    print '(I2,A)', rank,  ': STARTING MATRIX ASSEMBLY'
    call MatAssemblyBegin(A_Mat,MAT_FINAL_ASSEMBLY,ierr)
    call VecAssemblyBegin(B_Vec,ierr)
    call MatAssemblyEnd(A_Mat,MAT_FINAL_ASSEMBLY,ierr)
    call VecAssemblyEnd(B_Vec,ierr)

    !stop
    ! muss noch implementiert werden
    !call MatGetInfo(A_Mat,MAT_LOCAL,ierr)
    !stop
    !call MatView(A_Mat,PETSC_VIEWER_STDOUT_WORLD,ierr)
    !print *, ''
    !call VecView(B_Vec,PETSC_VIEWER_STDOUT_WORLD,ierr)
    !call VecView(SOL_Vec,PETSC_VIEWER_STDOUT_WORLD,ierr)
!
!.....SOLVE LINEAR SYSTEM
!
    print '(I2,A)', rank, ': SOLVING LINEAR SYSTEM'
    !call PetscGetCPUTime(time1,ierr)
    call solveSys(A_Mat,B_Vec,SOL_Vec,N,LS,res_Scalar)

    if (CONVERGED) then
        return
    end if

    !tges=time2-time1

    do B=1,NB
        call setBlockInd(B)
        do K=2,NKM
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            row=MIJK(IJK+IJKPROC)
            call VecGetValues(SOL_Vec,i1,row,valt,ierr)
            T(IJK)=valt
        end do
        end do
        end do
    end do

    !stop
    call calcErr

end subroutine calcSc

!===============================================================
!>  calculates the components of the gradient vector of a
!>  scalar FI at the CV center, using conservative scheme based
!>  on the Gauss theorem; FIE are values at east side, FIN at
!>  north side, FIT at top side. Contributions from boundary
!>  faces are calculated in separate loops. s.a. updateGrad
!################################################################
subroutine gradfi(FI,FIR,DFX,DFY,DFZ,DFX_vec,DFY_vec,DFZ_vec)
!################################################################

    use boundaryModule
    use controlModule
    use geoModule
    use gradModule
    use indexModule
    use parameterModule
    use scalarModule
    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>

    real(KIND=PREC), intent(in out) :: FI(NXYZA),FIR(NFACEAL),DFX(NXYZA),DFY(NXYZA),DFZ(NXYZA)
    Vec , intent(in out) :: DFX_vec,DFY_vec,DFZ_vec
    integer :: IJK1, IJK2, IJK3, IJK4
    PetscScalar, pointer :: DFX_Sca(:),DFY_Sca(:),DFZ_Sca(:)

    ! modify the local portion of the vector
    call VecGetArrayF90(DFX_Vec,DFX_Sca,ierr)
    call VecGetArrayF90(DFY_Vec,DFY_Sca,ierr)
    call VecGetArrayF90(DFZ_Vec,DFZ_Sca,ierr)

    do B=1,NB
        call setBlockInd(B)
!
!.....INITIALIZE FIELDS
!
        do IJK=IJKST+1,IJKST+NIJK
            DFX(IJK)=0.0d0
            DFY(IJK)=0.0d0
            DFZ(IJK)=0.0d0
        end do
!
!.....CONTRRIBUTION FROM INNER EAST SIDES
!
        do K=2,NKM
        do I=2,NIM-1
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            FIE=FI(IJK+NJ)*FX(IJK)+FI(IJK)*(1.0d0-FX(IJK))

            IJK4=IJK
            IJK3=IJK4-1
            IJK2=IJK3-NIJ
            IJK1=IJK4-NIJ
            
            call normalArea(IJK,IJK,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
            
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
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            FIN=FI(IJK+1)*FY(IJK)+FI(IJK)*(1.0d0-FY(IJK))

            IJK3=IJK
            IJK4=IJK3-NJ
            IJK2=IJK3-NIJ
            IJK1=IJK4-NIJ
            
            call normalArea(IJK,IJK,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
            
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
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            FIT=FI(IJK+NIJ)*FZ(IJK)+FI(IJK)*(1.0d0-FZ(IJK))

            IJK4=IJK
            IJK3=IJK4-NJ
            IJK1=IJK4-1
            IJK2=IJK3-1
            
            call normalArea(IJK,IJK,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
            
            SX=AR*NX
            SY=AR*NY
            SZ=AR*NZ

            DFXT=FIT*SX
            DFYT=FIT*SY
            DFZT=FIT*SZ

            DFX(IJK)=DFX(IJK)+DFXT
            DFY(IJK)=DFY(IJK)+DFYT
            DFZ(IJK)=DFZ(IJK)+DFZT
            DFX(IJK+1)=DFX(IJK+1)-DFXT
            DFY(IJK+1)=DFY(IJK+1)-DFYT
            DFZ(IJK+1)=DFZ(IJK+1)-DFZT
        end do
        end do
        end do
!
!.....CONTRIBUTION FROM DIRICHLET BOUNDARIES
!
        do IJKDIR=IJKDIRST+1,IJKDIRST+NDIR
            IJKB=IJKBDIR(IJKDIR)
            IJKP=IJKPDIR(IJKDIR)
            !IJK1=IJKDIR1(IJKDIR)
            IJK2=IJKDIR2(IJKDIR)
            IJK3=IJKDIR3(IJKDIR)
            IJK4=IJKDIR4(IJKDIR)
            
            call normalArea(IJKP,IJKB,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
            
            SX=AR*NX
            SY=AR*NY
            SZ=AR*NZ

            DFX(IJKP)=DFX(IJKP)+FI(IJKB)*SX
            DFY(IJKP)=DFY(IJKP)+FI(IJKB)*SY
            DFZ(IJKP)=DFZ(IJKP)+FI(IJKB)*SZ
        end do
!
!.....CONTRIBUTION FROM NEUMANN BOUNDARIES
!
        do IJKNEU=IJKNEUST+1,IJKNEUST+NNEU
            IJKB=IJKBNEU(IJKNEU)
            IJKP=IJKPNEU(IJKNEU)
            !IJK1=IJKNEU1(IJKNEU)
            IJK2=IJKNEU2(IJKNEU)
            IJK3=IJKNEU3(IJKNEU)
            IJK4=IJKNEU4(IJKNEU)
            
            call normalArea(IJKP,IJKB,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
            
            SX=AR*NX
            SY=AR*NY
            SZ=AR*NZ

            DFX(IJKPNEU(IJKNEU))=DFX(IJKPNEU(IJKNEU))+FI(IJKBNEU(IJKNEU))*SX
            DFY(IJKPNEU(IJKNEU))=DFY(IJKPNEU(IJKNEU))+FI(IJKBNEU(IJKNEU))*SY
            DFZ(IJKPNEU(IJKNEU))=DFZ(IJKPNEU(IJKNEU))+FI(IJKBNEU(IJKNEU))*SZ
        end do
!
!.....CONTRIBUTION FROM WALL BOUNDARIES
!
        do IJKWAL=IJKWALST+1,IJKWALST+NWAL
            IJKB=IJKBWAL(IJKWAL)
            IJKP=IJKPWAL(IJKWAL)
            !IJK1=IJKWAL1(IJKWAL)
            IJK2=IJKWAL2(IJKWAL)
            IJK3=IJKWAL3(IJKWAL)
            IJK4=IJKWAL4(IJKWAL)
            
            call normalArea(IJKP,IJKB,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
            
            SX=AR*NX
            SY=AR*NY
            SZ=AR*NZ

            DFX(IJKPWAL(IJKWAL))=DFX(IJKPWAL(IJKWAL))+FI(IJKBWAL(IJKWAL))*SX
            DFY(IJKPWAL(IJKWAL))=DFY(IJKPWAL(IJKWAL))+FI(IJKBWAL(IJKWAL))*SY
            DFZ(IJKPWAL(IJKWAL))=DFZ(IJKPWAL(IJKWAL))+FI(IJKBWAL(IJKWAL))*SZ
        end do
!
!.....CONTRIBUTION FROM BLOCK BOUNDARIES
!
        do F=FACEST+1,FACEST+NFACE
            FIF=FIR(F)*FF(F)+FI(L(F)-IJKPROC)*(1.0d0-FF(F))
            
            SX=ARF(F)*NXF(F)
            SY=ARF(F)*NYF(F)
            SZ=ARF(F)*NZF(F)
            
            DFXF=FIF*SX
            DFYF=FIF*SY
            DFZF=FIF*SZ
            
            DFX(L(F)-IJKPROC)=DFX(L(F)-IJKPROC)+DFXF
            DFY(L(F)-IJKPROC)=DFY(L(F)-IJKPROC)+DFYF
            DFZ(L(F)-IJKPROC)=DFZ(L(F)-IJKPROC)+DFZF
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

            DFX_Sca(MIJK(IJK+IJKPROC)-IJKPROC_GLO+1)=DFX(IJK)
            DFY_Sca(MIJK(IJK+IJKPROC)-IJKPROC_GLO+1)=DFY(IJK)
            DFZ_Sca(MIJK(IJK+IJKPROC)-IJKPROC_GLO+1)=DFZ(IJK)
        end do
        end do
        end do
    end do

    call VecRestoreArrayF90(DFX_Vec,DFX_Sca,ierr)
    call VecRestoreArrayF90(DFY_Vec,DFY_Sca,ierr)
    call VecRestoreArrayF90(DFZ_Vec,DFZ_Sca,ierr)

end subroutine gradfi

!===============================================================
!>  scatters the block neighbour components of the gradient
!>  vector of the scalar  FI at the CV center (both local
!>  and off-processor).
!################################################################
subroutine updateGrad
!################################################################

    use boundaryModule
    use indexModule
    use parameterModule
    use varModule
    implicit none

    call VecToArr(NFACEAL,MIJK(R),DTX_Vec,DTXR)
    call VecToArr(NFACEAL,MIJK(R),DTY_Vec,DTYR)
    call VecToArr(NFACEAL,MIJK(R),DTZ_Vec,DTZR)

end subroutine updateGrad

!==============================================================
!>   calculates fluxes (convective and diffusive) through the
!>   inner cell face between nodes IJKP and IJKN. IJK1...4 are the
!>   indices of CV corners defining the cell face. FM is the
!>   mass flux through the face, and FAC is the interpolation
!>   factor (distance from node IJKP to cell face center over
!>   the sum of this distance and the distance from cell face
!>   center to node IJN). CAP and CAN are the contributions to
!>   matrix coefficients in the transport equation at nodes IJKP
!>   and IJKN. Diffusive fluxes are discretized using central
!>   differences; for convective fluxes, linear interpolation
!>   can be blended (G) with upwind approximation. Note: cell
!>   face surface vector is directed from P to N.
!################################################################
subroutine fluxSc(IJKP,IJKN,IJK2,IJK3,IJK4,FM,CAP,CAN,FAC,G)
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
    real(KIND=PREC), intent(in) :: FM,FAC,G
    integer, intent(in) :: IJKP,IJKN,IJK2,IJK3,IJK4
    real(KIND=PREC), intent(out) :: CAN, CAP

    FACP=1.0d0-FAC
    FII=T(IJKN)*FAC+T(IJKP)*FACP
    DFXI=DTX(IJKN)*FAC+DTX(IJKP)*FACP
    DFYI=DTY(IJKN)*FAC+DTY(IJKP)*FACP
    DFZI=DTZ(IJKN)*FAC+DTZ(IJKP)*FACP    
!
!.....SURFACE AND DISTANCE VECTOR COMPONENTS, DIFFUSION COEFF.
!
    call normalArea(IJKP,IJKN,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)

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
    !print *, MIJK(IJKP),VSOL,FM
    !if (IJKN.eq.86) print *,'MASS FLUX: ', FM,VSOL 

end subroutine fluxsc

!===============================================================
!>  assembles the source terms (volume integrals) and applies
!>  wall boundary conditions (velocity equals 0) for the
!>  scalar transport equation.
!################################################################
subroutine walBdFlux
!################################################################

    use boundaryModule
    use coefModule
    use geoModule
    use indexModule
    use scalarModule
    use parameterModule
    use varModule
    implicit none

    real(KIND=PREC) :: COEFC,COEFD
    integer :: IJK1, IJK2, IJK3, IJK4

!
!.....WALL BOUNDARY CONDITION
!
    do IJKWAL=IJKWALST+1,IJKWALST+NWAL
        IJKP=IJKPWAL(IJKWAL)
        IJKB=IJKBWAL(IJKWAL)
        !IJK1=IJKWAL1(IJKWAL)
        IJK2=IJKWAL2(IJKWAL)
        IJK3=IJKWAL3(IJKWAL)
        IJK4=IJKWAL4(IJKWAL)
        !
        call normalArea(IJKB,IJKP,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
        !
        SX=AR*NX
        SY=AR*NY
        SZ=AR*NZ
        !
        !COEFC=RHO*(SX*VX+SY*VY+SZ*VZ)
        COEFD=ALPHA*SRDWAL(IJKWAL)
        AP(IJKP)=AP(IJKP)+COEFD
        Q(IJKP)=Q(IJKP)+COEFD*T(IJKB)
        !print *, MIJK(IJKP),AP(IJKP),Q(IJKP),COEFD,COEFC
    end do

end subroutine walBdFlux

!===============================================================
!>  assembles the source terms (volume integrals), and calculates
!>  the convective and diffusive fluxes resulting from the
!>  the proximity relationship of the current block to its 
!>  neighbours
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
    real(KIND=PREC) :: FAC,FM,G

    G=1.0d0

    do F=FACEST+1,FACEST+NFACE
        FAC=FF(F)
        FACP=1.0d0-FAC
        FII=TR(F)*FAC+T(L(F)-IJKPROC)*FACP
        DFXI=DTXR(F)*FAC+DTX(L(F)-IJKPROC)*FACP
        DFYI=DTYR(F)*FAC+DTY(L(F)-IJKPROC)*FACP
        DFZI=DTZR(F)*FAC+DTZ(L(F)-IJKPROC)*FACP

        FM=F1(L(F)-IJKPROC)*NXF(F)+F2(L(F)-IJKPROC)*NYF(F)+F3(L(F)-IJKPROC)*NZF(F)

        VSOL=ALPHA*ARF(F)/(DNF(F)+SMALL)
!
!.....EXPLICIT CONVECTIVE AND DIFFUSIVE FLUXES
!
        FCFIE=FM*FII
        FDFIE=ALPHA*(DFXI*ARF(F)*NXF(F)+DFYI*ARF(F)*NYF(F)+DFZI*ARF(F)*NZF(F))
!
!.....IMPLICIT CONVECTIVE AND DIFFUSIVE FLUXES
!
        FCFII=MIN(FM,ZERO)*TR(F)+MAX(FM,ZERO)*T(L(F)-IJKPROC)
        FDFII=VSOL*(DFXI*XPNF(F)*NXF(F)+DFYI*YPNF(F)*NYF(F)+DFZI*ZPNF(F)*NZF(F))
        !print *, MIJK(L(F)),VSOL,FM
!
!.....COEFFICIENTS, DEFERRED CORRECTION, SOURCE TERMS
!
        ! Boundary must be treated inoutflow like
        AF(F)=-VSOL+MIN(FM,ZERO)
        !print *, AP(L(F))
        AP(L(F)-IJKPROC)=AP(L(F)-IJKPROC)-AF(F)
        !print *, MIJK(L(F)),MIJK(R(F)),MIN(FM,ZERO),AF(F),AP(L(F))
        FFIC=G*(FCFIE-FCFII)
        Q(L(F)-IJKPROC)=Q(L(F)-IJKPROC)-FFIC+FDFIE-FDFII
    end do

end subroutine blockBdFlux

!================================================================
!> calculates the global mean error and outputs error, time and
!> iterations
!################################################################
subroutine calcErr
!################################################################

    use parameterModule
    use coefModule
    use controlModule
    use geoModule
    use indexModule
    use mmsModule
    use varModule
    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
!#include <finclude/petscvec.h90>

    real(KIND=PREC) :: E,ER
    PetscScalar :: ERR_Sca
    PetscInt :: N_GLO

    call VecWAXPY(ERR_Vec,-1.0d0,MMS_Vec,SOL_Vec,ierr)
    call VecNorm(ERR_Vec,NORM_1,ERR_Sca,ierr)
    call VecGetSize(ERR_Vec,N_GLO,ierr)

    !if (rank.eq. 0) print *,'ERROR ',ERR_Sca/N_GLO,'TGES ',tges,'ITS ',itsInt 
    if (rank.eq. 0) print *,'ERROR ',ERR_Sca/N_GLO
    if (rank.eq. 0) write(9,*),'ERROR ',ERR_Sca/N_GLO, 'RESNORM ', res_Scalar, 'ITS ', itsInt

end subroutine calcErr

!================================================================
!>  sets the manufactured solution for the current timestep
!################################################################
subroutine setSolution
!################################################################
    use parameterModule
    use controlModule
    use geoModule
    use indexModule
    use mmsModule
    use varModule
    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>

    PetscScalar, pointer :: TEMP_Sca(:)

    call VecGetArrayF90(MMS_Vec,TEMP_Sca,ierr)

    IJKP=0
    do B=1,NB
        call setBlockInd(B)
        do K=2,NKM
        do I=2,NIM
        do J=2,NJM
            IJKP=IJKP+1
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            TEMP_Sca(IJKP)=phi(XC(IJK),YC(IJK),ZC(IJK),TIME)
        end do
        end do
        end do
    end do

    call VecRestoreArrayF90(MMS_Vec,TEMP_Sca,ierr)

end subroutine setSolution
