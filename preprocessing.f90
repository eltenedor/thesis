!#########################################################
program main
!#########################################################

    use iso_c_binding
    implicit none
    real :: T1,T2

    interface cpp_interface

        subroutine commonFace(XYZL, XYZR, XYZCommon, AR) bind (C, name = "commonFace")

        use iso_c_binding
        implicit none

        real(kind = c_double),intent(inout),dimension(*) :: XYZL,XYZR,XYZCommon
        real(kind = c_double),intent(inout) :: AR

        end subroutine commonFace
    end interface cpp_interface

    call readData

    call cpu_time(T1)
    call findNeighbours
    call cpu_time(T2)
    print *, 'T2-T1 = ', T2-T1

    print *, 'CALCULATING MAXIMAL PROCESSOR LOAD'
    call calcProcessorLoad
    print *, 'WRITING parameterModule.f90'
    call writeParameterModule

end program main

!##########################################################
subroutine readData
!#########################################################

    use boundaryModule
    use charModule
    use geoModule
    use parameterModule
    use indexModule
    implicit none
    
    print *, 'ENTER NUMBER OF PROCS: '
    read *, NPROCSA
    
    NB=NBLOCKS
    BLOCKUNIT=OFFSET+1
    write(BLOCK_CH,*) (BLOCKUNIT-OFFSET)
    BLOCKFILE='grid_'//trim(adjustl(BLOCK_CH))//'.pre'
    
    open(UNIT=BLOCKUNIT,FILE=BLOCKFILE,FORM=FORM_CH,POSITION='REWIND')
    print *, 'READING FROM FILE: ', BLOCKFILE
    read(BLOCKUNIT) NI,NJ,NK,NIJK,NDIR,NNEU,NWAL,NBLO
    
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
    
    do B=2,NB
        BLOCKUNIT=OFFSET+B
        write(BLOCK_CH,*) (BLOCKUNIT-OFFSET)
        BLOCKFILE='grid_'//trim(adjustl(BLOCK_CH))//'.pre'
        open(UNIT=BLOCKUNIT,FILE=BLOCKFILE,FORM=FORM_CH,POSITION='REWIND')
        print *, 'READING FROM FILE: ', BLOCKFILE
        read(BLOCKUNIT) NI,NJ,NK,NIJK,NDIR,NNEU,NWAL,NBLO

        BB=B-1
        IBL(B)=IBL(BB)+NIBL(BB)
        JBL(B)=JBL(BB)+NJBL(BB)
        KBL(B)=KBL(BB)+NKBL(BB)
        IJKBL(B)=IJKBL(BB)+NIJKBL(BB)
        IJKDIRBL(B)=IJKDIRBL(BB)+NDIRBL(BB)
        IJKNEUBL(B)=IJKNEUBL(BB)+NNEUBL(BB)
        IJKWALBL(B)=IJKWALBL(BB)+NWALBL(BB)
        IJKBLOBL(B)=IJKBLOBL(BB)+NBLOBL(BB)
        NIBL(B)=NI
        NJBL(B)=NJ
        NKBL(B)=NK
        NIJKBL(B)=NIJK
        NDIRBL(B)=NDIR
        NNEUBL(B)=NNEU
        NWALBL(B)=NWAL
        NBLOBL(B)=NBLO
    end do

    call createMapping

    do B=1,NB
        call setBlockInd(B)
        BLOCKUNIT=OFFSET+B
        read(BLOCKUNIT)   (NEIGH(B,I),I=1,6)
        !print *, (NEIGH(B,I),I=1,6)
        !read(BLOCKUNIT)   (LK(KST+K),K=1,NK)
        !read(BLOCKUNIT)   (LI(IST+I),I=1,NI)
        
        read(BLOCKUNIT)   (X(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT)   (Y(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT)   (Z(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT)   (XC(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT)   (YC(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT)   (ZC(IJKST+IJK),IJK=1,NIJK)
        
        read(BLOCKUNIT)   (IJKBDIR(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT)   (IJKPDIR(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT)   (IJKDIR1(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT)   (IJKDIR2(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT)   (IJKDIR3(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT)   (IJKDIR4(IJKDIRST+IJK),IJK=1,NDIR)
        
        read(BLOCKUNIT)   (IJKBNEU(IJKNEUST+IJK),IJK=1,NNEU)
        read(BLOCKUNIT)   (IJKPNEU(IJKNEUST+IJK),IJK=1,NNEU)
        read(BLOCKUNIT)   (IJKNEU1(IJKNEUST+IJK),IJK=1,NNEU)
        read(BLOCKUNIT)   (IJKNEU2(IJKNEUST+IJK),IJK=1,NNEU)
        read(BLOCKUNIT)   (IJKNEU3(IJKNEUST+IJK),IJK=1,NNEU)
        read(BLOCKUNIT)   (IJKNEU4(IJKNEUST+IJK),IJK=1,NNEU)
        
        read(BLOCKUNIT)   (IJKBWAL(IJKWALST+IJK),IJK=1,NWAL)
        read(BLOCKUNIT)   (IJKPWAL(IJKWALST+IJK),IJK=1,NWAL)
        read(BLOCKUNIT)   (IJKWAL1(IJKWALST+IJK),IJK=1,NWAL)
        read(BLOCKUNIT)   (IJKWAL2(IJKWALST+IJK),IJK=1,NWAL)
        read(BLOCKUNIT)   (IJKWAL3(IJKWALST+IJK),IJK=1,NWAL)
        read(BLOCKUNIT)   (IJKWAL4(IJKWALST+IJK),IJK=1,NWAL)

        read(BLOCKUNIT)   (IJKBBLO(IJKBLOST+IJK),IJK=1,NBLO)
        read(BLOCKUNIT)   (IJKPBLO(IJKBLOST+IJK),IJK=1,NBLO)
        read(BLOCKUNIT)   (IJKBLO1(IJKBLOST+IJK),IJK=1,NBLO)
        read(BLOCKUNIT)   (IJKBLO2(IJKBLOST+IJK),IJK=1,NBLO)
        read(BLOCKUNIT)   (IJKBLO3(IJKBLOST+IJK),IJK=1,NBLO)
        read(BLOCKUNIT)   (IJKBLO4(IJKBLOST+IJK),IJK=1,NBLO)
        
    end do

    print *, 'REMAPPING VALUES'
    do B=1,NB
        call setBlockInd(B)
        do IJKBLO=IJKBLOST+1,IJKBLOST+NBLO
            IJKBBLO(IJKBLO)=IJKBBLO(IJKBLO)+IJKST
            IJKPBLO(IJKBLO)=IJKPBLO(IJKBLO)+IJKST
            IJKBLO1(IJKBLO)=IJKBLO1(IJKBLO)+IJKST
            IJKBLO2(IJKBLO)=IJKBLO2(IJKBLO)+IJKST
            IJKBLO3(IJKBLO)=IJKBLO3(IJKBLO)+IJKST
            IJKBLO4(IJKBLO)=IJKBLO4(IJKBLO)+IJKST
        end do
        !
        do IJKDIR=IJKDIRST+1,IJKDIRST+NDIR
            IJKBDIR(IJKDIR)=IJKBDIR(IJKDIR)+IJKST
            IJKPDIR(IJKDIR)=IJKPDIR(IJKDIR)+IJKST
            IJKDIR1(IJKDIR)=IJKDIR1(IJKDIR)+IJKST
            IJKDIR2(IJKDIR)=IJKDIR2(IJKDIR)+IJKST
            IJKDIR3(IJKDIR)=IJKDIR3(IJKDIR)+IJKST
            IJKDIR4(IJKDIR)=IJKDIR4(IJKDIR)+IJKST
        end do
        !
        do IJKNEU=IJKNEUST+1,IJKNEUST+NNEU
            IJKBNEU(IJKNEU)=IJKBNEU(IJKNEU)+IJKST
            IJKPNEU(IJKNEU)=IJKPNEU(IJKNEU)+IJKST
            IJKNEU1(IJKNEU)=IJKNEU1(IJKNEU)+IJKST
            IJKNEU2(IJKNEU)=IJKNEU2(IJKNEU)+IJKST
            IJKNEU3(IJKNEU)=IJKNEU3(IJKNEU)+IJKST
            IJKNEU4(IJKNEU)=IJKNEU4(IJKNEU)+IJKST
        end do
        !
        do IJKWAL=IJKWALST+1,IJKWALST+NWAL
            IJKBWAL(IJKWAL)=IJKBWAL(IJKWAL)+IJKST
            IJKPWAL(IJKWAL)=IJKPWAL(IJKWAL)+IJKST
            IJKWAL1(IJKWAL)=IJKWAL1(IJKWAL)+IJKST
            IJKWAL2(IJKWAL)=IJKWAL2(IJKWAL)+IJKST
            IJKWAL3(IJKWAL)=IJKWAL3(IJKWAL)+IJKST
            IJKWAL4(IJKWAL)=IJKWAL4(IJKWAL)+IJKST
        end do
    end do

end subroutine readData

!#########################################################
subroutine findNeighbours
!#########################################################

    use iso_c_binding
    use boundaryModule
    use geoModule
    use controlModule
    use indexModule
    implicit none
    integer :: iterationsCounter
    logical :: equalSize

    interface cpp_interface

        subroutine commonFace(XYZL,XYZR,XYZCommon,AR,XYZF) bind (C, name = "commonFace")

        use iso_c_binding
        implicit none

        real(kind = c_double),intent(inout),dimension(*) :: XYZL,XYZR,XYZCommon,XYZF
        real(kind = c_double),intent(inout) :: AR

        end subroutine commonFace
    end interface cpp_interface

    NF=0
    iterationsCounter=0

    ! use other block loop: B=1,NB
    do B=1,NB
        FACEBL(B)=NF
        call setBlockInd(B)
        IJKSTL=IJKBLOST+1
        IJKENL=IJKBLOST+NBLO
        print *, 'BLOCK: ',B
        neighbour: do INEIGH=1,6
            if (NEIGH(B,INEIGH).gt.0) then
                call setBlockInd(B,NEIGH(B,INEIGH))
                if (NIJKL.eq.NIJKR) then ! this works only with quadratic grids
                    equalSize=.true.
                    !equalSize=.false.
                else
                    equalSize=.false.
                end if
                IJKSTR=IJKBLOSTR+1
                IJKENR=IJKBLOSTR+NBLOR
                STARTED=.false.
                FOUND=.false.
                select case (INEIGH)
                    case (1)
                        !
                        !..........SOUTH..........
                        !
                        print *, 'SOUTH'
                        SouthOuter: do IJKL=IJKSTL,IJKENL
                            !
                            XYZL(1:3)=[X(IJKBLO2(IJKL)),Y(IJKBLO2(IJKL)),Z(IJKBLO2(IJKL))]
                            XYZL(4:6)=[X(IJKBLO1(IJKL)),Y(IJKBLO1(IJKL)),Z(IJKBLO1(IJKL))]
                            XYZL(7:9)=[X(IJKBLO3(IJKL)),Y(IJKBLO3(IJKL)),Z(IJKBLO3(IJKL))]
                            XYZL(10:12)=[X(IJKBLO4(IJKL)),Y(IJKBLO4(IJKL)),Z(IJKBLO4(IJKL))]
                            !
                            SouthInner: do IJKR=IJKSTR,IJKENR
                                iterationsCounter=iterationsCounter+1
                                !print *, iterationsCounter
                                !print *, IJKPBLO(IJKL),IJKPBLO(IJKR)
                                !
                                XYZR(1:3)=[X(IJKBLO1(IJKR)),Y(IJKBLO1(IJKR)),Z(IJKBLO1(IJKR))]
                                XYZR(4:6)=[X(IJKBLO2(IJKR)),Y(IJKBLO2(IJKR)),Z(IJKBLO2(IJKR))]
                                XYZR(7:9)=[X(IJKBLO4(IJKR)),Y(IJKBLO4(IJKR)),Z(IJKBLO4(IJKR))]
                                XYZR(10:12)=[X(IJKBLO3(IJKR)),Y(IJKBLO3(IJKR)),Z(IJKBLO3(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                !call reverseOrder(XYZCommon)
                                if (AR.gt.0.0d0) then
                                    STARTED=.true.
                                    FOUND=.true.
                                    NF=NF+1
                                    L(NF)=IJKPBLO(IJKL)
                                    R(NF)=IJKPBLO(IJKR)
                                    XF(NF)=XYZF(1)
                                    YF(NF)=XYZF(2)
                                    ZF(NF)=XYZF(3)
                                    XCF(NF)=XC(R(NF))
                                    YCF(NF)=YC(R(NF))
                                    ZCF(NF)=ZC(R(NF))
                                    call calcGrad(L(NF),R(NF),XF(NF),YF(NF),ZF(NF),FF(NF))
                                    call normalArea(&
                                        XYZCommon,L(NF),R(NF),ARF(NF),&
                                        DNF(NF),XPNF(NF),YPNF(NF),ZPNF(NF),&
                                        NXF(NF),NYF(NF),NZF(NF))
                                    if (IJKR.eq.IJKENR) then
                                        IJKSTL=IJKL+1
                                        exit SouthOuter
                                    end if
                                    !if (equalSize) then
                                    !    IJKSTR=IJKR+1
                                    !    cycle SouthOuter
                                    !end if
                                else if (.not.equalSize) then
                                    cycle SouthInner
                                else if (AR.le.0.0d0.and.STARTED) then
                                    IJKSTR=IJKR
                                    STARTED=.false.
                                    cycle SouthOuter
                                else if (FOUND) then
                                    IJKSTL=IJKL
                                    FOUND=.false.
                                    exit SouthOuter
                                else
                                    cycle SouthInner
                                end if
                            end do SouthInner
                        end do SouthOuter
                    case(2) 
                        !
                        !..........NORTH..........
                        !
                        print *, 'NORTH'
                        NorthOuter: do IJKL=IJKSTL,IJKENL
                            !
                            XYZL(1:3)=[X(IJKBLO1(IJKL)),Y(IJKBLO1(IJKL)),Z(IJKBLO1(IJKL))]
                            XYZL(4:6)=[X(IJKBLO2(IJKL)),Y(IJKBLO2(IJKL)),Z(IJKBLO2(IJKL))]
                            XYZL(7:9)=[X(IJKBLO4(IJKL)),Y(IJKBLO4(IJKL)),Z(IJKBLO4(IJKL))]
                            XYZL(10:12)=[X(IJKBLO3(IJKL)),Y(IJKBLO3(IJKL)),Z(IJKBLO3(IJKL))]
                            !
                            NorthInner: do IJKR=IJKSTR,IJKENR
                                iterationsCounter=iterationsCounter+1
                                !print *, iterationsCounter
                                !print *, IJKPBLO(IJKL),IJKPBLO(IJKR)
                                !
                                XYZR(1:3)=[X(IJKBLO2(IJKR)),Y(IJKBLO2(IJKR)),Z(IJKBLO2(IJKR))]
                                XYZR(4:6)=[X(IJKBLO1(IJKR)),Y(IJKBLO1(IJKR)),Z(IJKBLO1(IJKR))]
                                XYZR(7:9)=[X(IJKBLO3(IJKR)),Y(IJKBLO3(IJKR)),Z(IJKBLO3(IJKR))]
                                XYZR(10:12)=[X(IJKBLO4(IJKR)),Y(IJKBLO4(IJKR)),Z(IJKBLO4(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                call reverseOrder(XYZCommon)
                                if (AR.gt.0.0d0) then
                                    STARTED=.true.
                                    FOUND=.true.
                                    NF=NF+1
                                    L(NF)=IJKPBLO(IJKL)
                                    R(NF)=IJKPBLO(IJKR)
                                    XF(NF)=XYZF(1)
                                    YF(NF)=XYZF(2)
                                    ZF(NF)=XYZF(3)
                                    XCF(NF)=XC(R(NF))
                                    YCF(NF)=YC(R(NF))
                                    ZCF(NF)=ZC(R(NF))
                                    call calcGrad(L(NF),R(NF),XF(NF),YF(NF),ZF(NF),FF(NF))
                                    call normalArea(&
                                        XYZCommon,L(NF),R(NF),ARF(NF),&
                                        DNF(NF),XPNF(NF),YPNF(NF),ZPNF(NF),&
                                        NXF(NF),NYF(NF),NZF(NF))
                                    if (IJKR.eq.IJKENR) then
                                        IJKSTL=IJKL+1
                                        exit NorthOuter
                                    end if
                                    !if (equalSize) then
                                    !    IJKSTR=IJKR+1
                                    !    cycle NorthOuter
                                    !end if
                                else if (.not.equalSize) then
                                    cycle NorthInner
                                else if (AR.le.0.0d0.and.STARTED) then
                                    IJKSTR=IJKR
                                    STARTED=.false.
                                    IJKBLOSTR=IJKMARKR
                                    cycle NorthOuter
                                else if (FOUND) then
                                    IJKSTL=IJKL
                                    FOUND=.false.
                                    exit NorthOuter
                                else
                                    cycle NorthInner
                                end if
                            end do NorthInner
                        end do NorthOuter
                    case(3)
                        !
                        !..........WEST..........
                        !
                        print *, 'WEST'
                        WestOuter: do IJKL=IJKSTL,IJKENL
                            !
                            XYZL(1:3)=[X(IJKBLO1(IJKL)),Y(IJKBLO1(IJKL)),Z(IJKBLO1(IJKL))]
                            XYZL(4:6)=[X(IJKBLO2(IJKL)),Y(IJKBLO2(IJKL)),Z(IJKBLO2(IJKL))]
                            XYZL(7:9)=[X(IJKBLO4(IJKL)),Y(IJKBLO4(IJKL)),Z(IJKBLO4(IJKL))]
                            XYZL(10:12)=[X(IJKBLO3(IJKL)),Y(IJKBLO3(IJKL)),Z(IJKBLO3(IJKL))]
                            !
                            WestInner: do IJKR=IJKSTR,IJKENR
                                iterationsCounter=iterationsCounter+1
                                !print *, IJKPBLO(IJKL),IJKPBLO(IJKR)
                                !print *, iterationsCounter
                                !
                                XYZR(1:3)=[X(IJKBLO2(IJKR)),Y(IJKBLO2(IJKR)),Z(IJKBLO2(IJKR))]
                                XYZR(4:6)=[X(IJKBLO1(IJKR)),Y(IJKBLO1(IJKR)),Z(IJKBLO1(IJKR))]
                                XYZR(7:9)=[X(IJKBLO3(IJKR)),Y(IJKBLO3(IJKR)),Z(IJKBLO3(IJKR))]
                                XYZR(10:12)=[X(IJKBLO4(IJKR)),Y(IJKBLO4(IJKR)),Z(IJKBLO4(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                call reverseOrder(XYZCommon)
                                if (AR.gt.0.0d0) then
                                    STARTED=.true.
                                    FOUND=.true.
                                    NF=NF+1
                                    L(NF)=IJKPBLO(IJKL)
                                    R(NF)=IJKPBLO(IJKR)
                                    !print *, L(NF),R(NF)
                                    XF(NF)=XYZF(1)
                                    YF(NF)=XYZF(2)
                                    ZF(NF)=XYZF(3)
                                    XCF(NF)=XC(R(NF))
                                    YCF(NF)=YC(R(NF))
                                    ZCF(NF)=ZC(R(NF))
                                    call calcGrad(L(NF),R(NF),XF(NF),YF(NF),ZF(NF),FF(NF))
                                    call normalArea(&
                                        XYZCommon,L(NF),R(NF),ARF(NF),&
                                        DNF(NF),XPNF(NF),YPNF(NF),ZPNF(NF),&
                                        NXF(NF),NYF(NF),NZF(NF))
                                    if (IJKR.eq.IJKENR) then
                                        IJKSTL=IJKL+1
                                        exit WestOuter
                                    end if
                                    !if (equalSize) then
                                    !    IJKSTR=IJKR+1
                                    !    cycle WestOuter
                                    !end if
                                else if (.not.equalSize) then
                                    cycle WestInner
                                else if (AR.le.0.0d0.and.STARTED) then
                                    IJKSTR=IJKR
                                    STARTED=.false.
                                    cycle WestOuter
                                else if (FOUND) then
                                    IJKSTL=IJKL
                                    FOUND=.false.
                                    exit WestOuter
                                else
                                    cycle WestInner
                                end if
                            end do WestInner
                        end do WestOuter
                    case(4)
                        !
                        !..........EAST..........
                        !
                        print *, 'EAST'
                        EastOuter: do IJKL=IJKSTL,IJKENL
                            !
                            XYZL(1:3)=[X(IJKBLO2(IJKL)),Y(IJKBLO2(IJKL)),Z(IJKBLO2(IJKL))]
                            XYZL(4:6)=[X(IJKBLO1(IJKL)),Y(IJKBLO1(IJKL)),Z(IJKBLO1(IJKL))]
                            XYZL(7:9)=[X(IJKBLO3(IJKL)),Y(IJKBLO3(IJKL)),Z(IJKBLO3(IJKL))]
                            XYZL(10:12)=[X(IJKBLO4(IJKL)),Y(IJKBLO4(IJKL)),Z(IJKBLO4(IJKL))]
                            !
                            EastInner: do IJKR=IJKSTR,IJKENR
                                iterationsCounter=iterationsCounter+1
                                !print *, IJKPBLO(IJKL),IJKPBLO(IJKR)
                                !print *, iterationsCounter
                                !
                                XYZR(1:3)=[X(IJKBLO1(IJKR)),Y(IJKBLO1(IJKR)),Z(IJKBLO1(IJKR))]
                                XYZR(4:6)=[X(IJKBLO2(IJKR)),Y(IJKBLO2(IJKR)),Z(IJKBLO2(IJKR))]
                                XYZR(7:9)=[X(IJKBLO4(IJKR)),Y(IJKBLO4(IJKR)),Z(IJKBLO4(IJKR))]
                                XYZR(10:12)=[X(IJKBLO3(IJKR)),Y(IJKBLO3(IJKR)),Z(IJKBLO3(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                if (AR.gt.0.0d0) then
                                    STARTED=.true.
                                    FOUND=.true.
                                    NF=NF+1
                                    L(NF)=IJKPBLO(IJKL)
                                    R(NF)=IJKPBLO(IJKR)
                                    !print *, L(NF),R(NF)
                                    XF(NF)=XYZF(1)
                                    YF(NF)=XYZF(2)
                                    ZF(NF)=XYZF(3)
                                    XCF(NF)=XC(R(NF))
                                    YCF(NF)=YC(R(NF))
                                    ZCF(NF)=ZC(R(NF))
                                    call calcGrad(L(NF),R(NF),XF(NF),YF(NF),ZF(NF),FF(NF))
                                    call normalArea(&
                                        XYZCommon,L(NF),R(NF),ARF(NF),&
                                        DNF(NF),XPNF(NF),YPNF(NF),ZPNF(NF),&
                                        NXF(NF),NYF(NF),NZF(NF))
                                    if (IJKR.eq.IJKENR) then
                                        IJKSTL=IJKL+1
                                        exit EastOuter
                                    end if
                                    !if (equalSize) then
                                    !    IJKSTR=IJKR+1
                                    !    cycle EastOuter
                                    !end if
                                else if (.not.equalSize) then
                                    cycle EastInner
                                else if (AR.le.0.0d0.and.STARTED) then
                                    IJKSTR=IJKR
                                    STARTED=.false.
                                    cycle EastOuter
                                else if (FOUND) then
                                    FOUND=.false.
                                    IJKSTL=IJKL
                                    exit EastOuter
                                else
                                    cycle EastInner
                                end if
                            end do EastInner
                        end do EastOuter
                    case(5)
                        !
                        !..........BOTTOM..........
                        !
                        print *, 'BOTTOM'
                        BottomOuter: do IJKL=IJKSTL,IJKENL
                            !
                            XYZL(1:3)=[X(IJKBLO1(IJKL)),Y(IJKBLO1(IJKL)),Z(IJKBLO1(IJKL))]
                            XYZL(4:6)=[X(IJKBLO2(IJKL)),Y(IJKBLO2(IJKL)),Z(IJKBLO2(IJKL))]
                            XYZL(7:9)=[X(IJKBLO4(IJKL)),Y(IJKBLO4(IJKL)),Z(IJKBLO4(IJKL))]
                            XYZL(10:12)=[X(IJKBLO3(IJKL)),Y(IJKBLO3(IJKL)),Z(IJKBLO3(IJKL))]
                            !
                            BottomInner: do IJKR=IJKSTR,IJKENR
                                iterationsCounter=iterationsCounter+1
                                !print *, IJKPBLO(IJKL),IJKPBLO(IJKR)
                                !print *, iterationsCounter
                                !
                                XYZR(1:3)=[X(IJKBLO2(IJKR)),Y(IJKBLO2(IJKR)),Z(IJKBLO2(IJKR))]
                                XYZR(4:6)=[X(IJKBLO1(IJKR)),Y(IJKBLO1(IJKR)),Z(IJKBLO1(IJKR))]
                                XYZR(7:9)=[X(IJKBLO3(IJKR)),Y(IJKBLO3(IJKR)),Z(IJKBLO3(IJKR))]
                                XYZR(10:12)=[X(IJKBLO4(IJKR)),Y(IJKBLO4(IJKR)),Z(IJKBLO4(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                call reverseOrder(XYZCommon)
                                if (AR.gt.0.0d0) then
                                    STARTED=.true.
                                    FOUND=.true.
                                    NF=NF+1
                                    L(NF)=IJKPBLO(IJKL)
                                    R(NF)=IJKPBLO(IJKR)
                                    XF(NF)=XYZF(1)
                                    YF(NF)=XYZF(2)
                                    ZF(NF)=XYZF(3)
                                    XCF(NF)=XC(R(NF))
                                    YCF(NF)=YC(R(NF))
                                    ZCF(NF)=ZC(R(NF))
                                    call calcGrad(L(NF),R(NF),XF(NF),YF(NF),ZF(NF),FF(NF))
                                    call normalArea(&
                                        XYZCommon,L(NF),R(NF),ARF(NF),&
                                        DNF(NF),XPNF(NF),YPNF(NF),ZPNF(NF),&
                                        NXF(NF),NYF(NF),NZF(NF))
                                    if (IJKR.eq.IJKENR) then
                                        IJKSTL=IJKL+1
                                        exit BottomOuter
                                    end if
                                    !if (equalSize) then
                                    !    IJKSTR=IJKR+1
                                    !    cycle BottomOuter
                                    !end if
                                else if (.not.equalSize) then
                                    cycle BottomInner
                                else if (AR.le.0.0d0.and.STARTED) then
                                    IJKSTR=IJKR
                                    STARTED=.false.
                                    cycle BottomOuter
                                else if (FOUND) then
                                    FOUND=.false.
                                    IJKSTL=IJKL
                                    exit BottomOuter
                                else
                                    cycle BottomInner
                                end if
                            end do BottomInner
                        end do BottomOuter
                    case(6) 
                        !
                        !..........TOP..........
                        !
                        print *, 'TOP'
                        TopOuter: do IJKL=IJKSTL,IJKENL
                            !
                            XYZL(1:3)=[X(IJKBLO2(IJKL)),Y(IJKBLO2(IJKL)),Z(IJKBLO2(IJKL))]
                            XYZL(4:6)=[X(IJKBLO1(IJKL)),Y(IJKBLO1(IJKL)),Z(IJKBLO1(IJKL))]
                            XYZL(7:9)=[X(IJKBLO3(IJKL)),Y(IJKBLO3(IJKL)),Z(IJKBLO3(IJKL))]
                            XYZL(10:12)=[X(IJKBLO4(IJKL)),Y(IJKBLO4(IJKL)),Z(IJKBLO4(IJKL))]
                            !
                            TopInner: do IJKR=IJKSTR,IJKENR
                                iterationsCounter=iterationsCounter+1
                                !print *, IJKPBLO(IJKL),IJKPBLO(IJKR)
                                !print *, iterationsCounter
                                !
                                XYZR(1:3)=[X(IJKBLO1(IJKR)),Y(IJKBLO1(IJKR)),Z(IJKBLO1(IJKR))]
                                XYZR(4:6)=[X(IJKBLO2(IJKR)),Y(IJKBLO2(IJKR)),Z(IJKBLO2(IJKR))]
                                XYZR(7:9)=[X(IJKBLO4(IJKR)),Y(IJKBLO4(IJKR)),Z(IJKBLO4(IJKR))]
                                XYZR(10:12)=[X(IJKBLO3(IJKR)),Y(IJKBLO3(IJKR)),Z(IJKBLO3(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                if (AR.gt.0.0d0) then
                                    STARTED=.true.
                                    FOUND=.true.
                                    NF=NF+1
                                    L(NF)=IJKPBLO(IJKL)
                                    R(NF)=IJKPBLO(IJKR)
                                    XF(NF)=XYZF(1)
                                    YF(NF)=XYZF(2)
                                    ZF(NF)=XYZF(3)
                                    XCF(NF)=XC(R(NF))
                                    YCF(NF)=YC(R(NF))
                                    ZCF(NF)=ZC(R(NF))
                                    call calcGrad(L(NF),R(NF),XF(NF),YF(NF),ZF(NF),FF(NF))
                                    call normalArea(&
                                        XYZCommon,L(NF),R(NF),ARF(NF),&
                                        DNF(NF),XPNF(NF),YPNF(NF),ZPNF(NF),&
                                        NXF(NF),NYF(NF),NZF(NF))
                                    if (IJKR.eq.IJKENR) then
                                        IJKSTL=IJKL+1
                                        exit TopOuter
                                    end if
                                    !if (equalSize) then
                                    !    IJKSTR=IJKR+1
                                    !    cycle TopOuter
                                    !end if
                                else if (.not.equalSize) then
                                    cycle TopInner
                                else if (AR.le.0.0d0.and.STARTED) then
                                    IJKSTR=IJKR
                                    STARTED=.false.
                                    cycle TopOuter
                                else if (FOUND) then
                                    FOUND=.false.
                                    IJKSTL=IJKL
                                    exit TopOuter
                                else
                                    cycle TopInner
                                end if
                            end do TopInner
                        end do TopOuter
                end select
            end if
        end do neighbour
        NFACEBL(B)=NF-FACEBL(B)
        print *, NFACEBL(B)
        call writeBlockData(B)
    end do

    print *, 'No. iterations: ', iterationsCounter

end subroutine findNeighbours

!#####################################################
subroutine writeBlockData(IB)
!#####################################################

    use boundaryModule
    use charModule
    use geoModule
    use indexModule
    implicit none
    integer,intent(in) :: IB
    BLOCKUNIT=OFFSET+IB
    call setBlockInd(IB)
    
    ! read in rest of relevant data
    read(BLOCKUNIT)  (FX(I), I=1,NIJK)
    read(BLOCKUNIT)  (FY(I), I=1,NIJK)
    read(BLOCKUNIT)  (FZ(I), I=1,NIJK)
    read(BLOCKUNIT)  DX,DY,DZ,VOL
    !read(BLOCKUNIT)  (SRDDIR(I),I=1,NDIR)
    !read(BLOCKUNIT)  (SRDNEU(I),I=1,NNEU)
    read(BLOCKUNIT)  (SRDWAL(I),I=1,NWAL)
    !
    ! overwrite data back to input file
    !
    write(BLOCK_CH,*) IB
    BLOCKFILE='grid_'//trim(adjustl(BLOCK_CH))//'.out'
    
    open(UNIT=BLOCKUNIT,FILE=BLOCKFILE,FORM=FORM_CH,POSITION='REWIND')
    print *, 'WRITING TO FILE: ', BLOCKFILE
    N=(NK-2)*(NI-2)*(NJ-2)
    write(BLOCKUNIT) NI,NJ,NK,NIJK,NDIR,NNEU,NWAL,NBLO,NFACE,N,IJKST,IJK_GLOBL(B)
    !write(BLOCKUNIT) (NEIGH(B,I),I=1,6)
    !write(BLOCKUNIT) (LK(KST+K),K=1,NK)
    !write(BLOCKUNIT) (LI(IST+I),I=1,NI)
    write(BLOCKUNIT) (X(IJKST+IJK),IJK=1,NIJK)
    write(BLOCKUNIT) (Y(IJKST+IJK),IJK=1,NIJK)
    write(BLOCKUNIT) (Z(IJKST+IJK),IJK=1,NIJK)
    write(BLOCKUNIT) (XC(IJKST+IJK),IJK=1,NIJK)
    write(BLOCKUNIT) (YC(IJKST+IJK),IJK=1,NIJK)
    write(BLOCKUNIT) (ZC(IJKST+IJK),IJK=1,NIJK)

    write(BLOCKUNIT) (IJKBBLO(IJKBLOST+IJK),IJK=1,NBLO)
    write(BLOCKUNIT) (IJKPBLO(IJKBLOST+IJK),IJK=1,NBLO)
    write(BLOCKUNIT) (IJKBLO1(IJKBLOST+IJK),IJK=1,NBLO)
    write(BLOCKUNIT) (IJKBLO2(IJKBLOST+IJK),IJK=1,NBLO)
    write(BLOCKUNIT) (IJKBLO3(IJKBLOST+IJK),IJK=1,NBLO)
    write(BLOCKUNIT) (IJKBLO4(IJKBLOST+IJK),IJK=1,NBLO)

    write(BLOCKUNIT) (IJKBDIR(IJKDIRST+IJK),IJK=1,NDIR)
    write(BLOCKUNIT) (IJKPDIR(IJKDIRST+IJK),IJK=1,NDIR)
    write(BLOCKUNIT) (IJKDIR1(IJKDIRST+IJK),IJK=1,NDIR)
    write(BLOCKUNIT) (IJKDIR2(IJKDIRST+IJK),IJK=1,NDIR)
    write(BLOCKUNIT) (IJKDIR3(IJKDIRST+IJK),IJK=1,NDIR)
    write(BLOCKUNIT) (IJKDIR4(IJKDIRST+IJK),IJK=1,NDIR)

    write(BLOCKUNIT) (IJKBNEU(IJKNEUST+IJK),IJK=1,NNEU)
    write(BLOCKUNIT) (IJKPNEU(IJKNEUST+IJK),IJK=1,NNEU)
    write(BLOCKUNIT) (IJKNEU1(IJKNEUST+IJK),IJK=1,NNEU)
    write(BLOCKUNIT) (IJKNEU2(IJKNEUST+IJK),IJK=1,NNEU)
    write(BLOCKUNIT) (IJKNEU3(IJKNEUST+IJK),IJK=1,NNEU)
    write(BLOCKUNIT) (IJKNEU4(IJKNEUST+IJK),IJK=1,NNEU)

    write(BLOCKUNIT) (IJKBWAL(IJKWALST+IJK),IJK=1,NWAL)
    write(BLOCKUNIT) (IJKPWAL(IJKWALST+IJK),IJK=1,NWAL)
    write(BLOCKUNIT) (IJKWAL1(IJKWALST+IJK),IJK=1,NWAL)
    write(BLOCKUNIT) (IJKWAL2(IJKWALST+IJK),IJK=1,NWAL)
    write(BLOCKUNIT) (IJKWAL3(IJKWALST+IJK),IJK=1,NWAL)
    write(BLOCKUNIT) (IJKWAL4(IJKWALST+IJK),IJK=1,NWAL)
    !
    ! untouched variables
    write(BLOCKUNIT) (FX(I), I=1,NIJK)
    write(BLOCKUNIT) (FY(I), I=1,NIJK)
    write(BLOCKUNIT) (FZ(I), I=1,NIJK)
    write(BLOCKUNIT) DX,DY,DZ,VOL
    write(BLOCKUNIT)  (SRDWAL(I),I=1,NWAL)
    !
    write(BLOCKUNIT) (L(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (R(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (XF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (YF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (ZF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (XCF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (YCF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (ZCF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (FF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (ARF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (DNF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (XPNF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (YPNF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (ZPNF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (NXF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (NYF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT) (NZF(FACEST+I),I=1,NFACE)
    !
    ! Map
    !
    write(BLOCKUNIT) (MIJK(IJK),IJK=1,NXYZA)
    write(BLOCKUNIT) (RMIJK(IJK),IJK=0,NA-1)
    close(unit=BLOCKUNIT)
    
end subroutine writeBlockData

!#####################################################
subroutine createMapping
!#####################################################
    
    use indexModule
    implicit none

    IJK_GLO=-1
    do B=1,NB
        NBL(B)=0
        call setBlockInd(B)
        do K=2,NK-1
        do I=2,NI-1
        do J=2,NJ-1
            NBL(B)=NBL(B)+1
            IJK_GLO=IJK_GLO+1
            IJK_LOC=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            MIJK(IJK_LOC)=IJK_GLO
            !print *, IJK_GLO, IJK_LOC,MIJK(IJK_LOC)
            RMIJK(IJK_GLO)=IJK_LOC
        end do
        end do
        end do
        IJK_GLOBL(B)=IJK_GLO-NBL(B)+1
    end do

    NA = IJK_GLO+1

end subroutine createMapping

!#####################################################
subroutine writeParameterModule
!#####################################################

    use parameterModule
    use indexModule
    implicit none

    OPEN(UNIT=9,FILE='parameterModule.f90')
    REWIND 9
    read(9,*) ! 'module parameterModule'
    read(9,*) ! 'implicit none'
    write(9,'(4X, A22, A4, I12)') 'integer, parameter :: ', 'NXA=', NXMAX
    write(9,'(4X, A22, A4, I12)') 'integer, parameter :: ', 'NYA=', NYMAX
    write(9,'(4X, A22, A4, I12)') 'integer, parameter :: ', 'NZA=', NZMAX
    write(9,'(4X, A22, A6, I12)') 'integer, parameter :: ', 'NXYZA=', NXYZMAX
    write(9,'(4X, A22, A7, I12)') 'integer, parameter :: ', 'NDIRAL=', NDIRMAX
    write(9,'(4X, A22, A7, I12)') 'integer, parameter :: ', 'NNEUAL=', NNEUMAX
    write(9,'(4X, A22, A7, I12)') 'integer, parameter :: ', 'NWALAL=', NWALMAX
    write(9,'(4X, A22, A7, I12)') 'integer, parameter :: ', 'NBLOAL=', NBLOMAX
    write(9,'(4X, A22, A8, I12)') 'integer, parameter :: ', 'NBLOCKS=',NBLOCKSMAX
    write(9,'(4X, A22, A5, I1)') 'integer, parameter :: ', 'PREC=',PREC
    write(9,'(4X, A22, A4, I12)') 'integer, parameter :: ', 'NAL=', NA
    write(9,'(4X, A22, A8, I12)') 'integer, parameter :: ', 'NFACEAL=',NFMAX
    write(9,'(4X, A22, A8, I4)') 'integer, parameter :: ', 'NPROCS=',NPROCSA
    write(9,'(A)') 'end module parameterModule'
    print *, 'NFACEAL= ', NF, NPROCSA, NF/NPROCSA

end subroutine writeParameterModule

!#####################################################
subroutine calcProcessorLoad
!#####################################################

    use parameterModule
    use charModule
    use controlModule
    use indexModule
    implicit none

    integer :: NXMAX_TEMP,NYMAX_TEMP,NZMAX_TEMP,&
               NXYZMAX_TEMP,NDIRMAX_TEMP,NNEUMAX_TEMP,NWALMAX_TEMP,NBLOMAX_TEMP,&
               NBLOCKSMAX_TEMP,NFMAX_TEMP

    NXMAX=0
    NYMAX=0
    NXYZMAX=0
    NDIRMAX=0
    NWALMAX=0
    NNEUMAX=0
    NBLOMAX=0
    NBLOCKSMAX=0
    NFMAX=0

    do PROC=0,NPROCSA-1
        write(PROC_CH,*) PROC
        PROCUNIT=PROC+PROCOFFSET
        PROCFILE='proc_'//trim(adjustl(PROC_CH))//'.inp'

        open(UNIT=PROCUNIT,FILE=PROCFILE)
        print *, 'opening ... ', PROCFILE
        rewind PROCUNIT
        read(PROCUNIT,*) NB
        ! allocate statement could handle load imbalance
        read(PROCUNIT,*) (B_GLO(B),B=1,NB)
        NXMAX_TEMP=sum(NIBL(B_GLO(1:NB)))
        NYMAX_TEMP=sum(NJBL(B_GLO(1:NB)))
        NZMAX_TEMP=sum(NKBL(B_GLO(1:NB)))
        NXYZMAX_TEMP=sum(NIJKBL(B_GLO(1:NB)))
        NDIRMAX_TEMP=sum(NDIRBL(B_GLO(1:NB)))
        NNEUMAX_TEMP=sum(NNEUBL(B_GLO(1:NB)))
        NWALMAX_TEMP=sum(NWALBL(B_GLO(1:NB)))
        NBLOMAX_TEMP=sum(NBLOBL(B_GLO(1:NB)))
        NFMAX_TEMP=sum(NFACEBL(B_GLO(1:NB)))

        NXMAX=max(NXMAX,NXMAX_TEMP)
        NYMAX=max(NYMAX,NYMAX_TEMP)
        NZMAX=max(NZMAX,NZMAX_TEMP)
        NXYZMAX=max(NXYZMAX,NXYZMAX_TEMP)
        NDIRMAX=max(NDIRMAX,NDIRMAX_TEMP)
        NNEUMAX=max(NNEUMAX,NNEUMAX_TEMP)
        NWALMAX=max(NWALMAX,NWALMAX_TEMP)
        NBLOMAX=max(NBLOMAX,NBLOMAX_TEMP)
        NBLOCKSMAX=max(NBLOCKSMAX,NB)
        NFMAX=max(NFMAX,NFMAX_TEMP)
    end do

end subroutine calcProcessorLoad
