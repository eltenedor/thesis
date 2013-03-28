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
    
    NB=NBLOCKS
    OFFSET=20
    BLOCKUNIT=OFFSET+1
    write(BLOCK_CH,'(I1)') (BLOCKUNIT-OFFSET)
    BLOCKFILE='grid_'//trim(BLOCK_CH)//'.pre'
    
    open(UNIT=BLOCKUNIT,FILE=BLOCKFILE)
    rewind BLOCKUNIT
    read(BLOCKUNIT,*) NI,NJ,NK,NIJK,NBLOCK,NDIR
    
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
    
    do B=2,NB
        BLOCKUNIT=OFFSET+B
        write(BLOCK_CH,'(I1)') (BLOCKUNIT-OFFSET)
        BLOCKFILE='grid_'//trim(BLOCK_CH)//'.pre'
        open(UNIT=BLOCKUNIT,FILE=BLOCKFILE)
        rewind BLOCKUNIT
        read(BLOCKUNIT,*) NI,NJ,NK,NIJK,NBLOCK,NDIR

        BB=B-1
        IBL(B)=IBL(BB)+NIBL(BB)
        JBL(B)=JBL(BB)+NJBL(BB)
        KBL(B)=KBL(BB)+NKBL(BB)
        IJKBL(B)=IJKBL(BB)+NIJKBL(BB)
        IJKBLOCKBL(B)=IJKBLOCKBL(BB)+NBLOCKBL(BB)
        IJKDIRBL(B)=IJKDIRBL(BB)+NDIRBL(BB)
        NIBL(B)=NI
        NJBL(B)=NJ
        NKBL(B)=NK
        NIJKBL(B)=NIJK
        NBLOCKBL(B)=NBLOCK
        NDIRBL(B)=NDIR
    end do

    call createMapping

    do B=1,NB
        call setBlockInd(B)
        BLOCKUNIT=OFFSET+B
        read(BLOCKUNIT,*)   (NEIGH(B,I),I=1,6)
        !print *, (NEIGH(B,I),I=1,6)
        !read(BLOCKUNIT,*)   (LK(KST+K),K=1,NK)
        !read(BLOCKUNIT,*)   (LI(IST+I),I=1,NI)
        read(BLOCKUNIT,*)   (X(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT,*)   (Y(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT,*)   (Z(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT,*)   (XC(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT,*)   (YC(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT,*)   (ZC(IJKST+IJK),IJK=1,NIJK)
        read(BLOCKUNIT,*)   (IJKBBL(IJKBLOCKST+IJK),IJK=1,NBLOCK)
        read(BLOCKUNIT,*)   (IJKPBL(IJKBLOCKST+IJK),IJK=1,NBLOCK)
        read(BLOCKUNIT,*)   (IJKBL1(IJKBLOCKST+IJK),IJK=1,NBLOCK)
        read(BLOCKUNIT,*)   (IJKBL2(IJKBLOCKST+IJK),IJK=1,NBLOCK)
        read(BLOCKUNIT,*)   (IJKBL3(IJKBLOCKST+IJK),IJK=1,NBLOCK)
        read(BLOCKUNIT,*)   (IJKBL4(IJKBLOCKST+IJK),IJK=1,NBLOCK)
        read(BLOCKUNIT,*)   (IJKBDI(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*)   (IJKPDI(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*)   (IJKDI1(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*)   (IJKDI2(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*)   (IJKDI3(IJKDIRST+IJK),IJK=1,NDIR)
        read(BLOCKUNIT,*)   (IJKDI4(IJKDIRST+IJK),IJK=1,NDIR)
    end do

    ! Remap values for IJKBL1-4
    print *, 'REMAPPING VALUES'
    do B=1,NB
        call setBlockInd(B)
        do IJKBLOCK=IJKBLOCKST+1,IJKBLOCKST+NBLOCK
            IJKBBL(IJKBLOCK)=IJKBBL(IJKBLOCK)+IJKST
            IJKPBL(IJKBLOCK)=IJKPBL(IJKBLOCK)+IJKST
            IJKBL1(IJKBLOCK)=IJKBL1(IJKBLOCK)+IJKST
            IJKBL2(IJKBLOCK)=IJKBL2(IJKBLOCK)+IJKST
            IJKBL3(IJKBLOCK)=IJKBL3(IJKBLOCK)+IJKST
            IJKBL4(IJKBLOCK)=IJKBL4(IJKBLOCK)+IJKST
        end do
        !
        do IJKDIR=IJKDIRST+1,IJKDIRST+NDIR
            IJKBDI(IJKDIR)=IJKBDI(IJKDIR)+IJKST
            IJKPDI(IJKDIR)=IJKPDI(IJKDIR)+IJKST
            IJKDI1(IJKDIR)=IJKDI1(IJKDIR)+IJKST
            IJKDI2(IJKDIR)=IJKDI2(IJKDIR)+IJKST
            IJKDI3(IJKDIR)=IJKDI3(IJKDIR)+IJKST
            IJKDI4(IJKDIR)=IJKDI4(IJKDIR)+IJKST
        end do
        !do K=1,NK
        !    LK(KST+K)=LK(KST+K)+IJKST
        !end do
        !do I=1,NI
        !    LI(IST+I)=LI(IST+I)+IJKST
        !end do
    end do

    !print *, 'REMAPPED VALUES'
    !do B=1,NB
        !call setBlockInd(B)
        !print *, 'BLOCK: ', B
        !do IJKDIR=IJKDIRST+1,IJKDIRST+NDIR
        !    IJKBDI(IJKDIR)
        !print *, IJKPDI(IJKDIR)
        !    IJKDI1(IJKDIR)
        !    IJKDI2(IJKDIR)
        !    IJKDI3(IJKDIR)
        !    IJKDI4(IJKDIR)
    !    end do
    !end do

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
        IJKSTL=IJKBLOCKST+1
        IJKEL=IJKBLOCKST+NBLOCK
        print *, 'BLOCK: ',B
        neighbour: do INEIGH=1,6
            if (NEIGH(B,INEIGH).gt.0) then
                call setBlockInd(B,NEIGH(B,INEIGH))
                IJKSTR=IJKBLOCKSTR+1
                IJKER=IJKBLOCKSTR+NBLOCKR
                STARTED=.false.
                FOUND=.false.
                select case (INEIGH)
                    case (1)
                        !
                        !..........SOUTH..........
                        !
                        print *, 'SOUTH'
                        SouthOuter: do IJKL=IJKSTL,IJKEL
                            !
                            XYZL(1:3)=[X(IJKBL2(IJKL)),Y(IJKBL2(IJKL)),Z(IJKBL2(IJKL))]
                            XYZL(4:6)=[X(IJKBL1(IJKL)),Y(IJKBL1(IJKL)),Z(IJKBL1(IJKL))]
                            XYZL(7:9)=[X(IJKBL3(IJKL)),Y(IJKBL3(IJKL)),Z(IJKBL3(IJKL))]
                            XYZL(10:12)=[X(IJKBL4(IJKL)),Y(IJKBL4(IJKL)),Z(IJKBL4(IJKL))]
                            !
                            SouthInner: do IJKR=IJKSTR,IJKER
                                iterationsCounter=iterationsCounter+1
                                !print *, iterationsCounter
                                !print *, IJKPBL(IJKL),IJKPBL(IJKR)
                                !
                                XYZR(1:3)=[X(IJKBL1(IJKR)),Y(IJKBL1(IJKR)),Z(IJKBL1(IJKR))]
                                XYZR(4:6)=[X(IJKBL2(IJKR)),Y(IJKBL2(IJKR)),Z(IJKBL2(IJKR))]
                                XYZR(7:9)=[X(IJKBL4(IJKR)),Y(IJKBL4(IJKR)),Z(IJKBL4(IJKR))]
                                XYZR(10:12)=[X(IJKBL3(IJKR)),Y(IJKBL3(IJKR)),Z(IJKBL3(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                call reverseOrder(XYZCommon)
                                if (AR.gt.0.0d0) then
                                    STARTED=.true.
                                    FOUND=.true.
                                    NF=NF+1
                                    L(NF)=IJKPBL(IJKL)
                                    R(NF)=IJKPBL(IJKR)
                                    XF(NF)=XYZF(1)
                                    YF(NF)=XYZF(2)
                                    ZF(NF)=XYZF(3)
                                    call calcGrad(L(NF),R(NF),XF(NF),YF(NF),ZF(NF),FF(NF))
                                    call normalArea(&
                                        XYZCommon,L(NF),R(NF),ARF(NF),&
                                        DNF(NF),XPNF(NF),YPNF(NF),ZPNF(NF),&
                                        NXF(NF),NYF(NF),NZF(NF))
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
                        NorthOuter: do IJKL=IJKSTL,IJKEL
                            !
                            XYZL(1:3)=[X(IJKBL1(IJKL)),Y(IJKBL1(IJKL)),Z(IJKBL1(IJKL))]
                            XYZL(4:6)=[X(IJKBL2(IJKL)),Y(IJKBL2(IJKL)),Z(IJKBL2(IJKL))]
                            XYZL(7:9)=[X(IJKBL4(IJKL)),Y(IJKBL4(IJKL)),Z(IJKBL4(IJKL))]
                            XYZL(10:12)=[X(IJKBL3(IJKL)),Y(IJKBL3(IJKL)),Z(IJKBL3(IJKL))]
                            !
                            NorthInner: do IJKR=IJKSTR,IJKER
                                iterationsCounter=iterationsCounter+1
                                !print *, iterationsCounter
                                !print *, IJKPBL(IJKL),IJKPBL(IJKR)
                                !
                                XYZR(1:3)=[X(IJKBL2(IJKR)),Y(IJKBL2(IJKR)),Z(IJKBL2(IJKR))]
                                XYZR(4:6)=[X(IJKBL1(IJKR)),Y(IJKBL1(IJKR)),Z(IJKBL1(IJKR))]
                                XYZR(7:9)=[X(IJKBL3(IJKR)),Y(IJKBL3(IJKR)),Z(IJKBL3(IJKR))]
                                XYZR(10:12)=[X(IJKBL4(IJKR)),Y(IJKBL4(IJKR)),Z(IJKBL4(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                if (AR.gt.0.0d0) then
                                    STARTED=.true.
                                    FOUND=.true.
                                    NF=NF+1
                                    L(NF)=IJKPBL(IJKL)
                                    R(NF)=IJKPBL(IJKR)
                                    XF(NF)=XYZF(1)
                                    YF(NF)=XYZF(2)
                                    ZF(NF)=XYZF(3)
                                    call calcGrad(L(NF),R(NF),XF(NF),YF(NF),ZF(NF),FF(NF))
                                    call normalArea(&
                                        XYZCommon,L(NF),R(NF),ARF(NF),&
                                        DNF(NF),XPNF(NF),YPNF(NF),ZPNF(NF),&
                                        NXF(NF),NYF(NF),NZF(NF))
                                else if (AR.le.0.0d0.and.STARTED) then
                                    IJKSTR=IJKR
                                    STARTED=.false.
                                    IJKBLOCKSTR=IJKMARKR
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
                        WestOuter: do IJKL=IJKSTL,IJKEL
                            !
                            XYZL(1:3)=[X(IJKBL1(IJKL)),Y(IJKBL1(IJKL)),Z(IJKBL1(IJKL))]
                            XYZL(4:6)=[X(IJKBL2(IJKL)),Y(IJKBL2(IJKL)),Z(IJKBL2(IJKL))]
                            XYZL(7:9)=[X(IJKBL4(IJKL)),Y(IJKBL4(IJKL)),Z(IJKBL4(IJKL))]
                            XYZL(10:12)=[X(IJKBL3(IJKL)),Y(IJKBL3(IJKL)),Z(IJKBL3(IJKL))]
                            !
                            WestInner: do IJKR=IJKSTR,IJKER
                                iterationsCounter=iterationsCounter+1
                                !print *, IJKPBL(IJKL),IJKPBL(IJKR)
                                !print *, iterationsCounter
                                !
                                XYZR(1:3)=[X(IJKBL2(IJKR)),Y(IJKBL2(IJKR)),Z(IJKBL2(IJKR))]
                                XYZR(4:6)=[X(IJKBL1(IJKR)),Y(IJKBL1(IJKR)),Z(IJKBL1(IJKR))]
                                XYZR(7:9)=[X(IJKBL3(IJKR)),Y(IJKBL3(IJKR)),Z(IJKBL3(IJKR))]
                                XYZR(10:12)=[X(IJKBL4(IJKR)),Y(IJKBL4(IJKR)),Z(IJKBL4(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                call reverseOrder(XYZCommon)
                                if (AR.gt.0.0d0) then
                                    STARTED=.true.
                                    FOUND=.true.
                                    NF=NF+1
                                    L(NF)=IJKPBL(IJKL)
                                    R(NF)=IJKPBL(IJKR)
                                    XF(NF)=XYZF(1)
                                    YF(NF)=XYZF(2)
                                    ZF(NF)=XYZF(3)

                                    call calcGrad(L(NF),R(NF),XF(NF),YF(NF),ZF(NF),FF(NF))
                                    call normalArea(&
                                        XYZCommon,L(NF),R(NF),ARF(NF),&
                                        DNF(NF),XPNF(NF),YPNF(NF),ZPNF(NF),&
                                        NXF(NF),NYF(NF),NZF(NF))
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
                        EastOuter: do IJKL=IJKSTL,IJKEL
                            !
                            XYZL(1:3)=[X(IJKBL2(IJKL)),Y(IJKBL2(IJKL)),Z(IJKBL2(IJKL))]
                            XYZL(4:6)=[X(IJKBL1(IJKL)),Y(IJKBL1(IJKL)),Z(IJKBL1(IJKL))]
                            XYZL(7:9)=[X(IJKBL3(IJKL)),Y(IJKBL3(IJKL)),Z(IJKBL3(IJKL))]
                            XYZL(10:12)=[X(IJKBL4(IJKL)),Y(IJKBL4(IJKL)),Z(IJKBL4(IJKL))]
                            !
                            EastInner: do IJKR=IJKSTR,IJKER
                                iterationsCounter=iterationsCounter+1
                                !print *, IJKPBL(IJKL),IJKPBL(IJKR)
                                !print *, iterationsCounter
                                !
                                XYZR(1:3)=[X(IJKBL1(IJKR)),Y(IJKBL1(IJKR)),Z(IJKBL1(IJKR))]
                                XYZR(4:6)=[X(IJKBL2(IJKR)),Y(IJKBL2(IJKR)),Z(IJKBL2(IJKR))]
                                XYZR(7:9)=[X(IJKBL4(IJKR)),Y(IJKBL4(IJKR)),Z(IJKBL4(IJKR))]
                                XYZR(10:12)=[X(IJKBL3(IJKR)),Y(IJKBL3(IJKR)),Z(IJKBL3(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                if (AR.gt.0.0d0) then
                                    STARTED=.true.
                                    FOUND=.true.
                                    NF=NF+1
                                    L(NF)=IJKPBL(IJKL)
                                    R(NF)=IJKPBL(IJKR)
                                    XF(NF)=XYZF(1)
                                    YF(NF)=XYZF(2)
                                    ZF(NF)=XYZF(3)
                                    call calcGrad(L(NF),R(NF),XF(NF),YF(NF),ZF(NF),FF(NF))
                                    call normalArea(&
                                        XYZCommon,L(NF),R(NF),ARF(NF),&
                                        DNF(NF),XPNF(NF),YPNF(NF),ZPNF(NF),&
                                        NXF(NF),NYF(NF),NZF(NF))
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
                        BottomOuter: do IJKL=IJKSTL,IJKEL
                            !
                            XYZL(1:3)=[X(IJKBL1(IJKL)),Y(IJKBL1(IJKL)),Z(IJKBL1(IJKL))]
                            XYZL(4:6)=[X(IJKBL2(IJKL)),Y(IJKBL2(IJKL)),Z(IJKBL2(IJKL))]
                            XYZL(7:9)=[X(IJKBL4(IJKL)),Y(IJKBL4(IJKL)),Z(IJKBL4(IJKL))]
                            XYZL(10:12)=[X(IJKBL3(IJKL)),Y(IJKBL3(IJKL)),Z(IJKBL3(IJKL))]
                            !
                            BottomInner: do IJKR=IJKSTR,IJKER
                                iterationsCounter=iterationsCounter+1
                                !print *, IJKPBL(IJKL),IJKPBL(IJKR)
                                !print *, iterationsCounter
                                !
                                XYZR(1:3)=[X(IJKBL2(IJKR)),Y(IJKBL2(IJKR)),Z(IJKBL2(IJKR))]
                                XYZR(4:6)=[X(IJKBL1(IJKR)),Y(IJKBL1(IJKR)),Z(IJKBL1(IJKR))]
                                XYZR(7:9)=[X(IJKBL3(IJKR)),Y(IJKBL3(IJKR)),Z(IJKBL3(IJKR))]
                                XYZR(10:12)=[X(IJKBL4(IJKR)),Y(IJKBL4(IJKR)),Z(IJKBL4(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                call reverseOrder(XYZCommon)
                                if (AR.gt.0.0d0) then
                                    STARTED=.true.
                                    FOUND=.true.
                                    NF=NF+1
                                    L(NF)=IJKPBL(IJKL)
                                    R(NF)=IJKPBL(IJKR)
                                    XF(NF)=XYZF(1)
                                    YF(NF)=XYZF(2)
                                    ZF(NF)=XYZF(3)
                                    call calcGrad(L(NF),R(NF),XF(NF),YF(NF),ZF(NF),FF(NF))
                                    call normalArea(&
                                        XYZCommon,L(NF),R(NF),ARF(NF),&
                                        DNF(NF),XPNF(NF),YPNF(NF),ZPNF(NF),&
                                        NXF(NF),NYF(NF),NZF(NF))
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
                        TopOuter: do IJKL=IJKSTL,IJKEL
                            !
                            XYZL(1:3)=[X(IJKBL2(IJKL)),Y(IJKBL2(IJKL)),Z(IJKBL2(IJKL))]
                            XYZL(4:6)=[X(IJKBL1(IJKL)),Y(IJKBL1(IJKL)),Z(IJKBL1(IJKL))]
                            XYZL(7:9)=[X(IJKBL3(IJKL)),Y(IJKBL3(IJKL)),Z(IJKBL3(IJKL))]
                            XYZL(10:12)=[X(IJKBL4(IJKL)),Y(IJKBL4(IJKL)),Z(IJKBL4(IJKL))]
                            !
                            TopInner: do IJKR=IJKSTR,IJKER
                                iterationsCounter=iterationsCounter+1
                                !print *, IJKPBL(IJKL),IJKPBL(IJKR)
                                !print *, iterationsCounter
                                !
                                XYZR(1:3)=[X(IJKBL1(IJKR)),Y(IJKBL1(IJKR)),Z(IJKBL1(IJKR))]
                                XYZR(4:6)=[X(IJKBL2(IJKR)),Y(IJKBL2(IJKR)),Z(IJKBL2(IJKR))]
                                XYZR(7:9)=[X(IJKBL4(IJKR)),Y(IJKBL4(IJKR)),Z(IJKBL4(IJKR))]
                                XYZR(10:12)=[X(IJKBL3(IJKR)),Y(IJKBL3(IJKR)),Z(IJKBL3(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                if (AR.gt.0.0d0) then
                                    STARTED=.true.
                                    FOUND=.true.
                                    NF=NF+1
                                    L(NF)=IJKPBL(IJKL)
                                    R(NF)=IJKPBL(IJKR)
                                    XF(NF)=XYZF(1)
                                    YF(NF)=XYZF(2)
                                    ZF(NF)=XYZF(3)
                                    call calcGrad(L(NF),R(NF),XF(NF),YF(NF),ZF(NF),FF(NF))
                                    call normalArea(&
                                        XYZCommon,L(NF),R(NF),ARF(NF),&
                                        DNF(NF),XPNF(NF),YPNF(NF),ZPNF(NF),&
                                        NXF(NF),NYF(NF),NZF(NF))
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
        call writeBlockData(B)
    end do

    print *, 'No. iterations: ', iterationsCounter

end subroutine findNeighbours

!#####################################################
subroutine writeBlockData(IB)
!#####################################################

    use boundaryModule
    use geoModule
    use indexModule
    implicit none
    integer,intent(in) :: IB
    integer :: BLOCKUNIT,OFFSET
    character(len=20) :: UNIT_CH,BLOCKFILE
    OFFSET=20
    BLOCKUNIT=OFFSET+IB
    call setBlockInd(IB)
    
    !
    ! read in rest of relevant data
    read(BLOCKUNIT,*)  (FX(I), I=1,NIJK)
    read(BLOCKUNIT,*)  (FY(I), I=1,NIJK)
    read(BLOCKUNIT,*)  (FZ(I), I=1,NIJK)
    read(BLOCKUNIT,*)  DX,DY,DZ, VOL
    read(BLOCKUNIT,*)  (SRDDI(I),I=1,NDIR)
    !
    ! overwrite data back to input file
    !
    write(UNIT_CH,'(I1)') IB
    BLOCKFILE='grid_'//trim(UNIT_CH)//'.out'
    
    open(UNIT=BLOCKUNIT,FILE=BLOCKFILE)
    rewind BLOCKUNIT
    N=(NK-2)*(NI-2)*(NJ-2)
    write(BLOCKUNIT,*) NI,NJ,NK,NIJK,NBLOCK,NDIR,NFACE,N,IJKST
    !write(BLOCKUNIT,*) (NEIGH(B,I),I=1,6)
    !write(BLOCKUNIT,*) (LK(KST+K),K=1,NK)
    !write(BLOCKUNIT,*) (LI(IST+I),I=1,NI)
    write(BLOCKUNIT,*) (X(IJKST+IJK),IJK=1,NIJK)
    write(BLOCKUNIT,*) (Y(IJKST+IJK),IJK=1,NIJK)
    write(BLOCKUNIT,*) (Z(IJKST+IJK),IJK=1,NIJK)
    write(BLOCKUNIT,*) (XC(IJKST+IJK),IJK=1,NIJK)
    write(BLOCKUNIT,*) (YC(IJKST+IJK),IJK=1,NIJK)
    write(BLOCKUNIT,*) (ZC(IJKST+IJK),IJK=1,NIJK)
    write(BLOCKUNIT,*) (IJKBBL(IJKBLOCKST+IJK),IJK=1,NBLOCK)
    write(BLOCKUNIT,*) (IJKPBL(IJKBLOCKST+IJK),IJK=1,NBLOCK)
    write(BLOCKUNIT,*) (IJKBL1(IJKBLOCKST+IJK),IJK=1,NBLOCK)
    write(BLOCKUNIT,*) (IJKBL2(IJKBLOCKST+IJK),IJK=1,NBLOCK)
    write(BLOCKUNIT,*) (IJKBL3(IJKBLOCKST+IJK),IJK=1,NBLOCK)
    write(BLOCKUNIT,*) (IJKBL4(IJKBLOCKST+IJK),IJK=1,NBLOCK)
    write(BLOCKUNIT,*) (IJKBDI(IJKDIRST+IJK),IJK=1,NDIR)
    write(BLOCKUNIT,*) (IJKPDI(IJKDIRST+IJK),IJK=1,NDIR)
    write(BLOCKUNIT,*) (IJKDI1(IJKDIRST+IJK),IJK=1,NDIR)
    write(BLOCKUNIT,*) (IJKDI2(IJKDIRST+IJK),IJK=1,NDIR)
    write(BLOCKUNIT,*) (IJKDI3(IJKDIRST+IJK),IJK=1,NDIR)
    write(BLOCKUNIT,*) (IJKDI4(IJKDIRST+IJK),IJK=1,NDIR)
    !
    ! untouched variables
    write(BLOCKUNIT,*) (FX(I), I=1,NIJK)
    write(BLOCKUNIT,*) (FY(I), I=1,NIJK)
    write(BLOCKUNIT,*) (FZ(I), I=1,NIJK)
    write(BLOCKUNIT,*) DX,DY,DZ,VOL
    write(BLOCKUNIT,*) (SRDDI(I),I=1,NDIR)
    !
    write(BLOCKUNIT,*) (L(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT,*) (R(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT,*) (XF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT,*) (YF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT,*) (ZF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT,*) (FF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT,*) (ARF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT,*) (DNF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT,*) (XPNF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT,*) (YPNF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT,*) (ZPNF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT,*) (NXF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT,*) (NYF(FACEST+I),I=1,NFACE)
    write(BLOCKUNIT,*) (NZF(FACEST+I),I=1,NFACE)
    !
    ! Map
    !
    write(BLOCKUNIT,*) (MIJK(IJK),IJK=1,NXYZA)

    rewind BLOCKUNIT
    close (unit=BLOCKUNIT)
    
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
        end do
        end do
        end do
    end do

end subroutine createMapping

!#####################################################
subroutine writeParameterModule
!#####################################################

    use indexModule
    implicit none

    OPEN(UNIT=9,FILE='parameterModule.F90')
    REWIND 9
    read(9,*) !'module param'
    read(9,*) !'  implicit none'
    read(9,*) !'integer, parameter :: ', 'NXA=', NIA,'&' 
    read(9,*) !',NYA=', NJA, '&'
    read(9,*) !',NZA=', NKA, '&'
    read(9,*) !',NXYZA=', NIJKA, '&'
    read(9,*) !',NDIRAL=', NDIRA,'&'
    read(9,*) ! ',NBLOCKAL=', NBLOCKA, '&'
    read(9,*) !  ',NBLOCKS=',NB,'&'
    read(9,*) !',PREC=',PREC,'&'
    write(9,'(A9 I5 A1)') ',NFACEAL=',NF
    write(9,'(A)') 'end module parameterModule'

end subroutine writeParameterModule
