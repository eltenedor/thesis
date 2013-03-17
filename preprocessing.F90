!#########################################################
program main
!#########################################################

    use iso_c_binding
    implicit none

    interface cpp_interface

        subroutine commonFace(XYZL, XYZR, XYZCommon, AR) bind (C, name = "commonFace")

        use iso_c_binding
        implicit none

        real(kind = c_double),intent(inout),dimension(*) :: XYZL,XYZR,XYZCommon
        real(kind = c_double),intent(inout) :: AR

        end subroutine commonFace
    end interface cpp_interface

    call readData
    call findNeighbours

end program main

!##########################################################
subroutine readData
!#########################################################

    use bc
    use geo
    use param
    use preProcInd
    implicit none

    integer :: BLOCKUNIT,OFFSET
    character(len=20) :: UNIT_CH,BLOCKFILE
    
    NB=NBLOCKS
    OFFSET=20
    BLOCKUNIT=OFFSET+1
    write(UNIT_CH,'(I1)') (BLOCKUNIT-OFFSET)
    BLOCKFILE='grid_'//trim(UNIT_CH)//'.out'
    
    open(UNIT=BLOCKUNIT,FILE=BLOCKFILE)
    rewind BLOCKUNIT
    read(BLOCKUNIT,*)   NI,NJ,NK,NIJK,NBLOCK
    
    IJKBL(1)=0
    IJKBLOCKBL(1)=0
    NIJKBL(1)=NIJK
    NBLOCKBL(1)=NBLOCK
    !print *, IJKBLOCKBL(1),NBLOCK
    
    do B=2,NB
        BLOCKUNIT=OFFSET+B
        write(UNIT_CH,'(I1)') (BLOCKUNIT-OFFSET)
        BLOCKFILE='grid_'//trim(UNIT_CH)//'.out'
        open(UNIT=BLOCKUNIT,FILE=BLOCKFILE)
        rewind BLOCKUNIT
        read(BLOCKUNIT,*)   NI,NJ,NK,NIJK,NBLOCK
        !print *, NIJK

        BB=B-1
        IJKBL(B)=IJKBL(BB)+NIJKBL(BB)
        IJKBLOCKBL(B)=IJKBLOCKBL(BB)+NBLOCKBL(BB)
        NIJKBL(B)=NIJK
        NBLOCKBL(B)=NBLOCK
        !print *, IJKBLOCKBL(B),NBLOCK
    end do
    
    do B=1,NB
        call setBlockInd(B)
        BLOCKUNIT=OFFSET+B
        read(BLOCKUNIT,*)   (NEIGH(B,I),I=1,6)
        read(BLOCKUNIT,*)   !(LK(KST+K),K=1,NKL)
        read(BLOCKUNIT,*)   !(LI(IST+I),I=1,NIL)
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
    end do

    ! Remap values for IJKBL1-4
    print *, 'REMAPPING VALUES'
    do B=1,NB
        call setBlockInd(B)
        do IJK=1,NBLOCK
            IJKBBL(IJKBLOCKST+IJK)=IJKBBL(IJKBLOCKST+IJK)+IJKST
            IJKPBL(IJKBLOCKST+IJK)=IJKPBL(IJKBLOCKST+IJK)+IJKST
            IJKBL1(IJKBLOCKST+IJK)=IJKBL1(IJKBLOCKST+IJK)+IJKST
            IJKBL2(IJKBLOCKST+IJK)=IJKBL2(IJKBLOCKST+IJK)+IJKST
            IJKBL3(IJKBLOCKST+IJK)=IJKBL3(IJKBLOCKST+IJK)+IJKST
            IJKBL4(IJKBLOCKST+IJK)=IJKBL4(IJKBLOCKST+IJK)+IJKST
        end do
    end do

    !print *, 'REMAPPED VALUES'
    !do B=1,NB
    !    call setBlockInd(B)
    !    print *,IJKST, (IJKPBL(IJKBLOCKST+IJK),IJK=1,NBLOCK)
    !end do

end subroutine readData

!#########################################################
subroutine findNeighbours
!#########################################################

    use iso_c_binding
    use bc
    use geo
    use preProcInd
    implicit none

    interface cpp_interface

        subroutine commonFace(XYZL,XYZR,XYZCommon,AR,XYZF) bind (C, name = "commonFace")

        use iso_c_binding
        implicit none

        real(kind = c_double),intent(inout),dimension(*) :: XYZL,XYZR,XYZCommon,XYZF
        real(kind = c_double),intent(inout) :: AR

        end subroutine commonFace
    end interface cpp_interface


    NF=0
    ! use other block loop: B=1,NB
    do B=1,NB
        FACESTBL(B)=NF
        print *, 'BLOCK: ',B
        neighbour: do INEIGH=1,6
            if (NEIGH(B,INEIGH).gt.0) then
                call setBlockInd(B,NEIGH(B,INEIGH))
                select case (INEIGH)
                    case (1)
                        !
                        !..........SOUTH..........
                        !
                        !print *, IJKBLOCKSTL,NBLOCKL,IJKBLOCKSTR,NBLOCKR
                        do IJKL=IJKBLOCKSTL+1,IJKBLOCKSTL+NBLOCKL
                            !
                            XYZL(1:3)=[X(IJKBL2(IJKL)),Y(IJKBL2(IJKL)),Z(IJKBL2(IJKL))]
                            XYZL(4:6)=[X(IJKBL1(IJKL)),Y(IJKBL1(IJKL)),Z(IJKBL1(IJKL))]
                            XYZL(7:9)=[X(IJKBL3(IJKL)),Y(IJKBL3(IJKL)),Z(IJKBL3(IJKL))]
                            XYZL(10:12)=[X(IJKBL4(IJKL)),Y(IJKBL4(IJKL)),Z(IJKBL4(IJKL))]
                            !
                            do IJKR=IJKBLOCKSTR+1,IJKBLOCKSTR+NBLOCKR
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
                                    !print *,NF,L(NF),R(NF)
                                    !print *,XF(NF),YF(NF),ZF(NF)
                                    !print *, DNF(NF),XPNF(NF),YPNF(NF),ZPNF(NF)
                                    !print *, FF(NF)
                                end if
                            end do
                        end do
                    case(2) 
                        !
                        !..........NORTH..........
                        !
                        do IJKL=IJKBLOCKSTL+1,IJKBLOCKSTL+NBLOCKL
                            !
                            XYZL(1:3)=[X(IJKBL1(IJKL)),Y(IJKBL1(IJKL)),Z(IJKBL1(IJKL))]
                            XYZL(4:6)=[X(IJKBL2(IJKL)),Y(IJKBL2(IJKL)),Z(IJKBL2(IJKL))]
                            XYZL(7:9)=[X(IJKBL4(IJKL)),Y(IJKBL4(IJKL)),Z(IJKBL4(IJKL))]
                            XYZL(10:12)=[X(IJKBL3(IJKL)),Y(IJKBL3(IJKL)),Z(IJKBL3(IJKL))]
                            !
                            do IJKR=IJKBLOCKSTR+1,IJKBLOCKSTR+NBLOCKR
                                !
                                XYZR(1:3)=[X(IJKBL2(IJKR)),Y(IJKBL2(IJKR)),Z(IJKBL2(IJKR))]
                                XYZR(4:6)=[X(IJKBL1(IJKR)),Y(IJKBL1(IJKR)),Z(IJKBL1(IJKR))]
                                XYZR(7:9)=[X(IJKBL3(IJKR)),Y(IJKBL3(IJKR)),Z(IJKBL3(IJKR))]
                                XYZR(10:12)=[X(IJKBL4(IJKR)),Y(IJKBL4(IJKR)),Z(IJKBL4(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                if (AR.gt.0.0d0) then
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
                                end if
                            end do
                        end do
                    case(3)
                        !
                        !..........WEST..........
                        !
                        do IJKL=IJKBLOCKSTL+1,IJKBLOCKSTL+NBLOCKL
                            !
                            XYZL(1:3)=[X(IJKBL1(IJKL)),Y(IJKBL1(IJKL)),Z(IJKBL1(IJKL))]
                            XYZL(4:6)=[X(IJKBL2(IJKL)),Y(IJKBL2(IJKL)),Z(IJKBL2(IJKL))]
                            XYZL(7:9)=[X(IJKBL4(IJKL)),Y(IJKBL4(IJKL)),Z(IJKBL4(IJKL))]
                            XYZL(10:12)=[X(IJKBL3(IJKL)),Y(IJKBL3(IJKL)),Z(IJKBL3(IJKL))]
                            !
                            do IJKR=IJKBLOCKSTR+1,IJKBLOCKSTR+NBLOCKR
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
                                end if
                            end do
                        end do
                    case(4)
                        !
                        !..........EAST..........
                        !
                        do IJKL=IJKBLOCKSTL+1,IJKBLOCKSTL+NBLOCKL
                            !
                            XYZL(1:3)=[X(IJKBL2(IJKL)),Y(IJKBL2(IJKL)),Z(IJKBL2(IJKL))]
                            XYZL(4:6)=[X(IJKBL1(IJKL)),Y(IJKBL1(IJKL)),Z(IJKBL1(IJKL))]
                            XYZL(7:9)=[X(IJKBL3(IJKL)),Y(IJKBL3(IJKL)),Z(IJKBL3(IJKL))]
                            XYZL(10:12)=[X(IJKBL4(IJKL)),Y(IJKBL4(IJKL)),Z(IJKBL4(IJKL))]
                            !
                            do IJKR=IJKBLOCKSTR+1,IJKBLOCKSTR+NBLOCKR
                                !
                                XYZR(1:3)=[X(IJKBL1(IJKR)),Y(IJKBL1(IJKR)),Z(IJKBL1(IJKR))]
                                XYZR(4:6)=[X(IJKBL2(IJKR)),Y(IJKBL2(IJKR)),Z(IJKBL2(IJKR))]
                                XYZR(7:9)=[X(IJKBL4(IJKR)),Y(IJKBL4(IJKR)),Z(IJKBL4(IJKR))]
                                XYZR(10:12)=[X(IJKBL3(IJKR)),Y(IJKBL3(IJKR)),Z(IJKBL3(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                if (AR.gt.0.0d0) then
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
                                end if
                            end do
                        end do
                    case(5)
                        !
                        !..........BOTTOM..........
                        !
                        do IJKL=IJKBLOCKSTL+1,IJKBLOCKSTL+NBLOCKL
                            !
                            XYZL(1:3)=[X(IJKBL1(IJKL)),Y(IJKBL1(IJKL)),Z(IJKBL1(IJKL))]
                            XYZL(4:6)=[X(IJKBL2(IJKL)),Y(IJKBL2(IJKL)),Z(IJKBL2(IJKL))]
                            XYZL(7:9)=[X(IJKBL4(IJKL)),Y(IJKBL4(IJKL)),Z(IJKBL4(IJKL))]
                            XYZL(10:12)=[X(IJKBL3(IJKL)),Y(IJKBL3(IJKL)),Z(IJKBL3(IJKL))]
                            !
                            do IJKR=IJKBLOCKSTR+1,IJKBLOCKSTR+NBLOCKR
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
                                end if
                            end do
                        end do
                    case(6) 
                        !
                        !..........TOP..........
                        !
                        do IJKL=IJKBLOCKSTL+1,IJKBLOCKSTL+NBLOCKL
                            !
                            XYZL(1:3)=[X(IJKBL2(IJKL)),Y(IJKBL2(IJKL)),Z(IJKBL2(IJKL))]
                            XYZL(4:6)=[X(IJKBL1(IJKL)),Y(IJKBL1(IJKL)),Z(IJKBL1(IJKL))]
                            XYZL(7:9)=[X(IJKBL3(IJKL)),Y(IJKBL3(IJKL)),Z(IJKBL3(IJKL))]
                            XYZL(10:12)=[X(IJKBL4(IJKL)),Y(IJKBL4(IJKL)),Z(IJKBL4(IJKL))]
                            !
                            do IJKR=IJKBLOCKSTR+1,IJKBLOCKSTR+NBLOCKR
                                !
                                XYZR(1:3)=[X(IJKBL1(IJKR)),Y(IJKBL1(IJKR)),Z(IJKBL1(IJKR))]
                                XYZR(4:6)=[X(IJKBL2(IJKR)),Y(IJKBL2(IJKR)),Z(IJKBL2(IJKR))]
                                XYZR(7:9)=[X(IJKBL4(IJKR)),Y(IJKBL4(IJKR)),Z(IJKBL4(IJKR))]
                                XYZR(10:12)=[X(IJKBL3(IJKR)),Y(IJKBL3(IJKR)),Z(IJKBL3(IJKR))]
                                !
                                AR=0.0d0
                                call commonFace(XYZL,XYZR,XYZCommon,AR,XYZF)
                                if (AR.gt.0.0d0) then
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
                                end if
                            end do
                        end do
                end select
            end if
        end do neighbour
        NFACEBL(B)=NF-FACESTBL(B)
        call writeBlockData(B)
    end do

end subroutine findNeighbours

!#####################################################
subroutine writeBlockData(IB)
!#####################################################

    use bc
    use geo
    use preProcInd
    implicit none
    integer,intent(in) :: IB
    integer :: BLOCKUNIT,OFFSET
    OFFSET=20
    BLOCKUNIT=OFFSET+IB
    call setBlockInd(IB)
    
    read(BLOCKUNIT,*)  !(IJKBDI(I),I=1,NDIR)
    read(BLOCKUNIT,*)  !(IJKPDI(I),I=1,NDIR)
    read(BLOCKUNIT,*)  !(IJKDI1(I),I=1,NDIR)
    read(BLOCKUNIT,*)  !(IJKDI2(I),I=1,NDIR)
    read(BLOCKUNIT,*)  !(IJKDI3(I),I=1,NDIR)
    read(BLOCKUNIT,*)  !(IJKDI4(I),I=1,NDIR)

    read(BLOCKUNIT,*)  !(FX(I), I=1,NIJK)
    read(BLOCKUNIT,*)  !(FY(I), I=1,NIJK)
    read(BLOCKUNIT,*)  !(FZ(I), I=1,NIJK)

    read(BLOCKUNIT,*)  !DX,DY,DZ, VOL
    read(BLOCKUNIT,*)  !(SRDDI(I),I=1,NDIR)
    !
    ! write data back to file
    !
    write(BLOCKUNIT,*) NFACE
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
    ! update block data
    !
    rewind BLOCKUNIT
    close (unit=BLOCKUNIT)
    
end subroutine writeBlockData
