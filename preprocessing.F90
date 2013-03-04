!#########################################################
program main
!#########################################################

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
    
    IBL(1)=0
    JBL(1)=0
    KBL(1)=0
    IJKBL(1)=0
    IJKBLOCKBL(1)=0
    NIBL(1)=NI
    NJBL(1)=NJ
    NKBL(1)=NK
    NIJKBL(1)=NIJK
    NBLOCKBL(1)=NBLOCK
    print *, IJKBLOCKBL(1),NBLOCK
    
    do B=2,NB
        BLOCKUNIT=OFFSET+B
        write(UNIT_CH,'(I1)') (BLOCKUNIT-OFFSET)
        BLOCKFILE='grid_'//trim(UNIT_CH)//'.out'
        open(UNIT=BLOCKUNIT,FILE=BLOCKFILE)
        rewind BLOCKUNIT
        read(BLOCKUNIT,*)   NI,NJ,NK,NIJK,NBLOCK
        !print *, NIJK

        BB=B-1
        IBL(B)=IBL(BB)+NIBL(BB)
        JBL(B)=JBL(BB)+NJBL(BB)
        KBL(B)=KBL(BB)+NKBL(BB)
        IJKBL(B)=IJKBL(BB)+NIJKBL(BB)
        IJKBLOCKBL(B)=IJKBLOCKBL(BB)+NBLOCKBL(BB)
        NIBL(B)=NI
        NJBL(B)=NJ
        NKBL(B)=NK
        NIJKBL(B)=NIJK
        NBLOCKBL(B)=NBLOCK
        print *, IJKBLOCKBL(B),NBLOCK
    end do
    
    do B=1,NB
        call setBlockInd(B,B)
        BLOCKUNIT=OFFSET+B
        read(BLOCKUNIT,*)   (NEIGH(B,I),I=1,6)
        read(BLOCKUNIT,*)   (LK(KSTL+K),K=1,NKL)
        read(BLOCKUNIT,*)   (LI(ISTL+I),I=1,NIL)
        read(BLOCKUNIT,*)   (X(IJKSTL+IJK),IJK=1,NIJKL)
        read(BLOCKUNIT,*)   (Y(IJKSTL+IJK),IJK=1,NIJKL)
        read(BLOCKUNIT,*)   (Z(IJKSTL+IJK),IJK=1,NIJKL)
        read(BLOCKUNIT,*)   (XC(IJKSTL+IJK),IJK=1,NIJKL)
        read(BLOCKUNIT,*)   (YC(IJKSTL+IJK),IJK=1,NIJKL)
        read(BLOCKUNIT,*)   (ZC(IJKSTL+IJK),IJK=1,NIJKL)
        read(BLOCKUNIT,*)   (IJKBBL(IJKBLOCKSTL+IJK),IJK=1,NBLOCKL)
        read(BLOCKUNIT,*)   (IJKPBL(IJKBLOCKSTL+IJK),IJK=1,NBLOCKL)
        print *, (IJKPBL(IJKBLOCKSTL+IJK),IJK=1,NBLOCKL)
        read(BLOCKUNIT,*)   (IJKBL1(IJKBLOCKSTL+IJK),IJK=1,NBLOCKL)
        read(BLOCKUNIT,*)   (IJKBL2(IJKBLOCKSTL+IJK),IJK=1,NBLOCKL)
        read(BLOCKUNIT,*)   (IJKBL3(IJKBLOCKSTL+IJK),IJK=1,NBLOCKL)
        read(BLOCKUNIT,*)   (IJKBL4(IJKBLOCKSTL+IJK),IJK=1,NBLOCKL)
    end do

    ! Remap values for IJKBL1-4
    print *, 'REMAPPING VALUES'
    do B=1,NB
        call setBlockInd(B,B)
        do IJK=1,NBLOCKL
            IJKBBL(IJKBLOCKSTL+IJK)=IJKBBL(IJKBLOCKSTL+IJK)+IJKSTL
            IJKPBL(IJKBLOCKSTL+IJK)=IJKPBL(IJKBLOCKSTL+IJK)+IJKSTL
            IJKBL1(IJKBLOCKSTL+IJK)=IJKBL1(IJKBLOCKSTL+IJK)+IJKSTL
            IJKBL2(IJKBLOCKSTL+IJK)=IJKBL2(IJKBLOCKSTL+IJK)+IJKSTL
            IJKBL3(IJKBLOCKSTL+IJK)=IJKBL3(IJKBLOCKSTL+IJK)+IJKSTL
            print *, IJK,IJKBLOCKSTL,IJKBL4(IJKBLOCKSTL+IJK),IJKSTL
            IJKBL4(IJKBLOCKSTL+IJK)=IJKBL4(IJKBLOCKSTL+IJK)+IJKSTL
            print *, IJKBL4(IJKBLOCKSTL+IJK)
        end do
    end do

    print *, 'REMAPPED VALUES'
    do B=1,NB
        call setBlockInd(B,B)
        print *,IJKSTL, (IJKPBL(IJKBLOCKSTL+IJK),IJK=1,NBLOCKL)
    end do

end subroutine readData

!#########################################################
subroutine findNeighbours
!#########################################################

    use bc
    use geo
    use preProcInd
    implicit none

    NIB=1
    NJB=2
    NKB=1
    IJKPLE=0
    IJKPRE=0
    FE=0

    ! use other block loop: B=1,NB
    do KB=1,NKB
    do IB=1,NIB
    do JB=1,NJB
        IJKB=(KB-1)*NIB*NJB+(IB-1)*NIB+JB
        FBL(IJKB)=FE
        print *, 'BLOCK: ',IJKB
        neighbour: do INEIGH=1,6
            if (NEIGH(IJKB,INEIGH) .gt. 0) then
                call setBlockInd(IJKB,NEIGH(IJKB,INEIGH))
                IJKPLST=IJKSTL+IJKPLE
                IJKPRST=IJKSTR+IJKPRE
                FST=FST+FE
                select case (INEIGH)
                    case (1)
                        !
                        !..........SOUTH..........
                        !
                        K1=0
                        do KL=1,NKML
                            IJKL=IJKSTL+LK(KSTL+KL)+LI(ISTL+1)+1
                            K1=K1+1                            
                            ZL(K1)=Z(IJKL)
                        end do
                        I1=0
                        do IL=1,NIML
                            IJKL=IJKSTL+LK(KSTL+1)+LI(ISTL+IL)+1
                            I1=I1+1
                            XL(I1)=X(IJKL)
                        end do
                        K2=0
                        do KR=1,NKMR
                            IJKR=IJKSTR+LK(KSTR+KR)+LI(ISTR+1)+1
                            K2=K2+1                            
                            ZR(K2)=Z(IJKR)
                        end do
                        I2=0
                        do IR=1,NIMR
                            IJKR=IJKSTR+LK(KSTR+1)+LI(ISTR+IR)+1
                            I2=I2+1
                            XR(I2)=X(IJKR)
                        end do
                        call calcFace(XL,ZL,XR,ZR,I1,K1,I2,K2,XF,ZF,NIJF,NIMF,NJMF)
                        print *, NIJF
                        NFBL(NBLOCKS)=NFBL(NBLOCKS)+NIJF
                        IJKPLE=IJKPLST+(K1-1)*(I1-1)
                        NIJKPL=IJKPLE-IJKPLST
                        IJKPRE=IJKPRST+(K2-1)*(I2-1)
                        NIJKPR=IJKPRE-IJKPRST
                        print *,IJKPLST,IJKPLE
                        do IJK=IJKBLOCKSTL+1,IJKBLOCKSTL+NBLOCKL
                            print *, IJKPBL(IJK),X(IJKPBL(IJK))
                        end do
                            
                        call findConnectivity &
                            (IJKPBL(IJKPLST+1:IJKPLE),&
                            IJKBL1(IJKPLST+1:IJKPLE),IJKBL4(IJKPLST+1:IJKPLE),IJKBL3(IJKPLST+1:IJKPLE),NIJKPL,&
                            IJKPBL(IJKPRST+1:IJKPRE),&
                            IJKBL2(IJKPRST+1:IJKPRE),IJKBL3(IJKPRST+1:IJKPRE),IJKBL4(IJKPRST+1:IJKPRE),NIJKPR,&
                            X,Z,XF,ZF,NIMF,NJMF,NIJK,L(FST+1:FE),R(FST+1:FE),NIJF)
                        do IR=1,NIJF
                            !print *, L(IR)
                        end do
                    case(2) 
                        !
                        !..........NORTH..........
                        !
                        print *, 'NORTH'
                        K1=0
                        do KL=1,NKML
                            IJKL=IJKSTL+LK(KSTL+KL)+LI(ISTL+1)+NJML
                            ZL(K1)=Z(IJKL)
                        end do
                        I1=0
                        do IL=1,NIML
                            IJKL=IJKSTL+LK(KSTL+1)+LI(ISTL+IL)+NJML
                            XL(I1)=X(IJKL)
                        end do
                        K2=0
                        do KR=1,NKMR
                            IJKR=IJKSTR+LK(KSTR+KR)+LI(ISTR+1)+NJMR
                            ZR(K2)=Z(IJKR)
                        end do
                        I2=0
                        do IR=1,NIMR
                            IJKR=IJKSTR+LK(KSTR+1)+LI(ISTR+IR)+NJMR
                            XR(I2)=X(IJKR)
                        end do
                        call calcFace(XL,ZL,XR,ZR,NIML,NKML,NIMR,NKMR,XF,ZF,NIJF,NIMF,NJMF)
                        print *, NIJF
                        NFBL(NBLOCKS)=NFBL(NBLOCKS)+NIJF
                        IJKPLE=IJKPLST+(K1-1)*(I1-1)
                        NIJKPL=IJKPLE-IJKPLST
                        IJKPRE=IJKPRST+(K2-1)*(I2-1)
                        NIJKPR=IJKPRE-IJKPRST
                        print *,IJKPLST,IJKPLE
                        
                        do IJK=IJKBLOCKSTL+1,IJKBLOCKSTL+NBLOCKL
                            print *, IJKPBL(IJK),X(IJKPBL(IJK))
                        end do

                        call findConnectivity &
                            (IJKPBL(IJKPLST+1:IJKPLE),&
                            IJKBL2(IJKPLST+1:IJKPLE),IJKBL3(IJKPLST+1:IJKPLE),IJKBL4(IJKPLST+1:IJKPLE),NIJKPL,&
                            IJKPBL(IJKPRST+1:IJKPRE),&
                            IJKBL1(IJKPRST+1:IJKPRE),IJKBL4(IJKPRST+1:IJKPRE),IJKBL3(IJKPRST+1:IJKPRE),NIJKPR,&
                            X,Z,XF,ZF,NIMF,NJMF,NIJK,L(FST+1:FE),R(FST+1:FE),NIJF)
                        do IR=1,NIJF
                            print *, L(IR)
                        end do
                    case(3)
                        print *, 'WEST'
                        !
                        !..........WEST..........
                        !
                        K1=0
                        do KL=1,NKML
                            IJKL=IJKSTL+LK(KSTL+KL)+LI(ISTL+1)+1
                            K1=K1+1                            
                            ZL(K1)=Z(IJKL)
                        end do
                        J1=0
                        do JL=1,NJML
                            IJKL=IJKSTL+LK(KSTL+1)+LI(ISTL+1)+JL
                            J1=J1+1
                            YL(J1)=Y(IJKL)
                        end do
                        K2=0
                        do KR=1,NKMR
                            IJKR=IJKSTR+LK(KSTR+KR)+LI(ISTR+1)+1
                            K2=K2+1                            
                            ZR(K2)=Z(IJKR)
                        end do
                        J2=0
                        do JR=1,NJMR
                            IJKR=IJKSTR+LK(KSTR+1)+LI(ISTR+1)+JR
                            J2=J2+1
                            YR(J2)=Y(IJKR)
                        end do
                        call calcFace(YL,ZL,YR,ZR,J1,K1,J2,K2,YF,ZF,NIJF,NIMF,NJMF)
                        NFBL(NBLOCKS)=NFBL(NBLOCKS)+NIJF
                        IJKPLE=IJKPLST+(K1-1)*(J1-1)
                        NIJKPL=IJKPLE-IJKPLST
                        IJKPRE=IJKPRST+(K2-1)*(J2-1)
                        NIJKPR=IJKPRE-IJKPRST
                        call findConnectivity &
                            (IJKPBL(IJKPLST+1:IJKPLE),&
                            IJKBL2(IJKPLST+1:IJKPLE),IJKBL3(IJKPLST+1:IJKPLE),IJKBL4(IJKPLST+1:IJKPLE),NIJKPL,&
                            IJKPBL(IJKPRST+1:IJKPRE),&
                            IJKBL1(IJKPRST+1:IJKPRE),IJKBL4(IJKPRST+1:IJKPRE),IJKBL3(IJKPRST+1:IJKPRE),NIJKPR,&
                            Y,Z,YF,ZF,NIMF,NJMF,NIJK,L(FST+1:FE),R(FST+1:FE),NIJF)
                    case(4)
                        print *, 'EAST' 
                        !
                        !..........EAST..........
                        !
                        K1=0
                        do KL=1,NKML
                            IJKL=IJKSTL+LK(KSTL+KL)+LI(ISTL+1)+NJML
                            K1=K1+1                            
                            ZL(K1)=Z(IJKL)
                        end do
                        J1=0
                        do JL=1,NJML
                            IJKL=IJKSTL+LK(KSTL+1)+LI(ISTL+1)+JL
                            J1=J1+1
                            YL(J1)=Y(IJKL)
                        end do
                        K2=0
                        do KR=1,NKMR
                            IJKR=IJKSTR+LK(KSTR+KR)+LI(ISTR+1)+NJMR
                            K2=K2+1                            
                            ZR(K2)=Z(IJKR)
                        end do
                        J2=0
                        do JR=1,NJMR
                            IJKR=IJKSTR+LK(KSTR+1)+LI(ISTR+1)+JR
                            J2=J2+1
                            YR(J2)=Y(IJKR)
                        end do
                        call calcFace(YL,ZL,YR,ZR,J1,K1,J2,K2,YF,ZF,NIJF,NIMF,NJMF)
                        NFBL(NBLOCKS)=NFBL(NBLOCKS)+NIJF
                        IJKPLE=IJKPLST+(K1-1)*(J1-1)
                        NIJKPL=IJKPLE-IJKPLST
                        IJKPRE=IJKPRST+(K2-1)*(J2-1)
                        NIJKPR=IJKPRE-IJKPRST
                        call findConnectivity &
                            (IJKPBL(IJKPLST+1:IJKPLE),&
                            IJKBL1(IJKPLST+1:IJKPLE),IJKBL4(IJKPLST+1:IJKPLE),IJKBL3(IJKPLST+1:IJKPLE),NIJKPL,&
                            IJKPBL(IJKPRST+1:IJKPRE),&
                            IJKBL2(IJKPRST+1:IJKPRE),IJKBL3(IJKPRST+1:IJKPRE),IJKBL4(IJKPRST+1:IJKPRE),NIJKPR,&
                            Y,Z,YF,ZF,NIMF,NJMF,NIJK,L(FST+1:FE),R(FST+1:FE),NIJF)
                    case(5)
                        !
                        !..........BOTTOM..........
                        !
                        I1=0
                        do IL=1,NIML
                            IJKL=IJKSTL+LK(KSTL+1)+LI(ISTL+IL)+1
                            I1=I1+1
                            XL(I1)=X(IJKL)
                        end do
                        J1=0
                        do JL=1,NJML
                            IJKL=IJKSTL+LK(KSTL+1)+LI(ISTL+1)+JL
                            J1=J1+1
                            YL(J1)=Y(IJKL)
                        end do
                        I2=0
                        do IR=1,NIMR
                            IJKR=IJKSTR+LK(KSTR+1)+LI(ISTR+IR)+1
                            I2=I2+1
                            XR(I2)=X(IJKR)
                        end do
                        J2=0
                        do JR=1,NJMR
                            IJKR=IJKSTR+LK(KSTR+1)+LI(ISTR+1)+JR
                            J2=J2+1
                            YR(J2)=Y(IJKR)
                        end do
                        call calcFace(XL,YL,XR,YR,I1,J1,I2,J2,XF,YF,NIJF,NIMF,NJMF)
                        NFBL(NBLOCKS)=NFBL(NBLOCKS)+NIJF
                        IJKPLE=IJKPLST+(K1-1)*(J1-1)
                        NIJKPL=IJKPLE-IJKPLST
                        IJKPRE=IJKPRST+(K2-1)*(J2-1)
                        NIJKPR=IJKPRE-IJKPRST
                        call findConnectivity &
                            (IJKPBL(IJKPLST+1:IJKPLE),&
                            IJKBL2(IJKPLST+1:IJKPLE),IJKBL3(IJKPLST+1:IJKPLE),IJKBL4(IJKPLST+1:IJKPLE),NIJKPL,&
                            IJKPBL(IJKPRST+1:IJKPRE),&
                            IJKBL1(IJKPRST+1:IJKPRE),IJKBL4(IJKPRST+1:IJKPRE),IJKBL3(IJKPRST+1:IJKPRE),NIJKPR,&
                            Y,Z,YF,ZF,NIMF,NJMF,NIJK,L(FST+1:FE),R(FST+1:FE),NIJF)
                    case(6) 
                        !
                        !..........TOP..........
                        !
                        I1=0
                        do IL=1,NIML
                            IJKL=IJKSTL+LK(KSTL+NKML)+LI(ISTL+IL)+1
                            I1=I1+1
                            XL(I1)=X(IJKL)
                        end do
                        J1=0
                        do JL=1,NJML
                            IJKL=IJKSTL+LK(KSTL+NKML)+LI(ISTL+1)+JL
                            J1=J1+1
                            YL(J1)=Y(IJKL)
                        end do
                        I2=0
                        do IR=1,NIMR
                            IJKR=IJKSTR+LK(KSTR+NKMR)+LI(ISTR+IR)+1
                            I2=I2+1
                            XR(I2)=X(IJKR)
                        end do
                        J2=0
                        do JR=1,NJMR
                            IJKR=IJKSTR+LK(KSTR+NKMR)+LI(ISTR+1)+JR
                            J2=J2+1
                            YR(J2)=Y(IJKR)
                        end do
                        call calcFace(XL,YL,XR,YR,I1,J1,I2,J2,XF,YF,NIJF,NIMF,NJMF)
                        NFBL(NBLOCKS)=NFBL(NBLOCKS)+NIJF
                        IJKPLE=IJKPLST+(K1-1)*(J1-1)
                        NIJKPL=IJKPLE-IJKPLST
                        IJKPRE=IJKPRST+(K2-1)*(J2-1)
                        NIJKPR=IJKPRE-IJKPRST
                        call findConnectivity &
                            (IJKPBL(IJKPLST+1:IJKPLE),&
                            IJKBL1(IJKPLST+1:IJKPLE),IJKBL4(IJKPLST+1:IJKPLE),IJKBL3(IJKPLST+1:IJKPLE),NIJKPL,&
                            IJKPBL(IJKPRST+1:IJKPRE),&
                            IJKBL2(IJKPRST+1:IJKPRE),IJKBL3(IJKPRST+1:IJKPRE),IJKBL4(IJKPRST+1:IJKPRE),NIJKPR,&
                            X,Y,XF,YF,NIMF,NJMF,NIJK,L(FST+1:FE),R(FST+1:FE),NIJF)
                end select
            end if
        end do neighbour
    end do
    end do
    end do

end subroutine findNeighbours

!########################################################
subroutine calcFace(XL,YL,XR,YR,NIL,NJL,NIR,NJR,XF,YF,NIJF,NIMF,NJMF)
!########################################################

    implicit none
    real*8,intent(in) :: XL(NIL),YL(NJL),XR(NIR),YR(NIL)
    integer,intent(in) :: NIL,NJL,NIR,NJR
    real*8,intent(out) :: XF((NIL+NIR)*(NJL+NJR)),YF((NIL+NIR)*(NJL+NJR))
    integer,intent(out) :: NIJF,NIMF,NJMF

    real*8 :: XB(NIL+NIR),YB(NJL+NJR)
    integer :: I,IL,IR,J,JL,JR,IJ,NIF,NJF

    XB(1)=XL(1)
    I=1
    IL=2
    IR=2
    do
        if (XL(IL).lt.XR(IR)) then
            I=I+1
            XB(I)=XL(IL)
            IL=IL+1
        else if (XL(IL).gt.XR(IR)) then
            I=I+1
            XB(I)=XR(IR)
            IR=IR+1
        else
            I=I+1
            XB(I)=XL(IL)
            IL=IL+1
            IR=IR+1
        end if

        if (IR.eq.NIR+1.or.IL.eq.NIL+1) exit
    end do

    YB(1)=YL(1)
    J=1
    JL=2
    JR=2
    do
        if (YL(JL).lt.YR(JR)) then
            J=J+1
            YB(J)=YL(JL)
            JL=JL+1
        else if (YL(JL).gt.YR(JR)) then
            J=J+1
            YB(J)=YR(JR)
            JR=JR+1
        else
            J=J+1
            YB(J)=YL(JL)
            JL=JL+1
            JR=JR+1
        end if

        if ((JR.eq.NJR+1).or.(JL.eq.NJL+1)) exit
    end do

    IJ=0
    NIMF=I
    NJMF=J
    NIF=I+1
    NJF=J+1
    NIJF=(NIMF-1)*(NJMF-1)
    do I=1,NIMF
    do J=1,NJMF
        IJ=(I-1)*NJF+J
        XF(IJ)=XB(I)
        YF(IJ)=YB(J)
    end do
    end do

end subroutine calcFace

!########################################################
subroutine findConnectivity(IJKPL,IJK2L,IJK3L,IJK4L,NIJKL,IJKPR,IJK2R,IJK3R,IJK4R,NIJKR,X,Y,XF,YF,NIMF,NJMF,NIJK,L,R,NIJF)
!########################################################
! Beispiel: NORTH FACE -> 2 recht unten, 3 rechts oben, 4 links oben

    implicit none
    integer,intent(in)  :: NIJKL,NIJKR,NIMF,NJMF,NIJK,NIJF
    integer,intent(in)  :: IJKPL(NIJKL),IJK2L(NIJKL),IJK3L(NIJKL),IJK4L(NIJKL),IJKPR(NIJKR),IJK2R(NIJKR),IJK3R(NIJKR),IJK4R(NIJKR)
    integer,intent(in out) :: L(NIJF),R(NIJF)
    integer :: F,IJKL,IJKR,I,J,IJ,NJ,IJKPLL,IJKPRR,IJK2LL,IJK3LL,IJK4LL,IJK2RR,IJK3RR,IJK4RR
    real,intent(in) :: X(NIJK),Y(NIJK),XF((NIMF+1)*(NJMF+1)),YF((NIMF+1)*(NJMF+1))

    F=0
    NJ=NJMF+1
    do I=2,NIMF
    do J=2,NJMF
        IJ=(I-1)*NJ+J
        F=F+1
        !
        R(F)=IJKPRR
        !
        do IJKL=1,NIJKL
            IJKPLL=IJKPL(IJKL)
            IJK2LL=IJK2L(IJKL)
            IJK3LL=IJK3L(IJKL)
            IJK4LL=IJK4L(IJKL)
            print *, IJK3LL,X(IJK3LL),Y(IJK3LL)
            !
            if ((XF(IJ).ge.X(IJK2LL).and.YF(IJ).le.Y(IJK2LL)).and. &
                (XF(IJ).le.X(IJK3LL).and.YF(IJ).le.Y(IJK3LL)).and. &
                (XF(IJ).le.X(IJK4LL).and.YF(IJ).ge.Y(IJK4LL)))  then
            !
                L(F)=IJKPLL
                exit
            end if
        end do

        do IJKR=1,NIJKR
            IJKPRR=IJKPR(IJKR)
            IJK2RR=IJK2R(IJKR)
            IJK3RR=IJK3R(IJKR)
            IJK4RR=IJK4R(IJKR)
            !
            if ((XF(IJ).ge.X(IJK2RR).and.YF(IJ).le.Y(IJK2RR)).and. &
                (XF(IJ).le.X(IJK3RR).and.YF(IJ).le.Y(IJK3RR)).and. &
                (XF(IJ).le.X(IJK4RR).and.YF(IJ).ge.Y(IJK4RR)))  then
            !
                R(F)=IJKPRR
                exit
            end if
        end do
    end do
    end do

end subroutine findConnectivity
