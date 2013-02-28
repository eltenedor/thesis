!#########################################################
program main
!#########################################################

    call readData
    call findNeighbours

end program main

!##########################################################
subroutine readData
!#########################################################

    use geo
    use preProcInd
    use param
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
    read(BLOCKUNIT,*)   NI,NJ,NK,NIJK
    
    IBL(1)=0
    JBL(1)=0
    KBL(1)=0
    IJKBL(1)=0
    NIBL(1)=NI
    NJBL(1)=NJ
    NKBL(1)=NK
    NIJKBL(1)=NIJK
    
    do B=2,NB
        BLOCKUNIT=OFFSET+B
        write(UNIT_CH,'(I1)') (BLOCKUNIT-OFFSET)
        BLOCKFILE='grid_'//trim(UNIT_CH)//'.out'
        open(UNIT=BLOCKUNIT,FILE=BLOCKFILE)
        rewind BLOCKUNIT
        read(BLOCKUNIT,*)   NI,NJ,NK,NIJK

        BB=B-1
        IBL(B)=IBL(BB)+NIBL(BB)
        JBL(B)=JBL(BB)+NJBL(BB)
        KBL(B)=KBL(BB)+NKBL(BB)
        IJKBL(B)=IJKBL(BB)+NIJKBL(BB)
        NIBL(B)=NI
        NJBL(B)=NJ
        NKBL(B)=NK
        NIJKBL(B)=NIJK
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
    end do

    !do B=1,2
    !do I=1,6
        !print *, B,I,NEIGH(B,I)
    !end do
    !end do
!
    !print *, 1,NIJKBL(2)
    !do B=1,NB
    !call setBlockInd(B,B)
    !do IJK=1,NIJKL
        !print *, IJKSTL+IJK, X(IJKSTL+IJK),Y(IJKSTL+IJK),Z(IJKSTL+IJK)
    !end do
    !end do
    !stop

end subroutine readData

!#########################################################
subroutine findNeighbours
!#########################################################

    use geo
    use preProcInd
    implicit none

    NIB=1
    NJB=2
    NKB=1

    do KB=1,NKB
    do IB=1,NIB
    do JB=1,NJB
        IJKB=(KB-1)*NIB*NJB+(IB-1)*NIB+JB
        print *, 'BLOCK: ',IJKB
        neighbour: do INEIGH=1,6
            if (NEIGH(IJKB,INEIGH) .gt. 0) then
                call setBlockInd(IJKB,NEIGH(IJKB,INEIGH))
                select case (INEIGH)
                    case (1)
                        !
                        !..........SOUTH..........
                        !
                        IJKPLST=IJKSTL+IJKPLE
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
                        IJKPLE=IJKPLST+(K1-1)*(I1-1)
                        NIJKPL=IJKPLE-IJKPLST
                        IJKPRST=IJKSTR+IJKPRE
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
                        IJKPRE=IJKPRST+(K2-1)*(I2-1)
                        NIJKPR=IJKPRE-IJKPRST
                        call calcFace(XL,ZL,XR,ZR,I1,K1,I2,K2,XF,ZF,NIJF)
                        call findConnectivity(IJKPBL(IJKPLST+1:IJKPLE),IJKBL3(IJKPLST+1:IJKPLE),NIJKPL,IJKPBL(IJKPRST+1:IJKPRE),IJKBL4(IJKPRST+1:IJKPRE),NIJKPR,X,Z,XF,ZF,NIMF,NJMF,NIJK,L(FST+1:FE),R(FST+1:FE),NIJF)
                    case(2) 
                        !
                        !..........NORTH..........
                        !
                        K1=0
                        do KL=1,NKML
                            IJKL=IJKSTL+LK(KSTL+KL)+LI(ISTL+1)+NJML
                            K1=K1+1                            
                            ZL(K1)=Z(IJKL)
                        end do
                        I1=0
                        do IL=1,NIML
                            IJKL=IJKSTL+LK(KSTL+1)+LI(ISTL+IL)+NJML
                            I1=I1+1
                            XL(I1)=X(IJKL)
                        end do
                        K2=0
                        do KR=1,NKMR
                            IJKR=IJKSTR+LK(KSTR+KR)+LI(ISTR+1)+NJMR
                            K2=K2+1                            
                            ZR(K2)=Z(IJKR)
                        end do
                        I2=0
                        do IR=1,NIMR
                            IJKR=IJKSTR+LK(KSTR+1)+LI(ISTR+IR)+NJMR
                            I2=I2+1
                            XR(I2)=X(IJKR)
                        end do
                        call calcFace(XL,ZL,XR,ZR,I1,K1,I2,K2,XF,ZF,NIJF)
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
                        call calcFace(YL,ZL,YR,ZR,J1,K1,J2,K2,YF,ZF,NIJF)
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
                        call calcFace(YL,ZL,YR,ZR,J1,K1,J2,K2,YF,ZF,NIJF)
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
                        call calcFace(XL,YL,XR,YR,I1,J1,I2,J2,XF,YF,NIJF)
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
                        call calcFace(XL,YL,XR,YR,I1,J1,I2,J2,XF,YF,NIJF)
                end select
            end if
        end do neighbour
    end do
    end do
    end do

end subroutine findNeighbours

!########################################################
subroutine calcFace(XL,YL,XR,YR,NIL,NJL,NIR,NJR,XF,YF,NIJF)
!########################################################

    implicit none
    real*8,intent(in) :: XL(NIL),YL(NJL),XR(NIR),YR(NIL)
    integer,intent(in) :: NIL,NJL,NIR,NJR
    real*8,intent(out) :: XF((NIL+NIR)*(NJL+NJR)),YF((NIL+NIR)*(NJL+NJR))
    integer,intent(out) :: NIJF

    real*8 :: XB(NIL+NIR),YB(NJL+NJR)
    integer :: I,IL,IR,J,JL,JR,IJ,NIF,NJF,NIMF,NJMF

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
            print *, YL(JL)
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
    NIJF=(NIF)*(NJF)
    do I=1,NIMF
    do J=1,NJMF
        IJ=(I-1)*NJF+J
        XF(IJ)=XB(I)
        YF(IJ)=YB(J)
        print *, IJ, XF(IJ),YF(IJ)
    end do
    end do

end subroutine calcFace

!########################################################
subroutine findConnectivity(IJKPL,IJK3L,NIJKL,IJKPR,IJK4R,NIJKR,X,Y,XF,YF,NIMF,NJMF,NIJK,L,R,NIJF)
!########################################################

    implicit none
    integer,intent(in)  :: NIJKL,NIJKR,NIMF,NJMF,NIJK,NIJF
    integer,intent(in)  :: IJKPL(NIJKL),IJK3L(NIJKL),IJKPR(NIJKR),IJK4R(NIJKR)
    integer,intent(out) :: L(NIJF),R(NIJF)
    integer :: F,IJKL,IJKR,I,J,IJ,IJKPLL,IJKPRR,IJK3LL,IJK4RR
    real,intent(in) :: X(NIJK),Y(NIJK),XF((NIMF+1)*(NJMF+1)),YF((NIMF+1)*(NJMF+1))

    F=0
    IJKL=1
    IJKR=1
    do I=2,NIMF
    do J=2,NJMF
        IJ=(I-1)*NJ+J
        F=F+1
        !
        IJKPLL=IJKPL(IJKL)
        IJKPRR=IJKPR(IJKR)
        IJK3LL=IJK3L(IJKL)
        IJK4RR=IJK4R(IJKR)
        !
        L(F)=IJKPLL
        R(F)=IJKPRR
        !
        if (XF(IJ).eq.X(IJK3LL).and.YF(IJ).eq.Y(IJK3LL)) IJKL=IJKL+1
        if (XF(IJ).eq.X(IJK4RR).and.YF(IJ).eq.Y(IJK4RR)) IJKR=IJKR+1
    end do
    end do

end subroutine findConnectivity
