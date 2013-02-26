!########################################################
module preProcInd
!########################################################

    use param
    implicit none
    integer ::  I,J,K,IJ,IJK,NI,NJ,NK,NIJ,NIJK,NB,&
                B,BB,IB,JB,KB,IJB,IJKB,NIB,NJB,NKB,&
                INEIGH,NEIGH(NBLOCKS,6),&
                I1,J1,K1,I2,J2,K2,&
                IL,JL,KL,IJKL,IR,JR,KR,IJKR,&
                NIL,NJL,NKL,NIJKL,ISTL,JSTL,KSTL,IJKSTL,NIML,NJML,NKML,&
                NIR,NJR,NKR,NIJKR,ISTR,JSTR,KSTR,IJKSTR,NIMR,NJMR,NKMR,&
                LK(100), LI(100),&
                IBL(NBLOCKS),JBL(NBLOCKS),KBL(NBLOCKS),IJKBL(NBLOCKS),&
                NIBL(NBLOCKS),NJBL(NBLOCKS),NKBL(NBLOCKS),NIJKBL(NBLOCKS),&
                NIJF
                
contains

!########################################################
subroutine setBlockInd(BL,BR)
!########################################################

    implicit none
    integer, intent(in) :: BL,BR

    NIL=NIBL(BL)
    NJL=NJBL(BL)
    NKL=NKBL(BL)
    NIJKL=NIJKBL(BL)
    ISTL=IBL(BL)
    JSTL=JBL(BL)
    KSTL=KBL(BL)
    IJKSTL=IJKBL(BL)
    NIML=NIL-1
    NJML=NJL-1
    NKML=NKL-1

    NIR=NIBL(BR)
    NJR=NJBL(BR)
    NKR=NKBL(BR)
    NIJKR=NIJKBL(BR)
    ISTR=IBL(BR)
    JSTR=JBL(BR)
    KSTR=KBL(BR)
    IJKSTR=IJKBL(BR)
    NIMR=NIR-1
    NJMR=NJR-1
    NKMR=NKR-1

end subroutine setBlockInd
                
end module preProcInd
