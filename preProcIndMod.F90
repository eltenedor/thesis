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
                NIMF,NJMF,NIJF,&
                IJKPLST,IJKPLE,IJKPRST,IJKPRE,NIJKPL,NIJKPR,&
                FBL(NBLOCKS),NFBL(NBLOCKS),FST,FE,&
                IJKBLOCKSTL,IJKBLOCKSTR,NBLOCKL,NBLOCKR,NBLOCK,&
                LK(100), LI(100),&
                IBL(NBLOCKS),JBL(NBLOCKS),KBL(NBLOCKS),IJKBL(NBLOCKS),&
                NIBL(NBLOCKS),NJBL(NBLOCKS),NKBL(NBLOCKS),NIJKBL(NBLOCKS),&
                IJKBLOCKBL(NBLOCKS),NBLOCKBL(NBLOCKS),&
                L(1000),R(1000)
                
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
    IJKBLOCKSTL=IJKBLOCKBL(BL)
    NIML=NIL-1
    NJML=NJL-1
    NKML=NKL-1
    NBLOCKL=NBLOCKBL(BL)

    NIR=NIBL(BR)
    NJR=NJBL(BR)
    NKR=NKBL(BR)
    NIJKR=NIJKBL(BR)
    ISTR=IBL(BR)
    JSTR=JBL(BR)
    KSTR=KBL(BR)
    IJKSTR=IJKBL(BR)
    IJKBLOCKSTR=IJKBLOCKBL(BR)
    NIMR=NIR-1
    NJMR=NJR-1
    NKMR=NKR-1
    NBLOCKR=NBLOCKBL(BR)

end subroutine setBlockInd
                
end module preProcInd
