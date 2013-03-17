!########################################################
module preProcInd
!########################################################

    use param
    implicit none
    integer ::  I,J,K,IJK,NI,NJ,NK,NIJ,NIJK,B,BB,NB,&
                IJKL,IJKR,&
                INEIGH,NEIGH(NBLOCKS,6),&
                IJKBL(NBLOCKS),NIJKBL(NBLOCKS),&
                IJKST,&
                IJKBLOCKBL(NBLOCKS),NBLOCKBL(NBLOCKS),&
                IJKBLOCKST,NBLOCK,IJKBLOCKSTL,IJKBLOCKSTR,NBLOCKL,NBLOCKR,&
                FACEST,NFACE,NF,FACESTBL(NBLOCKS),NFACEBL(NBLOCKS)
                
    public :: setBlockInd
    private :: setBlockInd2Int,setBlockInd1Int
    
    interface setBlockInd
        module procedure setBlockInd2Int,setBlockInd1Int
    end interface
                
contains

!########################################################
subroutine setBlockInd2Int(BL,BR)
!########################################################

    implicit none
    integer, intent(in) :: BL,BR

    IJKBLOCKSTL=IJKBLOCKBL(BL)
    NBLOCKL=NBLOCKBL(BL)
    
    IJKBLOCKSTR=IJKBLOCKBL(BR)
    NBLOCKR=NBLOCKBL(BR)

end subroutine setBlockInd2Int

!########################################################
subroutine setBlockInd1Int(B)
!########################################################

    implicit none
    integer, intent(in) :: B

    IJKST=IJKBL(B)
    NIJK=NIJKBL(B)

    IJKBLOCKST=IJKBLOCKBL(B)
    NBLOCK=NBLOCKBL(B)
    
    FACEST=FACESTBL(B)
    NFACE=NFACEBL(B)

end subroutine setBlockInd1Int
                
end module preProcInd
