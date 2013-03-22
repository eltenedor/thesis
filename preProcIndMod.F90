!########################################################
module preProcInd
!########################################################

    use param
    implicit none
    integer ::  &
                ! Regular Indices block independent
                I,J,K,IJK,NI,NJ,NK,NIJ,NIJK,&
                ! Indices for complete Blocks
                B,BB,NB,&
                ! Indices inside Blocks
                IJKL,IJKR,&
                ! Neighbour Block Indices
                INEIGH,NEIGH(NBLOCKS,6),&
                ! block dependent regular indices
                IST,IBL(NBLOCKS),NIBL(NBLOCKS),&
                JST,JBL(NBLOCKS),NJBL(NBLOCKS),&
                KST,KBL(NBLOCKS),NKBL(NBLOCKS),&
                IJKST,IJKBL(NBLOCKS),NIJKBL(NBLOCKS),&
                LI(1000),LK(1000),&
                ! block dependent indices (dirichlet boundary)
                IJKDIRST,IJKDIRBL(NBLOCKS),NDIRBL(NBLOCKS),NDIR,&
                ! block dependent indices (block boundary)
                IJKBLOCKBL(NBLOCKS),NBLOCKBL(NBLOCKS),&
                IJKBLOCKST,NBLOCK,IJKBLOCKSTL,NBLOCKL,IJKBLOCKSTR,NBLOCKR,&
                ! block dependent indices (boundary faces)
                NF,&
                FACEST,FACESTBL(NBLOCKS),NFACEBL(NBLOCKS),NFACE
                
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
    
    IST=IBL(B)
    JST=JBL(B)
    KST=KBL(B)
    NI=NIBL(B)
    NJ=NJBL(B)
    NK=NKBL(B)

    IJKST=IJKBL(B)
    NIJK=NIJKBL(B)

    IJKBLOCKST=IJKBLOCKBL(B)
    NBLOCK=NBLOCKBL(B)
    
    IJKDIRST=IJKDIRBL(B)
    NDIR=NDIRBL(B)
    
    FACEST=FACESTBL(B)
    NFACE=NFACEBL(B)

end subroutine setBlockInd1Int
                
end module preProcInd
