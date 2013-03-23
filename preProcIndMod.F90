!########################################################
module preProcInd
!########################################################

    use param
    implicit none
    integer ::  &
                ! Time stepping indices
                ITIM,ITIMS,ITIME,&
                ! Regular Indices block independent
                I,J,K,IJK,NI,NIM,NJ,NJM,NK,NKM,NIJ,NIJK,N,IDI,IJKB,IJKP,&
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
                IJKST,IJKBL(NBLOCKS),NIJKBL(NBLOCKS),NBL(NBLOCKS),&
                !LI(1000),LK(1000),&
                ! block dependent indices (dirichlet boundary)
                IJKDIRST,IJKDIRBL(NBLOCKS),NDIRBL(NBLOCKS),NDIR,&
                ! block dependent indices (block boundary)
                IJKBLOCKBL(NBLOCKS),NBLOCKBL(NBLOCKS),&
                IJKBLOCKST,NBLOCK,IJKBLOCKSTL,NBLOCKL,IJKBLOCKSTR,NBLOCKR,&
                ! block dependent indices (boundary faces)
                NF,&
                FACEST,FACEBL(NBLOCKS),NFACEBL(NBLOCKS),NFACE,&
                ! mapping related indices
                IJK_GLO,IJK_LOC,MIJK(NXYZA),&
                B_GLO(NBLOCKS),& !NBL(NBLOCKS),
                IJKBL_GLO(NBLOCKS)
                
                
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
    NIM=NI-1
    NJ=NJBL(B)
    NJM=NJ-1
    NK=NKBL(B)
    NKM=NK-1

    IJKST=IJKBL(B)
    NIJK=NIJKBL(B)

    IJKBLOCKST=IJKBLOCKBL(B)
    NBLOCK=NBLOCKBL(B)
    
    IJKDIRST=IJKDIRBL(B)
    NDIR=NDIRBL(B)
    
    FACEST=FACEBL(B)
    NFACE=NFACEBL(B)

end subroutine setBlockInd1Int

end module preProcInd
