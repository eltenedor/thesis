!########################################################
module indMod
!########################################################

    use param
    implicit none
    integer ::  &
                ! Time stepping indices
                ITIM,ITIMS,ITIME,&
                ! Regular Indices block independent (IJKB needed?)
                I,J,K,IJK,NI,NIM,NJ,NJM,NK,NKM,NIJ,NIJK,NICV,NJCV,NKCV,NIJCV,N,IJKDIR,IJKB,IJKBLOCK,IJKP,&
                IJKSTL,IJKEL,IJKSTR,IJKER,&
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
                !NICVBL(NBLOCKS),NJCVBL(NBLOCKS),NKCVBL(NBLOCKS),&
                !LI(1000),LK(1000),&
                ! block dependent indices (dirichlet boundary)
                IJKDIRST,IJKDIRBL(NBLOCKS),NDIRBL(NBLOCKS),NDIR,&
                ! block dependent indices (block boundary)
                IJKBLOCKBL(NBLOCKS),NBLOCKBL(NBLOCKS),&
                IJKBLOCKST,NBLOCK,IJKBLOCKSTL,NBLOCKL,IJKBLOCKSTR,NBLOCKR,&
                IJKMARKL,IJKMARKR,&
                ! block dependent indices (boundary faces)
                F,NF,&
                FACEST,FACEBL(NBLOCKS),NFACEBL(NBLOCKS),NFACE,&
                ! mapping related indices
                IJK_GLO,IJK_LOC,MIJK(NXYZA),&
                B_GLO(NBLOCKS),& !NBL(NBLOCKS),
                IJKPROC,&
                ! indices for outer iterations
                LS,LSG
                
                
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

    use geo, only : DX,DY,DZ,VOL,DXBL,DYBL,DZBL,VOLBL
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
    NIJ=NI*NJ
    NICV=NIM-1
    NJCV=NJM-1
    NKCV=NKM-1
    NIJCV=NICV*NJCV

    IJKST=IJKBL(B)
    NIJK=NIJKBL(B)

    IJKBLOCKST=IJKBLOCKBL(B)
    NBLOCK=NBLOCKBL(B)
    
    IJKDIRST=IJKDIRBL(B)
    NDIR=NDIRBL(B)
    
    FACEST=FACEBL(B)
    NFACE=NFACEBL(B)

    ! GEOMETRY RELATED DATA
    DX=DXBL(B)
    DY=DYBL(B)
    DZ=DZBL(B)
    VOL=VOLBL(B)

end subroutine setBlockInd1Int

end module indMod
