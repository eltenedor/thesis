!########################################################
module indexModule
!########################################################

    use parameterModule
    implicit none
    integer ::  &
                ! Time stepping indices
                ITIM,ITIMS,ITIME,&
                ! Global Size Variables
                NIA,NJA,NKA,NIJA,NIJKA,NDIRA,NNEUA,NWALA,NBLOA,NA,NPROCSA,&
                ! Regular Indices block independent (IJKB needed?)
                I,J,K,IJ,IK,JK,IJK,NI,NIM,NJ,NJM,NK,NKM,NIJ,NIJK,NICV,NJCV,NKCV,NIJCV,N,IJKDIR,IJKNEU,IJKWAL,IJKB,IJKBLO,IJKP,&
                IJKSTL,IJKENL,IJKSTR,IJKENR,&
                ! Indices for complete Blocks, Processor dependent
                B,BB,NB,&
                NIJKPROC,NDIRPROC,NNEUPROC,NWALPROC,NBLOPROC,NFACEPROC,&
                ! Indices inside Blocks
                IJKL,IJKR,&
                ! Neighbour Block Indices, only needed in grgen and preprocessing
                INEIGH,NEIGH(NBLOCKS,6),&
                ! Indices storing boundary type, only needed in grgen
                P,ITB(2,NXA*NZA),JTB(2,NYA*NZA),KTB(2,NXA*NZA),&
                ! block dependent regular indices
                IST,IBL(NBLOCKS),NIBL(NBLOCKS),&
                JST,JBL(NBLOCKS),NJBL(NBLOCKS),&
                KST,KBL(NBLOCKS),NKBL(NBLOCKS),&
                IJKST,IJKBL(NBLOCKS),NIJKBL(NBLOCKS),NBL(NBLOCKS),&
                !NICVBL(NBLOCKS),NJCVBL(NBLOCKS),NKCVBL(NBLOCKS),&
                !LI(1000),LK(1000),&
                ! block dependent indices (inlet boundary - DIR)
                IJKDIRST,IJKDIRBL(NBLOCKS),NDIRBL(NBLOCKS),NDIR,&
                ! block dependent indices (outlet boundary - NEU)
                IJKNEUST,IJKNEUBL(NBLOCKS),NNEUBL(NBLOCKS),NNEU,&
                ! block dependent indices (wall boundary - WAL)
                IJKWALST,IJKWALBL(NBLOCKS),NWALBL(NBLOCKS),NWAL,&
                ! block dependent indices (block boundary - BLO, left - L, right - R)
                IJKBLOBL(NBLOCKS),NBLOBL(NBLOCKS),&
                IJKBLOST,NBLO,IJKBLOSTL,NBLOL,NIJKL,IJKBLOSTR,NBLOR,NIJKR,&
                IJKMARKL,IJKMARKR,&
                ! block dependent indices (boundary faces)
                F,NF,&
                FACEST,FACEBL(NBLOCKS),NFACEBL(NBLOCKS),NFACE,&
                ! mapping related indices
                IJK_GLO,IJK_GLOBL(NBLOCKS),IJK_LOC,MIJK(NXYZA*NPROCS),RMIJK(0:NAL-1),&
                B_GLO(NBLOCKS),& !NBL(NBLOCKS),
                IJKPROC,IJKPROC_GLO,&
                ! indices for outer iterations
                LS,LSG
               ! integer values used in preprocessing to determine maximal
               ! processor load
    integer :: NXMAX,NYMAX,NZMAX,NXYZMAX,&
               NDIRMAX,NNEUMAX,NWALMAX,NBLOMAX,&
               NBLOCKSMAX,NFMAX
                
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

    IJKBLOSTL=IJKBLOBL(BL)
    NBLOL=NBLOBL(BL)
    NIJKL=NIJKBL(BL)
    
    IJKBLOSTR=IJKBLOBL(BR)
    NBLOR=NBLOBL(BR)
    NIJKR=NIJKBL(BR)

end subroutine setBlockInd2Int

!########################################################
subroutine setBlockInd1Int(B)
!########################################################

    use geoModule, only : DX,DY,DZ,VOL,DXBL,DYBL,DZBL,VOLBL
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

    IJKDIRST=IJKDIRBL(B)
    NDIR=NDIRBL(B)
    
    IJKNEUST=IJKNEUBL(B)
    NNEU=NNEUBL(B)

    IJKWALST=IJKWALBL(B)
    NWAL=NWALBL(B)
    
    IJKBLOST=IJKBLOBL(B)
    NBLO=NBLOBL(B)
    
    FACEST=FACEBL(B)
    NFACE=NFACEBL(B)

    ! GEOMETRY RELATED DATA
    DX=DXBL(B)
    DY=DYBL(B)
    DZ=DZBL(B)
    VOL=VOLBL(B)

end subroutine setBlockInd1Int

end module indexModule
