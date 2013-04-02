module boundaryModule

    use parameterModule
    implicit none
    integer ::  IJKBDIR(NDIRAL), IJKPDIR(NDIRAL), IJKDIR1(NDIRAL), IJKDIR2(NDIRAL), IJKDIR3(NDIRAL), IJKDIR4(NDIRAL)
    integer ::  IJKBNEU(NNEUAL), IJKPNEU(NNEUAL), IJKNEU1(NNEUAL), IJKNEU2(NNEUAL), IJKNEU3(NNEUAL), IJKNEU4(NNEUAL)
    integer ::  IJKBWAL(NWALAL), IJKPWAL(NWALAL), IJKWAL1(NWALAL), IJKWAL2(NWALAL), IJKWAL3(NWALAL), IJKWAL4(NWALAL)
    
    !
    ! block boundary related data (only needed for grid generation and preprocessing)
    !
    integer ::  IJKBBLO(NBLOAL), IJKPBLO(NBLOAL), IJKBLO1(NBLOAL), IJKBLO2(NBLOAL), IJKBLO3(NBLOAL), IJKBLO4(NBLOAL)

    integer ::  BTYP(6),&
    !
    ! face index related data
    !
                L(NFACEAL),R(NFACEAL)
                
    !real(KIND=PREC) :: SRDDIR(NDIRAL)
    !real(KIND=PREC) :: SRDNEU(NNEUAL)
    real(KIND=PREC) :: SRDWAL(NWALAL)
    !
    ! face geometry related data
    !
    real(KIND=PREC) ::  XF(NFACEAL),YF(NFACEAL),ZF(NFACEAL),&
                        FF(NFACEAL),ARF(NFACEAL),DNF(NFACEAL),&
                        XPNF(NFACEAL),YPNF(NFACEAL),ZPNF(NFACEAL),&
                        NXF(NFACEAL),NYF(NFACEAL),NZF(NFACEAL)

end module boundaryModule
