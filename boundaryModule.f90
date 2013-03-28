module boundaryModule

    use parameterModule
    implicit none
    integer ::  IJKBDI(NDIRAL), IJKPDI(NDIRAL), IJKDI1(NDIRAL), IJKDI2(NDIRAL), IJKDI3(NDIRAL), IJKDI4(NDIRAL), &
    
    !
    ! block boundary related data (only needed for grid generation and preprocessing)
    !
                IJKBBL(NBLOCKAL), IJKPBL(NBLOCKAL), IJKBL1(NBLOCKAL), IJKBL2(NBLOCKAL), IJKBL3(NBLOCKAL), IJKBL4(NBLOCKAL), &
                BTYP(6),&
    !
    ! face index related data
    !
                L(NFACEAL),R(NFACEAL)
                
    real(KIND=PREC) :: SRDDI(NDIRAL)
    !
    ! face geometry related data
    !
    real(KIND=PREC) ::  XF(NFACEAL),YF(NFACEAL),ZF(NFACEAL),&
                        FF(NFACEAL),ARF(NFACEAL),DNF(NFACEAL),&
                        XPNF(NFACEAL),YPNF(NFACEAL),ZPNF(NFACEAL),&
                        NXF(NFACEAL),NYF(NFACEAL),NZF(NFACEAL)

end module boundaryModule
