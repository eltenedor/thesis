module boundaryModule

    use parameterModule
    implicit none
    integer ::  IJKBINL(NINLAL), IJKPINL(NINLAL), IJKINL1(NINLAL), IJKINL2(NINLAL), IJKINL3(NINLAL), IJKINL4(NINLAL)
    integer ::  IJKBOUT(NOUTAL), IJKPOUT(NOUTAL), IJKOUT1(NOUTAL), IJKOUT2(NOUTAL), IJKOUT3(NOUTAL), IJKOUT4(NOUTAL)
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
                
    real(KIND=PREC) :: SRDINL(NINLAL)
    real(KIND=PREC) :: SRDOUT(NOUTAL)
    real(KIND=PREC) :: SRDWAL(NWALAL)
    !
    ! face geometry related data
    !
    real(KIND=PREC) ::  XF(NFACEAL),YF(NFACEAL),ZF(NFACEAL),&
                        FF(NFACEAL),ARF(NFACEAL),DNF(NFACEAL),&
                        XPNF(NFACEAL),YPNF(NFACEAL),ZPNF(NFACEAL),&
                        NXF(NFACEAL),NYF(NFACEAL),NZF(NFACEAL)

end module boundaryModule
