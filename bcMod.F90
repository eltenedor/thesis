module bc

    use param
    implicit none
    integer ::  IJKDI(NDIRAL), IJKPDI(NDIRAL), IJKDI1(NDIRAL), IJKDI2(NDIRAL), IJKDI3(NDIRAL), IJKDI4(NDIRAL), &
                IJKBL(NBLOCKAL), IJKPBL(NBLOCKAL), IJKBL1(NBLOCKAL), IJKBL2(NBLOCKAL), IJKBL3(NBLOCKAL), IJKBL4(NBLOCKAL), &
                BTYP(6)
    real(KIND=PREC) :: SRDDI(NDIRAL)

end module bc
