module bc

    use param
    implicit none
    integer ::  IJKD(NDIR), IJKPD(NDIR), IJKD1(NDIR), IJKD2(NDIR), IJKD3(NDIR), IJKD4(NDIR), &
                IJKB(NBLOCK), IJKPB(NBLOCK), IJKB1(NBLOCK), IJKB2(NBLOCK), IJKB3(NBLOCK), IJKB4(NBLOCK), &
                BTYP(6)
    real(KIND=PREC) :: SRDD(NDIR),SRDB(NBLOCK)

end module bc
