module bc

    use param
    implicit none
    integer ::  IJW(NWA), IJPW(NWA), IJW1(NWA), IJW2(NWA), &
                BTYP(4)
    real(KIND=PREC) :: XTW(NWA), YTW(NWA), SRDW(NWA), FWT(NWA)

end module bc
