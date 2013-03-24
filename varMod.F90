module var

    use param
    implicit none

    real(KIND=PREC) ::  T(NXYZA),TO(NXYZA),TIME,&
                        DTX(NXYZA),DTY(NXYZA),DTZ(NXYZA),&
                        F1(NXYZA),F2(NXYZA),F3(NXYZA)

end module var
