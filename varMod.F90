module var

    use param
    implicit none

    real(KIND=PREC) ::  T(NXYZA),TO(NXYZA),TR(NFACEAL),TIME,&
                        DTX(NXYZA),DTY(NXYZA),DTZ(NXYZA),&
                        DTXR(NFACEAL),DTYR(NFACEAL),DTZR(NFACEAL),&
                        F1(NXYZA),F2(NXYZA),F3(NXYZA)

end module var
