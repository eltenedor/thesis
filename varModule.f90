! this module contains the arrays for storing the scalar T, its respective
! cell center gradient components and the mass fluxes F1...3
module varModule

    use parameterModule
    implicit none

    real(KIND=PREC) ::  T(NXYZA),TO(NXYZA),TR(NFACEAL),TIME,&
                        DTX(NXYZA),DTY(NXYZA),DTZ(NXYZA),&
                        DTXR(NFACEAL),DTYR(NFACEAL),DTZR(NFACEAL),&
                        F1(NXYZA),F2(NXYZA),F3(NXYZA)

end module varModule
