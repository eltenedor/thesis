module geo

    use param
    implicit none
    real(KIND=PREC) ::  X(NXYZA),Y(NXYZA),Z(NXYZA), &
                        XC(NXYZA),YC(NXYZA),ZC(NXYZA), &
                        XXS,XXE,YYS,YYE,ZZS,ZZE, &
                        XN,YN,ZN, &
                        NX,NY,NZ, &
                        DX,DY,DZ, &
                        VOL,XPN,YPN,SX,SY,&
                        FX(NXYZA),FY(NXYZA),FZ(NXYZA)

end module geo
