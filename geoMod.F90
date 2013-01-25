module geo

    use param
    implicit none
    real(KIND=PREC) :: X(NXYA),Y(NXYA),XC(NXYA),YC(NXYA),&
            XXS,XXE,YYS,YYE, &
            DX,DY,VOL,XPN,YPN,SX,SY,&
            FX(NXYA),FY(NXYA)

end module geo
