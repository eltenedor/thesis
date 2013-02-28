!#########################################################
module geo
!#########################################################

    use param
    implicit none
    real(KIND=PREC) ::  X(NXYZA),Y(NXYZA),Z(NXYZA), &
                        XC(NXYZA),YC(NXYZA),ZC(NXYZA), &
                        ! Anstatt der 100 ein variabler Maximalwert oder allocatable
                        XL(100),YL(100),ZL(100),XR(100),YR(100),ZR(100),&
                        XF(100),YF(100),ZF(100),&
                        XXS,XXE,YYS,YYE,ZZS,ZZE, &
                        DX,DY,DZ, &
                        AR,DN,NX,NY,NZ,&
                        VOL,XPN,YPN,ZPN,SX,SY,SZ,&
                        FX(NXYZA),FY(NXYZA),FZ(NXYZA)
                        
contains

!#########################################################
subroutine normalArea(IJKP,IJKN,IJK2,IJK3,IJK4,ARR,DNN,XPNN,YPNN,ZPNN,NXX,NYY,NZZ)
!#########################################################
! THIS ROUTINE CALCULATES THE AREA SPANNED BY TWO VECTORS
! AND THE COMPONENTS OF THE RESPECTIVE NORMAL VECTOR USING
! A CROSS PRODUCT
    
    implicit none
    integer, intent(in) :: IJKP,IJKN,IJK2,IJK3,IJK4
    real(KIND=PREC), intent(in out) :: ARR,DNN,XPNN,YPNN,ZPNN,NXX,NYY,NZZ
    real(KIND=PREC) :: X2,X3,X4,Y2,Y3,Y4,Z2,Z3,Z4
    
    ! BOUNDARY NODES OF CELL FACE
    X2=X(IJK2)
    Y2=Y(IJK2)
    Z2=Z(IJK2)
    !
    X3=X(IJK3)
    Y3=Y(IJK3)
    Z3=Z(IJK3)
    !
    X4=X(IJK4)
    Y4=Y(IJK4)
    Z4=Z(IJK4)
    !
    ! CROSS PRODUCT
    !
    NXX = -(Y2-Y3)*(Z3-Z4)+(Y3-Y4)*(Z2-Z3)
    NYY = (X2-X3)*(Z3-Z4)-(X3-X4)*(Z2-Z3)
    NZZ = -(X2-X3)*(Y3-Y4)+(X3-X4)*(Y2-Y3)
    !
    ARR = SQRT(NXX**2+NYY**2+NZZ**2)
    !            
    ! UNITY NORMAL VECTOR OF FACE
    !
    NXX = NXX/AR
    NYY = NYY/AR
    NZZ = NZZ/AR
    !
    ! CALCULATE NORMAL DISTANCE BETWEEN ADJACENT CVS
    !              
    XPNN=XC(IJKN)-XC(IJKP)
    YPNN=YC(IJKN)-YC(IJKP)
    ZPNN=ZC(IJKN)-ZC(IJKP)
    !
    DNN=SQRT(XPNN**2+YPNN**2+ZPNN**2)
    !
end subroutine normalArea

end module geo
