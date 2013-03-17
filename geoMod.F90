!#########################################################
module geo
!#########################################################

    use param
    implicit none
    real(KIND=PREC) ::  X(NXYZA),Y(NXYZA),Z(NXYZA), &
                        XC(NXYZA),YC(NXYZA),ZC(NXYZA), &
                        ! Anstatt der 100 ein variabler Maximalwert oder allocatable
                        XYZL(12),XYZR(12),XYZCommon(12),XYZF(3),&
                        XXS,XXE,YYS,YYE,ZZS,ZZE, &
                        DX,DY,DZ, &
                        AR,DN,NX,NY,NZ,&
                        VOL,XPN,YPN,ZPN,SX,SY,SZ,&
                        FX(NXYZA),FY(NXYZA),FZ(NXYZA)
                        
    public :: normalArea,calcGrad
    private :: normalAreaIndex, normalAreaFace
    
    interface normalArea
        module procedure normalAreaIndex, normalAreaFace
    end interface
                        
contains

!#########################################################
subroutine normalAreaIndex(IJKP,IJKN,IJK2,IJK3,IJK4,ARR,DNN,XPNN,YPNN,ZPNN,NXX,NYY,NZZ)
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
    ! CALCULATE DISTANCE VECTOR BETWEEN ADJACENT CVS
    !              
    XPNN=XC(IJKN)-XC(IJKP)
    YPNN=YC(IJKN)-YC(IJKP)
    ZPNN=ZC(IJKN)-ZC(IJKP)
    !
    DNN=SQRT(XPNN**2+YPNN**2+ZPNN**2)
    !
end subroutine normalAreaIndex

!#########################################################
subroutine normalAreaFace(XYZC,L,R,ARR,DNN,XPNN,YPNN,ZPNN,NXX,NYY,NZZ)
!#########################################################
! THIS ROUTINE CALCULATES THE AREA SPANNED BY TWO VECTORS
! AND THE COMPONENTS OF THE RESPECTIVE NORMAL VECTOR USING
! A CROSS PRODUCT
    
    implicit none
    integer,intent(in) :: L,R
    real(KIND=PREC), intent(in) :: XYZC(12)
    real(KIND=PREC), intent(in out) :: ARR,DNN,XPNN,YPNN,ZPNN,NXX,NYY,NZZ
    real(KIND=PREC) :: X2,X3,X4,Y2,Y3,Y4,Z2,Z3,Z4
    
    ! BOUNDARY NODES OF CELL FACE
    X2=XYZC(1)
    Y2=XYZC(2)
    Z2=XYZC(3)
    !
    X3=XYZC(4)
    Y3=XYZC(5)
    Z3=XYZC(6)
    !
    X4=XYZC(7)
    Y4=XYZC(8)
    Z4=XYZC(9)
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
    ! CALCULATE DISTANCE VECTOR BETWEEN ADJACENT CVS
    !              
    XPNN=abs(XC(R)-XC(L))
    YPNN=abs(YC(R)-YC(L))
    ZPNN=abs(ZC(R)-ZC(L))
    !
    DNN=SQRT(XPNN**2+YPNN**2+ZPNN**2)
    !
end subroutine normalAreaFace

!#########################################################
subroutine calcGrad(L,R,XFF,YFF,ZFF,FF)
!#########################################################

    implicit none
    integer, intent(in) :: L,R
    real(KIND=PREC), intent(in) :: XFF,YFF,ZFF
    real(KIND=PREC), intent(inout) :: FF
    real(KIND=PREC) :: DLPN,DLNN,SMALL

    SMALL=1.e-20

    DLPN=sqrt((XFF-XC(L))**2+(YFF-YC(L))**2+(ZFF-ZC(L))**2)
    DLNN=sqrt((XC(R)-XFF)**2+(YC(R)-YFF)**2+(ZC(R)-ZFF)**2)

    FF=DLPN/(DLPN+DLNN+SMALL)
    
end subroutine calcGrad

!#########################################################
subroutine reverseOrder(XYZC)
!#########################################################
! Used to reverse order of face vertices
! needed for correct calculation of normal vector

    implicit none
    real(KIND=PREC), intent(inout) :: XYZC(12)
    real(KIND=PREC) :: XYZTEMP(12)
    
    XYZTEMP=XYZC
    XYZC(1:3)=XYZTEMP(4:6)
    XYZC(4:6)=XYZTEMP(1:3)
    XYZC(7:9)=XYZTEMP(10:12)
    XYZC(10:12)=XYZTEMP(7:9)
    
end subroutine reverseOrder
    
end module geo
