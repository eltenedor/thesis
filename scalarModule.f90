module scalarModule

    use parameterModule
    implicit none
    real(KIND=PREC) :: RHO = 1.0d0, ALPHA = 1.0d0 ! material constants
    real(KIND=PREC) :: VX = 1.0d0, VY = 1.0d0, VZ = 1.0d0,V_0=0.0d0 ! velocity
    real(KIND=PREC) :: T_0 = 0.0d0,DT = 0.125d0
    !DT=1.0d0
    !DT=0.5d0
    !DT=0.25d0
    !DT=0.125d0
    !DT=0.0625d0
    !DT=0.0625d0/3.8258604184154104d0
    !DT=0.03125d0
    !DT=0.015625d0
    !DT=0.0078125d0
    !DT=0.00390625d0
    !DT=0.000244140625d0
    real(KIND=PREC) :: PHI_0 = 1.0d0
    real(KIND=PREC) :: pi=3.1415926535897932d0,SMALL=1.0E-20,ZERO=0.0d0,ONE=1.0d0 ! mathematical constants

end module scalarModule
