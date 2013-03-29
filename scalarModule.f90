module scalarModule

    use parameterModule
    implicit none
    real(KIND=PREC) :: RHO = 1.0d0, ALPHA = 1.0d0, &
        VX = 1.0d0, VY = 1.0d0, VZ = 1.0d0,V_0=0.0d0, &
        T_0 = 1.0d0, PHI_0 = 1.0d0, &
        DT = 0.015625d0,pi=3.1415926535897932d0,&
        SMALL = 1.0E-20, ZERO = 0.0d0

end module scalarModule
