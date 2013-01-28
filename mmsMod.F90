module mms
    implicit none
    
contains

!####################################################
function phi(x,y,z,t) result(phi_res)
!####################################################

    use sc
    implicit none
    
    real*8, intent(in) :: x, y, z, t
    real*8 :: phi_res 

	! instationary
    !PHI_RES = PHI_0*exp(-(T-T_0)/T_0)*(X**2+Y**2)
    
    ! stationary
    !PHI_RES = PHI_0*(X**2+Y**2)
    
    ! stationary
    PHI_RES = PHI_0*(cos((Y*pi)/L)+sin((X*pi)/L))
    

end function phi

!####################################################
function src(x,y,z,t) result(src_res)
!####################################################

    use sc
    implicit none
    
    real*8, intent(in) :: x,y,z,t
    real*8 :: src_res
    
    ! instationary
    !SRC_RES = -(PHI_0*exp(-(T-T_0)/T_0)*(RHO*X**2+RHO*Y**2+ALPHA*T_0*4 &
    ! &-RHO*T_0*V_0*X*2-RHO*T_0*V_0*Y*2))/T_0
    
    ! stationary
    !SRC_RES = ALPHA*PHI_0*(-4)+PHI_0*RHO*V_0*X*2+PHI_0*RHO*V_0*Y*2
    
    ! stationary
    SRC_RES = 1/L**2*PHI_0*pi*(ALPHA*pi*cos((Y*pi)/L)*4+ALPHA*pi*sin((&
     &X*pi)/L)*4+L*RHO*V_0*cos((pi*(X*2-Y))/L)+L*RHO*V_0*cos((Y*pi)/L)*2&
     &-L*RHO*V_0*cos((pi*(X*2+Y))/L)*3+L*RHO*V_0*sin((pi*(X-Y*2))/L)+L*R&
     &HO*V_0*sin((pi*(X+Y*2))/L)*3)*(1.0D0/4.0D0)

end function src

!####################################################
function vel(x,y,z,t) result(vel_res)
!####################################################

    use sc
    implicit none
    
    real*8, intent(in) :: x,y,z,t
    real*8 :: vel_res

    VEL_RES = V_0*sin((X*pi)/L)*sin((Y*pi)/L)
    
end function vel

end module mms
