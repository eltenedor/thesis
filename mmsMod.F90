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
    PHI_RES = PHI_0*(X**2+Y**2)

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
    SRC_RES = ALPHA*PHI_0*(-4)+PHI_0*RHO*V_0*X*2+PHI_0*RHO*V_0*Y*2


    
end function src

end module mms
