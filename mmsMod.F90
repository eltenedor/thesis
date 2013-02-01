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
    
    ! stationary
    !PHI_RES = PHI_0*(dcos((Y*pi)/L)+dsin((X*pi)/L))
    

end function phi

!####################################################
function src(x,y,z,t) result(src_res)
!####################################################

    use sc
    implicit none
    
    real*8, intent(in) :: x,y,z,t
    real*8 :: src_res
    
    src_res = ALPHA*PHI_0*(-4.0d0)

end function src

!####################################################
function vel(x,y,z,t) result(vel_res)
!####################################################

    use sc
    implicit none
    
    real*8, intent(in) :: x,y,z,t
    real*8 :: vel_res

    !VEL_RES = V_0*dsin((X*pi)/L)*dsin((Y*pi)/L)
    vel_res = 0.0d0
    
    
end function vel

end module mms
