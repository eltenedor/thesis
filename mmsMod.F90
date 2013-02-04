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

    PHI_RES = X**2+Y**2

end function phi

!####################################################
function src(x,y,z,t) result(src_res)
!####################################################

    use sc
    implicit none
    
    real*8, intent(in) :: x,y,z,t
    real*8 :: src_res
    
    SRC_RES = x*2.0d0+y*2.0d0-4.0d0


end function src

!####################################################
function vel(x,y,z,t) result(vel_res)
!####################################################

    use sc
    implicit none
    
    real*8, intent(in) :: x,y,z,t
    real*8 :: vel_res

    
    VEL_RES = sin(X*pi)*sin(Y*pi)
    !vel_res = 0.0d0
    
    
end function vel

!####################################################
function velx(x,y) result(velx_res)
!####################################################

    use sc
    implicit none
    
    real*8, intent(in) :: x,y
    real*8 :: velx_res
    
    VELX_RES = -y/(x**2+y**2)
    
end function velx

!####################################################
function vely(x,y) result(vely_res)
!####################################################

    use sc
    implicit none
    
    real*8, intent(in) :: x,y
    real*8 :: vely_res
    
    VELY_RES = x/(x**2+y**2)
    
end function vely



end module mms
