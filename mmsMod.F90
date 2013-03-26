!#########################################################
module mms
!#########################################################

    use logic, only: LTIME
    implicit none
    
contains

!####################################################
function phi(x,y,z,t) result(phi_res)
!####################################################

    use sc
    implicit none
    
    real*8, intent(in) :: x, y, z, t
    real*8 :: phi_res 

    if (LTIME) then
        PHI_RES = exp(-t)*(x**2+y**2+z**2)
    else
        PHI_RES = X**2+Y**2+Z**2
    end if

end function phi

!####################################################
function src(x,y,z,t) result(src_res)
!####################################################

    use sc
    implicit none
    
    real*8, intent(in) :: x,y,z,t
    real*8 :: src_res
    
    if (LTIME) then
        SRC_RES = exp(-t)*(-6)+x*exp(-t)*2+y*exp(-t)*2+z*exp(-t)*2-exp(-t) &
                    *(x**2+y**2+z**2)
    else
        SRC_RES = (X+Y+Z)*2.0d0-6.0d0
    end if

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
