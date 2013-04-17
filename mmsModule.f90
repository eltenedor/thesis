!#########################################################
module mmsModule
!#########################################################

    use controlModule, only: LTIME
    implicit none
    
contains

!####################################################
function phi(x,y,z,t) result(phi_res)
!####################################################

    use scalarModule
    implicit none
    
    real*8, intent(in) :: x, y, z, t
    real*8 :: phi_res 

    if (LTIME) then
        PHI_RES = exp(-t)*(x**2+y**2+z**2)
        !PHI_RES = exp(-t)*(cos(pi*y**2)*3+sin(pi*x**2)*2+sin(pi*z**2)*5)
    else
        !PHI_RES = x**2.0d0*2.0d0+y**2.0d0*3.0d0+z**2.0d0*5.0d0
        !PHI_RES = x**2+y**2+z**2
        PHI_RES = cos(pi*y**2)*3+sin(pi*x**2)*2+sin(pi*z**2)*5
    end if

end function phi

!####################################################
function src(x,y,z,t) result(src_res)
!####################################################

    use scalarModule
    implicit none
    
    real*8, intent(in) :: x,y,z,t
    real*8 :: src_res
    
    if (LTIME) then
        SRC_RES = exp(-t)*(-6)+x*exp(-t)*2+y*exp(-t)*2+z*exp(-t)*2-exp(-t) &
                    *(x**2+y**2+z**2)
        !SRC_RES = -exp(-t)*(cos(pi*y**2)*3+sin(pi*x**2)*2+sin(pi*z**2)*5)- &
        !pi*cos(pi*x**2)*exp(-t)*4-pi*cos(pi*z**2)*exp(-t)*10+pi*sin(pi*y**   &
        !2)*exp(-t)*6+pi**2*y**2*cos(pi*y**2)*exp(-t)*12+pi**2*x**2*sin(pi*   &
        !x**2)*exp(-t)*8+pi**2*z**2*sin(pi*z**2)*exp(-t)*20+pi*x*cos(pi*x**   &
        !2)*exp(-t)*4+pi*z*cos(pi*z**2)*exp(-t)*10-pi*y*sin(pi*y**2)*exp(-t   &
        !)*6
    else
        !SRC_RES = x*4.0d0+y*6.0d0+z*10.0d0-20.0d0
        !SRC_RES = (X+Y+Z)*2.0d0-6.0d0
        SRC_RES = pi*cos(pi*x**2)*(-4)-pi*cos(pi*z**2)*10+pi*sin(pi*y**2)* &
            6+pi**2*y**2*cos(pi*y**2)*12+pi**2*x**2*sin(pi*x**2)*8+pi**2*z**2* &
            sin(pi*z**2)*20+pi*x*cos(pi*x**2)*4+pi*z*cos(pi*z**2)*10-pi*y*sin( &
            pi*y**2)*6
    end if

end function src

end module mmsModule
