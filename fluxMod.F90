module flux
    
    use param
    implicit none
    
    real(KIND=PREC) ::  VSOL,FACP,G,&
                        FII,DFXI,DFYI,DFZI,&
                        FCFIE,FDFIE,FCFII,FDFII,FFIC
                        
    real(KIND=PREC),parameter :: ZERO=0.0d0

end module flux
