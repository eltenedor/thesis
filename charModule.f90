module charModule

    implicit none
    character(len=20) :: FILIN, FILOUT,VTKFILE,&
                        BLOCKFILE,BLOCK_CH,PROCFILE,PROC_CH,&
                        TIME_CH
    character(len=11) :: FORM_CH='UNFORMATTED'
    integer :: BLOCKUNIT,PROC,PROCUNIT
    integer, parameter :: OFFSET=103

end module charModule
