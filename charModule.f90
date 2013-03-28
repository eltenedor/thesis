module charModule

    implicit none
    character(len=20) :: FILIN, FILOUT,VTKFILE,&
                        BLOCKFILE,BLOCK_CH,PROCFILE,PROC_CH,&
                        TIME_CH
    integer :: BLOCKUNIT,PROC,PROCUNIT,OFFSET

end module charModule
