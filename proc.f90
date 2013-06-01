program proc

    implicit none

    integer :: NPROCS,NBLOCKS,I,J,IBLOCK
    character(len=20) :: FILOUT,BLOCKNUMBER,PROCNUMBER

    print *, 'ENTER NUMBER OF BLOCKS'
    read (*,*) NBLOCKS
    print *, 'ENTER NUMBER OF PROCS'
    read (*,*) NPROCS

    IBLOCK=0
    do I=0,NPROCS-1
        write(PROCNUMBER,*) I
        FILOUT='proc_'//trim(adjustl(PROCNUMBER))//'.inp'
        open(UNIT=20,FILE=FILOUT)
        rewind 20
        write(20,*) NBLOCKS/NPROCS
        write(20,*) (J,J=IBLOCK,IBLOCK+NBLOCKS/NPROCS)
        IBLOCK=IBLOCK+NBLOCKS/NPROCS+1
    end do

end program proc
