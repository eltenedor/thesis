program proc

    implicit none

    integer :: NPROCS,NBLOCKS,I,J,IBLOCK,STRIDE
    character(len=20) :: FILOUT,BLOCKNUMBER,PROCNUMBER

    print *, 'ENTER NUMBER OF BLOCKS'
    read (*,*) NBLOCKS
    print *, 'ENTER NUMBER OF PROCS'
    read (*,*) NPROCS

    STRIDE=NBLOCKS/NPROCS

    IBLOCK=0
    do I=0,NPROCS-1
        print *, I
        write(PROCNUMBER,*) I
        FILOUT='proc_'//trim(adjustl(PROCNUMBER))//'.inp'
        open(UNIT=20,FILE=FILOUT)
        rewind 20
        write(20,*) STRIDE
        print *, IBLOCK,STRIDE,(J,J=IBLOCK+1,IBLOCK+STRIDE)
        write(20,*) (J,J=IBLOCK,IBLOCK+STRIDE)
        IBLOCK=IBLOCK+STRIDE
    end do

end program proc
