program grid

    implicit none
    
    integer :: N,NICV,NJCV,NKCV
    integer :: I,J,K,IJK,RES
    integer, dimension(:,:),allocatable :: A,B
    character(len=20) :: FILOUT,BLOCKNUMBER,FILEIN
    
    print *, 'ENTER RES'
    read (*,*) RES
    print *, 'ENTER N'
    read (*,*) N
    print *, 'ENTER NICV'
    read (*,*) NICV
    print *, 'ENTER NJCV'
    read (*,*) NJCV
    print *, 'ENTER NKCV'
    read (*,*) NKCV
    
    allocate(A(N,6))
    allocate(B(N,6))
    
    A=1
    B=-1

    !East and west neighbour
    do k=1,NKCV
    do i=1,NICV-1
    do j=1,NJCV
        ijk=j+NJCV*(i-1)+NJCV*NICV*(k-1)
        A(ijk,4) = 4
        B(ijk,4) = ijk+NJCV
        A(ijk+NJCV,3) = 4
        B(ijk+NJCV,3) = ijk
    enddo
    enddo
    enddo


    !North and south neighbour
    do k=1,NKCV
    do i=1,NICV
    do j=1,NJCV-1
        ijk=j+NJCV*(i-1)+NJCV*NICV*(k-1)
        A(ijk,2) = 4
        B(ijk,2) = ijk+1
        A(ijk+1,1) = 4
        B(ijk+1,1) = ijk
    enddo
    enddo
    enddo

    !Top and bottom neighbour
    do k=1,NKCV-1
    do i=1,NICV
    do j=1,NJCV
        ijk=j+NJCV*(i-1)+NJCV*NICV*(k-1)
        A(ijk,6)=4
        B(ijk,6) = ijk+NICV*NJCV
        A(ijk+NICV*NJCV,5)=4
        B(ijk+NICV*NJCV,5) = ijk
    enddo
    enddo
    enddo
    
    do K=1,NKCV
    do I=1,NICV
    do J=1,NJCV
        IJK=J+NJCV*(I-1)+NICV*NJCV*(K-1)
        write(BLOCKNUMBER,*) IJK
        FILOUT='grid_'//trim(adjustl(BLOCKNUMBER))//'.inp'
        open(UNIT=20,FILE=FILOUT)
        FILEIN='grid_'//trim(adjustl(BLOCKNUMBER))//'.pre'
        rewind 20
        write(20,*) FILEIN
        write(20,*) (-1.0d0+(I-1)*2.0d0/NICV), (-1.0d0+I*2.0d0/NICV), RES/NICV
        write(20,*) (-1.0d0+(J-1)*2.0d0/NJCV), (-1.0d0+J*2.0d0/NJCV), RES/NJCV
        write(20,*) (-1.0d0+(K-1)*2.0d0/NKCV), (-1.0d0+K*2.0d0/NKCV), RES/NKCV
        write(20,*) A(IJK,:)
        write(20,*) B(IJK,:)
        close (UNIT=20)
    enddo
    enddo
    enddo

end program grid
