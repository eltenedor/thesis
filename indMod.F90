module ind

    use param
    implicit none
    integer ::  I,J,K,IJK, &
                NI,NJ,NK, &
                NIJ,NIJK,NIM,NJM,NKM, &
                NICV,NJCV,NKCV,NWALI,&
                L, N,&
                ITIM,ITIMS,ITIME, & 
                IJST, IJEN, IP, IW, IJKB, IJKP, &
                LI(NXA), LK(NZA),&
                ITB(2,NXA),JTB(2,NYA),KTB(2,NZA), &
                ! PETSc routine indices start counting from 0!!
                CTD(0:NXYZA-1),DTC(NXYZA),&
                LS,LSG

end module ind
