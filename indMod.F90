module ind

    use param
    implicit none
    integer ::  I,J,K,IJ,IK,JK,IJK, &
                NI,NJ,NK, &
                NIJ,NIJK,NIM,NJM,NKM, &
                NICV,NJCV,NKCV,NWALI,NIJCV,&
                L, N,&
                ITIM,ITIMS,ITIME, & 
                IJST, IJEN, IP, IW, IJKB, IJKP, &
                LI(NXA), LK(NZA),&
                ITB(2,NXA*NZA),JTB(2,NYA*NZA),KTB(2,NXA*NZA), &
                ! PETSc routine indices start counting from 0!!
                CTD(0:NXYZA-1),DTC(NXYZA),&
                LS,LSG

end module ind
