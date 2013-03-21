module ind

    use param
    implicit none
    integer ::  I,J,K,B,IJ,IK,JK,IJK, &
                NI,NJ,NK,NB, &
                NIA,NJA,NKA,NIJA,NIJKA,NDIRA,NBLOCKA, &
                NIJ,NIJK,NIM,NJM,NKM, &
                NICV,NJCV,NKCV,NDIR,NBLOCK,NIJCV,&
                P,N,&
                ITIM,ITIMS,ITIME, & 
                IJST,IJEN,IP,IDI,IBL,IJKB,IJKP, &
                LI(NXA), LK(NZA),&
                ITB(2,NXA*NZA),JTB(2,NYA*NZA),KTB(2,NXA*NZA), &
                NEIGH(NBLOCKS,6),&
                LS,LSG

end module ind
