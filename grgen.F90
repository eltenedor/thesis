!########################################################
program grgen
!########################################################

    use ind
    implicit none

    print *, ' TOTAL NUMBER OF BLOCKS: '
    read(*,*) NB

    ! Create one grid.out and grid.vtk for each block
    do B=1,NB
        call readData
        call cartesian
        call setBc
        call calcG
        call gridExport
    end do
    call writeParamMod

end program grgen

!########################################################
subroutine readData
!########################################################

    use bc
    use ch
    use geo
    use ind
    implicit none

    integer :: ITYP
    
    PRINT *, ' INPUT FILE NAME (* - KEYBOARD):  '
    READ(*,1) FILIN
  1 FORMAT(A12)
    IF(FILIN.NE.'*') THEN
        OPEN (UNIT=2,FILE=FILIN)
        REWIND 2
        ITYP=0
    ELSE
        OPEN (UNIT=1,FILE='grid.inp')
        REWIND 1
        ITYP=1
    ENDIF

    PRINT *, ' OUTPUT FILE NAME:  '
    IF(ITYP.EQ.1) THEN
        READ(*,1) FILOUT
        WRITE(1,1) FILOUT
    ELSE
        READ(2,1) FILOUT
    ENDIF

    OPEN (UNIT=3,FILE=FILOUT)
    REWIND 3

    PRINT *, ' ENTER> XSTART, XEND, NUMBER OF CVS:  '
    IF(ITYP.EQ.1) THEN
        READ(*,*) XXS,XXE,NICV
        WRITE(1,*) XXS,XXE,NICV,'   XS,XE,NICV '
    ELSE
        READ(2,*) XXS,XXE,NICV
    END IF

    PRINT *, ' ENTER> YSTART, YEND, NUMBER OF CVS:  '
    IF(ITYP.EQ.1) THEN
        READ(*,*) YYS,YYE,NJCV
        WRITE(1,*) YYS,YYE,NJCV,'   YS,YE,NJCV '
    ELSE
        READ(2,*) YYS,YYE,NJCV
    END IF

    PRINT *, ' ENTER> ZSTART, ZEND, NUMBER OF CVS:  '
    IF(ITYP.EQ.1) THEN
        READ(*,*) ZZS,ZZE,NKCV
        WRITE(1,*) ZZS,ZZE,NKCV,'   ZS,ZE,NKCV '
    ELSE
        READ(2,*) ZZS,ZZE,NKCV
    END IF

    PRINT *, ' ENTER> BOUNDARY TYPE S N W E B T:  (1 - DIRICHLET, 2 - BLOCK)'
    IF(ITYP.EQ.1) THEN
        READ(*,*) BTYP(1:6)
        WRITE(1,*) BTYP(1:6),  '   SBTYP, NBTYP, WBTYP, EBTYP,  BBTYP ,TBTYP '
    ELSE
        READ(2,*) BTYP(1:6)
    END IF
    
    PRINT *, ' ENTER> NEIGHBOUR INDEX N S W E T B:  (-1 - NO NEIGHBOUR)'
    IF(ITYP.EQ.1) THEN
        READ(*,*) NEIGH(B,1:6)
        WRITE(1,*) NEIGH(B,1:6),  '   SNEIGH, NNEIGH, WNEIGH, ENEIGH, BNEIGH, TNEIGH '
    ELSE
        READ(2,*) NEIGH(B,1:6)
    END IF
       
end subroutine readData

!########################################################
subroutine cartesian
!########################################################

    use geo
    use ind
    implicit none
    
    ! Initialize all values with zero
    X=0.0d0
    Y=0.0d0
    Z=0.0d0
    XC=0.0d0
    YC=0.0d0
    ZC=0.0d0

    DX=(XXE-XXS)/dble(NICV)
    DY=(YYE-YYS)/dble(NJCV)
    DZ=(ZZE-ZZS)/dble(NKCV)
    VOL=DX*DY*DZ

    NI=NICV+2
    NJ=NJCV+2
    NK=NKCV+2
    NIJ=NI*NJ
    NIJK=NI*NJ*NK
    NIM=NI-1
    NJM=NJ-1
    NKM=NK-1
    N=NICV*NJCV*NKCV

    !Update global values
    NIA=NIA+NI
    NJA=NJA+NJ
    NKA=NKA+NK
    NIJA=NIJA+NIJ
    NIJKA=NIJKA+NIJK

    do I=1,NI
        LI(I)=(I-1)*NJ
    end do

    do K=1,NK
        LK(K)=(K-1)*NJ*NI
    end do

    do K=1,NKM
    do I=1,NIM
    do J=1,NJM
        IJK=LK(K)+LI(I)+J
        X(IJK)=XXS+dble(I-1)*DX
        Y(IJK)=YYS+dble(J-1)*DY 
        Z(IJK)=ZZS+dble(K-1)*DZ
    end do
    end do
    end do

end subroutine cartesian

!########################################################
subroutine gridExport
!########################################################

    use bc
    use ch
    use geo
    use ind
    implicit none

    character(10) :: OUT_CH

    DX=(XXE-XXS)/dble(NICV)
    DY=(YYE-YYS)/dble(NJCV)
    DZ=(ZZE-ZZS)/dble(NKCV)
    
    write(3,*)  NI,NJ,NK,NIJK,NBLOCK,NDIR
    write(3,*)  (NEIGH(B,I),I=1,6)
    print *, (NEIGH(B,I),I=1,6)
    write(3,*)  (LK(K),K=1,NK)
    write(3,*)  (LI(I),I=1,NI)

    write(3,*)  (X(I),I=1,NIJK)
    write(3,*)  (Y(I),I=1,NIJK)
    write(3,*)  (Z(I),I=1,NIJK)
    write(3,*)  (XC(I),I=1,NIJK)
    write(3,*)  (YC(I),I=1,NIJK)
    write(3,*)  (ZC(I),I=1,NIJK)

    write(3,*)  (IJKBL(I),I=1,NBLOCK)
    write(3,*)  (IJKPBL(I),I=1,NBLOCK)
    write(3,*)  (IJKBL1(I),I=1,NBLOCK)
    write(3,*)  (IJKBL2(I),I=1,NBLOCK)
    write(3,*)  (IJKBL3(I),I=1,NBLOCK)
    write(3,*)  (IJKBL4(I),I=1,NBLOCK)
    
    write(3,*)  (IJKDI(I),I=1,NDIR)
    write(3,*)  (IJKPDI(I),I=1,NDIR)
    write(3,*)  (IJKDI1(I),I=1,NDIR)
    write(3,*)  (IJKDI2(I),I=1,NDIR)
    write(3,*)  (IJKDI3(I),I=1,NDIR)
    write(3,*)  (IJKDI4(I),I=1,NDIR)

    write(3,*)  (FX(I), I=1,NIJK)
    write(3,*)  (FY(I), I=1,NIJK)
    write(3,*)  (FZ(I), I=1,NIJK)

    write(3,*)  DX,DY,DZ, VOL
    write(3,*)  (SRDDI(I),I=1,NDIR)

    !do I=1,NIJK
    !    write(3,*) XC(I), YC(I), ZC(I)
    !end do

    write(OUT_CH,'(I1)') B
    VTKFILE='grid_'//trim(OUT_CH)//'.vtk'

!...Create .vtk file
    print *, ' *** GENERATING .VTK *** '
    open (unit=4,FILE=VTKFILE)
    write(4,'(A)') '# vtk DataFile Version 3.0'
    write(4,'(A)') 'grid'
    write(4,'(A)') 'ASCII'
    write(4,'(A)') 'DATASET STRUCTURED_GRID'
    write(4,'(A I6 I6 I6)') 'DIMENSIONS', NIM,NJM,NKM
    write(4,'(A I10 A)') 'Points ', NIM*NJM*NKM, ' float'
    do K=1,NKM
    do I=1,NIM
    do J=1,NJM
        IJK=LK(K)+LI(I)+J
        write(4,'(E20.10,1X,E20.10,1X,E20.10)'), X(IJK), Y(IJK),Z(IJK)
    end do
    end do
    end do

end subroutine gridExport

!########################################################
subroutine setBc
!########################################################

    use bc
    use ind
    implicit none

    do L=1,6
        if(L.LE.2) then
            IK=0
            do K=2,NKM
            do I=2,NIM
                IK=IK+1
                ITB(L,IK)=BTYP(L)
            end do
            end do
        elseif(L.LE.4) then
            JK=0
            do K=2,NKM
            do J=2,NJM
                JK=JK+1
                JTB(L-2,JK)=BTYP(L)
            end do
            end do
        else
            IJ=0
            do I=2,NIM
            do J=2,NJM
                IJ=IJ+1
                KTB(L-4,IJ)=BTYP(L)
            end do
            end do
        end if
    end do

    call defbc(1,NDIR,IJKDI,IJKPDI,IJKDI1,IJKDI2,IJKDI3,IJKDI4)
    NDIRA=NDIRA+NDIR
    call defbc(2,NBLOCK,IJKBL,IJKPBL,IJKBL1,IJKBL2,IJKBL3,IJKBL4)
    NBLOCKA=NBLOCKA+NBLOCK

end subroutine setBc

!########################################################
subroutine defBc(LT,NBCF,IJKBB,IJKBP,IJK1,IJK2,IJK3,IJK4)
!########################################################

    use bc
    use ind
    implicit none
    integer, intent(in) :: LT
    integer, dimension(N), intent(inout) :: IJKBB,IJKBP,IJK1,IJK2,IJK3,IJK4
    integer, intent(out) :: NBCF
!
!.....COLLECT BOUNDARY CELL FACES OF TYPE 'LT' IN A LIST FOR EACH GRID
!
      NBCF=0
!
!.....SOUTH SIDE
!
    IK=0
    do K=2,NKM
    do I=2,NIM
        IK=IK+1
        IF(ITB(1,IK).EQ.LT) THEN
          NBCF=NBCF+1
          IJKBB(NBCF)=LK(K)+LI(I)+1
          IJKBP(NBCF)=IJKBB(NBCF)+1
          IJK4(NBCF)=IJKBB(NBCF)
          IJK3(NBCF)=IJK4(NBCF)-NJ
          IJK2(NBCF)=IJK3(NBCF)-NIJ
          IJK1(NBCF)=IJK4(NBCF)-NIJ
        end if
    end do
    end do
!
!.....NORTH SIDE
!
    IK=0
    do K=2,NKM
    do I=2,NIM
        IK=IK+1
        IF(ITB(2,IK).EQ.LT) THEN
            NBCF=NBCF+1
            IJKBB(NBCF)=LK(K)+LI(I)+NJ
            IJKBP(NBCF)=IJKBB(NBCF)-1
            IJK3(NBCF)=IJKBP(NBCF)
            IJK4(NBCF)=IJK3(NBCF)-NJ
            IJK2(NBCF)=IJK3(NBCF)-NIJ
            IJK1(NBCF)=IJK4(NBCF)-NIJ
        ENDIF
    end do
    end do

!
!.....WEST SIDE
!
    JK=0
    do K=2,NKM
    do J=2,NJM
        JK=JK+1
        IF(JTB(1,JK).EQ.LT) THEN
            NBCF=NBCF+1
            IJKBB(NBCF)=LK(K)+LI(1)+J
            IJKBP(NBCF)=IJKBB(NBCF)+NJ
            IJK3(NBCF)=IJKBB(NBCF)
            IJK4(NBCF)=IJK3(NBCF)-1
            IJK2(NBCF)=IJK3(NBCF)-NIJ
            IJK1(NBCF)=IJK4(NBCF)-NIJ
        ENDIF
    end do
    end do
!
!.....EAST SIDE
!
    JK=0
    do K=2,NKM
    do J=2,NJM
        JK=JK+1
        IF(JTB(2,JK).EQ.LT) THEN
            NBCF=NBCF+1
            IJKBB(NBCF)=LK(K)+LI(NI)+J
            IJKBP(NBCF)=IJKBB(NBCF)-NJ
            IJK4(NBCF)=IJKBP(NBCF)
            IJK3(NBCF)=IJK4(NBCF)-1
            IJK2(NBCF)=IJK3(NBCF)-NIJ
            IJK1(NBCF)=IJK4(NBCF)-NIJ
        end if
    end do
    end do
!
!......BOTTOM SIDE
!
    IJ=0
    do I=2,NIM
    do J=2,NJM
        IJ=IJ+1
        if (KTB(1,IJ).EQ.LT) then
            NBCF=NBCF+1
            IJKBB(NBCF)=LK(1)+LI(I)+J
            IJKBP(NBCF)=IJKBB(NBCF)+NIJ
            IJK3(NBCF)=IJKBB(NBCF)
            IJK4(NBCF)=IJK3(NBCF)-NJ
            IJK1(NBCF)=IJK4(NBCF)-1
            IJK2(NBCF)=IJK3(NBCF)-1
        end if
    end do
    end do
!
!......TOP SIDE
!
    IJ=0
    do I=2,NIM
    do J=2,NJM
        IJ=IJ+1
        if (KTB(2,I).EQ.LT) then
            NBCF=NBCF+1
            IJKBB(NBCF)=LK(NK)+LI(I)+J
            IJKBP(NBCF)=IJKBB(NBCF)-NIJ
            IJK4(NBCF)=IJKBP(NBCF)
            IJK3(NBCF)=IJK4(NBCF)-NJ
            IJK1(NBCF)=IJK4(NBCF)-1
            IJK2(NBCF)=IJK4(NBCF)-1
        end if
    end do
    end do

!
end subroutine defBc

!########################################################
subroutine calcG
!########################################################

    use bc
    use geo
    use ind
    implicit none
    real*8 ::   DLPE,DLEE,DLPN,DLNN,DLPT,DLTT,&
                XE,YE,ZE,XN,YN,ZN,XT,YT,ZT, &
                X1,X2,X4,Y1,Y2,Y4,Z1,Z2,Z4
    real*8,parameter :: SMALL=1.0E-20
    integer :: IJK1, IJK2, IJK3, IJK4

!.....CALCULATION OF NODE COORDINATES: CORNER (DUMMY) NODES
!
    IJK=LK(1)+LI(1)+1
    XC(IJK)=X(IJK)
    YC(IJK)=Y(IJK)
    ZC(IJK)=Z(IJK)
!
    IJK=LK(1)+LI(NIM)+1
    XC(IJK+NJ)=X(IJK)
    YC(IJK+NJ)=Y(IJK)
    ZC(IJK+NJ)=Z(IJK)
!
    IJK=LK(1)+LI(1)+NJM
    XC(IJK+1)=X(IJK)
    YC(IJK+1)=Y(IJK)
    ZC(IJK+1)=Z(IJK)
!
    IJK=LK(1)+LI(NIM)+NJM
    XC(IJK+NJ+1)=X(IJK)
    YC(IJK+NJ+1)=Y(IJK)
    ZC(IJK+NJ+1)=Z(IJK)
!
    IJK=LK(NKM)+LI(1)+1
    XC(IJK+NIJ)=X(IJK)
    YC(IJK+NIJ)=Y(IJK)
    ZC(IJK+NIJ)=Z(IJK)
!
    IJK=LK(NKM)+LI(NIM)+1
    XC(IJK+NIJ+NJ)=X(IJK)
    YC(IJK+NIJ+NJ)=Y(IJK)
    ZC(IJK+NIJ+NJ)=Z(IJK)
!
    IJK=LK(NKM)+LI(1)+NJM
    XC(IJK+NIJ+1)=X(IJK)
    YC(IJK+NIJ+1)=Y(IJK)
    ZC(IJK+NIJ+1)=Z(IJK)
!
    IJK=LK(NKM)+LI(NIM)+NJM
    XC(IJK+NIJ+NJ+1)=X(IJK)
    YC(IJK+NIJ+NJ+1)=Y(IJK)
    ZC(IJK+NIJ+NJ+1)=Z(IJK)

!
!.....CALCULATION OF NODE COORDINATES: BOUNDARY EDGE NODES
!               
    do I=2,NIM
        IJK=LK(1)+LI(I)+1
        XC(IJK)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK)=0.5*(Y(IJK)+Y(IJK-NJ))
        ZC(IJK)=0.5*(Z(IJK)+Z(IJK-NJ))
        IJK=LK(1)+LI(I)+NJM
        XC(IJK+1)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK+1)=0.5*(Y(IJK)+Y(IJK-NJ))
        ZC(IJK+1)=0.5*(Z(IJK)+Z(IJK-NJ))
        IJK=LK(NKM)+LI(I)+1
        XC(IJK+NIJ)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK+NIJ)=0.5*(Y(IJK)+Y(IJK-NJ))
        ZC(IJK+NIJ)=0.5*(Z(IJK)+Z(IJK-NJ))
        IJK=LK(NKM)+LI(I)+NJM
        XC(IJK+NIJ+1)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK+NIJ+1)=0.5*(Y(IJK)+Y(IJK-NJ))
        ZC(IJK+NIJ+1)=0.5*(Z(IJK)+Z(IJK-NJ))
    end do

    do J=2,NJM
        IJK=LK(1)+LI(NIM)+J
        XC(IJK+NJ)=0.5*(X(IJK)+X(IJK-1))
        YC(IJK+NJ)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK+NJ)=0.5*(Z(IJK)+Z(IJK-1))
        IJK=LK(1)+LI(1)+J
        XC(IJK)=0.5*(X(IJK)+X(IJK-1))
        YC(IJK)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK)=0.5*(Z(IJK)+Z(IJK-1))
        IJK=LK(NKM)+LI(NIM)+J
        XC(IJK+NIJ+NJ)=0.5*(X(IJK)+X(IJK-1))
        YC(IJK+NIJ+NJ)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK+NIJ+NJ)=0.5*(Z(IJK)+Z(IJK-1))
        IJK=LK(NKM)+LI(1)+J
        XC(IJK+NIJ)=0.5*(X(IJK)+X(IJK-1))
        YC(IJK+NIJ)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK+NIJ)=0.5*(Z(IJK)+Z(IJK-1))
    end do

    do K=2,NKM
        IJK=LK(K)+LI(NIM)+1
        XC(IJK+NJ)=0.5*(X(IJK)+X(IJK-NIJ))
        YC(IJK+NJ)=0.5*(Y(IJK)+Y(IJK-NIJ))
        ZC(IJK+NJ)=0.5*(Z(IJK)+Z(IJK-NIJ))
        IJK=LK(K)+LI(1)+1
        XC(IJK)=0.5*(X(IJK)+X(IJK+NIJ))
        YC(IJK)=0.5*(Y(IJK)+Y(IJK+NIJ))
        ZC(IJK)=0.5*(Z(IJK)+Z(IJK-NIJ))
        IJK=LK(K)+LI(NIM)+NJM
        XC(IJK+NJ+1)=0.5*(X(IJK)+X(IJK-NIJ))
        YC(IJK+NJ+1)=0.5*(Y(IJK)+Y(IJK-NIJ))
        ZC(IJK+NJ+1)=0.5*(Z(IJK)+Z(IJK-NIJ))
        IJK=LK(K)+LI(1)+NJM
        XC(IJK+1)=0.5*(X(IJK)+X(IJK-NIJ))
        YC(IJK+1)=0.5*(Y(IJK)+Y(IJK-NIJ))
        ZC(IJK+1)=0.5*(Z(IJK)+Z(IJK-NIJ))
    end do
!
!.....CALCULATION OF NODE COORDINATES: BOUNDARY FACE NODES
!
    do I=2,NIM
    do J=2,NJM
        IJK=LK(1)+LI(I)+J
        XC(IJK)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK)=0.5*(Z(IJK)+Z(IJK-NJ))
        IJK=LK(NKM)+LI(I)+J
        XC(IJK+NIJ)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK+NIJ)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK+NIJ)=0.5*(Z(IJK)+Z(IJK-NJ))
    end do
    end do
    
    do J=2,NJM
    do K=2,NKM
        IJK=LK(K)+LI(NIM)+J
        XC(IJK+NJ)=0.5*(X(IJK)+X(IJK-1))
        YC(IJK+NJ)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK+NJ)=0.5*(Z(IJK)+Z(IJK-NIJ))
        IJK=LK(K)+LI(1)+J
        XC(IJK)=0.5*(X(IJK)+X(IJK-1))
        YC(IJK)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK)=0.5*(Z(IJK)+Z(IJK-NIJ))
    end do
    end do
            
    do K=2,NKM
    do I=2,NIM
        IJK=LK(K)+LI(I)+1
        XC(IJK)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK)=0.5*(Y(IJK)+Y(IJK-NJ))
        ZC(IJK)=0.5*(Z(IJK)+Z(IJK-NIJ))
        IJK=LK(K)+LI(I)+NJM
        XC(IJK+1)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK+1)=0.5*(Y(IJK)+Y(IJK-NJ))
        ZC(IJK+1)=0.5*(Z(IJK)+Z(IJK-NIJ))
    end do
    end do
!
!.....CALCULATION OF NODE COORDINATES: CELL CENTERS
!
    do K=2,NKM
    do I=2,NIM
    do J=2,NJM
          IJK=LK(K)+LI(I)+J
          XC(IJK)=0.5d0*(X(IJK)+X(IJK-NJ))
          YC(IJK)=0.5d0*(Y(IJK)+Y(IJK-1))
          ZC(IJK)=0.5d0*(Z(IJK)+Z(IJK-NIJ))
    end do
    end do
    end do
!
!......CALCULATION OF INTERPOLATION FACTORS
!
    do K=2,NKM
    do I=2,NIM
    do J=2,NJM
        IJK=LK(K)+LI(I)+J
!
!.....INTERPOLATION IN I-DIRECTION: FX = Pe/PE
!
        XE=0.5d0*(X(IJK)+X(IJK-1))
        YE=0.5d0*(Y(IJK)+Y(IJK-1))
        ZE=0.5d0*(Z(IJK)+Z(IJK-NIJ))
        DLPE=SQRT((XE-XC(IJK))**2+(YE-YC(IJK))**2+(ZE-ZC(IJK))**2)
        DLEE=SQRT((XC(IJK+NJ)-XE)**2+(YC(IJK+NJ)-YE)**2+(ZC(IJK+NJ)-ZE)**2)
        FX(IJK)=DLPE/(DLPE+DLEE+SMALL)
!
!.....INTERPOLATION IN J-DIRECTION: FY = Pn/PN
!
        XN=0.5d0*(X(IJK)+X(IJK-NJ))
        YN=0.5d0*(Y(IJK)+Y(IJK-NJ))
        ZN=0.5d0*(Z(IJK)+Z(IJK-NIJ))
        DLPN=SQRT((XN-XC(IJK))**2+(YN-YC(IJK))**2+(ZN-ZC(IJK))**2)
        DLNN=SQRT((XC(IJK+1)-XN)**2+(YC(IJK+1)-YN)**2+(ZC(IJK+1)-ZN)**2)
        FY(IJK)=DLPN/(DLPN+DLNN+SMALL)
!
!.....INTERPOLATION IN K-DIRECTION: FZ = Pt/PT
!
        XT=0.5d0*(X(IJK)+X(IJK-NJ))
        YT=0.5d0*(Y(IJK)+Y(IJK-1))
        ZT=0.5d0*(Z(IJK)+Z(IJK-1))
        DLPT=SQRT((XT-XC(IJK))**2+(YT-YC(IJK))**2+(ZT-ZC(IJK))**2)
        DLTT=SQRT((XC(IJK+NIJ)-XT)**2+(YC(IJK+NIJ)-YT)**2+(ZC(IJK+NIJ)-ZT)**2)
        FZ(IJK)=DLPT/(DLPT+DLTT+SMALL)
    end do
    end do
    end do

!
!....Normal distance from cell face center to cell center
!
    do IDI=1,NDIR
        IJKB=IJKDI(IDI)
        IJKP=IJKPDI(IDI)
        !IJK1=IJKDI1(IDI)
        IJK2=IJKDI2(IDI)
        IJK3=IJKDI3(IDI)
        IJK4=IJKDI4(IDI)
        !
        call normalArea(IJKP,IJKB,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
        !
        SRDDI(IDI)=AR/(DN+SMALL)
    end do

end subroutine calcG

!########################################################
subroutine writeParamMod
!########################################################

    use bc
    use geo
    use ind
    implicit none

!...Create solver file
    !OPEN(UNIT=9,FILE='../../../pet_src/paramMod.F90')
    OPEN(UNIT=9,FILE='paramMod.F90')
    REWIND 9
    write(9,'(A16)') 'module param'
    write(9,'(A18)') '  implicit none'
    write(9,'(A22 A4 I6 A1 A5 I9 A5 I9 A7 I9 A8 I9 A10 I9 A9 I9 A6 I9)') 'integer, parameter :: ', 'NXA=', NIA,'&' 
    write(9,'(A5 I6 A1)') ',NYA=', NJA, '&'
    write(9,'(A5 I6 A1)') ',NZA=', NKA, '&'
    write(9,'(A7 I6 A1)') ',NXYZA=', NIJKA, '&'
    write(9,'(A8 I6 A1)') ',NDIRAL=', NDIRA,'&'
    write(9,'(A10 I6 A1)') ',NBLOCKAL=', NBLOCKA, '&'
    write(9,'(A9 I6 A1)')   ',NBLOCKS=',NB,'&'
    write(9,'(A6 I1)') ',PREC=',PREC
    write(9,'(A)') 'end module param'

end subroutine writeParamMod
