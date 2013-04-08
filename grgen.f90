!########################################################
program grgen
!########################################################

    use indexModule
    use charModule
    implicit none

    print *, ' TOTAL NUMBER OF BLOCKS: '
    read(*,*) NB
    print *, ' "NAME" OF INPUT FILE (NAME_BLOCK.inp, * - KEYBOARD)'
    read *, FILIN

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

    use boundaryModule
    use charModule
    use geoModule
    use indexModule
    implicit none

    integer :: ITYP
    
    BLOCKUNIT=OFFSET+B

    write(BLOCK_CH,'(I1)') B
    !PRINT *, ' INPUT FILE NAME (* - KEYBOARD):  '
    !READ(*,1) FILIN
    IF(FILIN.NE.'*') THEN
        BLOCKFILE=trim(FILIN)//'_'//trim(BLOCK_CH)//'.inp'
        OPEN (UNIT=BLOCKUNIT,FILE=BLOCKFILE)
        REWIND BLOCKUNIT
        ITYP=0
    ELSE
        BLOCKFILE='grid_'//trim(BLOCK_CH)//'.inp'
        OPEN (UNIT=BLOCKUNIT,FILE=BLOCKFILE)
        REWIND BLOCKUNIT
        ITYP=1
    ENDIF

    PRINT *, ' OUTPUT FILE NAME:  '
    IF(ITYP.EQ.1) THEN
        READ(*,1) FILOUT
        WRITE(BLOCKUNIT,1) FILOUT
    ELSE
        READ(BLOCKUNIT,1) FILOUT
        print *, FILOUT
    ENDIF
    1 FORMAT(A12)

    PRINT *, ' ENTER> XSTART, XEND, NUMBER OF CVS:  '
    IF(ITYP.EQ.1) THEN
        READ(*,*) XXS,XXE,NICV
        WRITE(BLOCKUNIT,*) XXS,XXE,NICV,'   XS,XE,NICV '
    ELSE
        READ(BLOCKUNIT,*) XXS,XXE,NICV
        print *, XXS,XXE,NICV
    END IF

    PRINT *, ' ENTER> YSTART, YEND, NUMBER OF CVS:  '
    IF(ITYP.EQ.1) THEN
        READ(*,*) YYS,YYE,NJCV
        WRITE(BLOCKUNIT,*) YYS,YYE,NJCV,'   YS,YE,NJCV '
    ELSE
        READ(BLOCKUNIT,*) YYS,YYE,NJCV
        print *, YYS,YYE,NJCV
    END IF

    PRINT *, ' ENTER> ZSTART, ZEND, NUMBER OF CVS:  '
    IF(ITYP.EQ.1) THEN
        READ(*,*) ZZS,ZZE,NKCV
        WRITE(BLOCKUNIT,*) ZZS,ZZE,NKCV,'   ZS,ZE,NKCV '
    ELSE
        READ(BLOCKUNIT,*) ZZS,ZZE,NKCV
        print *, ZZS,ZZE,NKCV
    END IF

    PRINT *, ' ENTER> BOUNDARY TYPE S N W E B T:'
    print *, '(1 - DIRICHLET, 2 - NEUMANN ZERO GRADIENT, 3 WALL (DIRICHLET), 4 BLOCK)'
    IF(ITYP.EQ.1) THEN
        READ(*,*) BTYP(1:6)
        WRITE(BLOCKUNIT,*) BTYP(1:6),  '   SBTYP, NBTYP, WBTYP, EBTYP,  BBTYP ,TBTYP '
    ELSE
        READ(BLOCKUNIT,*) BTYP(1:6)
        print *, BTYP(1:6)
    END IF
    
    !if(NB.gt.1) then
        PRINT *, ' ENTER> NEIGHBOUR INDEX S N W E T B:  (-1 - NO NEIGHBOUR)'
        IF(ITYP.EQ.1) THEN
            READ(*,*) NEIGH(B,1:6)
            WRITE(BLOCKUNIT,*) NEIGH(B,1:6),  '   SNEIGH, NNEIGH, WNEIGH, ENEIGH, BNEIGH, TNEIGH '
        ELSE
            READ(BLOCKUNIT,*) NEIGH(B,1:6)
            print *, NEIGH(B,1:6)
        END IF
    !else
    !    NEIGH(1,:)=-1
    !end if

    close(unit=BLOCKUNIT)
    OPEN (UNIT=BLOCKUNIT,FILE=FILOUT)
    REWIND BLOCKUNIT

       
end subroutine readData

!========================================================
!>   cartesian orthogonal 3d grid
!########################################################
subroutine cartesian
!########################################################

    use geoModule
    use indexModule
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

    !do I=1,NI
    !    !LI(I)=(I-1)*NJ
    !end do

    !do K=1,NK
    !    LK(K)=(K-1)*NJ*NI
    !end do

    do K=1,NKM
    do I=1,NIM
    do J=1,NJM
        !IJK=LK(K)+LI(I)+J
        IJK=(K-1)*NI*NJ+(I-1)*NJ+J
        X(IJK)=XXS+dble(I-1)*DX
        Y(IJK)=YYS+dble(J-1)*DY 
        Z(IJK)=ZZS+dble(K-1)*DZ
        !print *, X(IJK), Y(IJK), Z(IJK)
    end do
    end do
    end do

end subroutine cartesian

!########################################################
subroutine gridExport
!########################################################

    use boundaryModule
    use charModule
    use geoModule
    use indexModule
    implicit none

    DX=(XXE-XXS)/dble(NICV)
    DY=(YYE-YYS)/dble(NJCV)
    DZ=(ZZE-ZZS)/dble(NKCV)

    
    write(BLOCKUNIT,*)  NI,NJ,NK,NIJK,NDIR,NNEU,NWAL,NBLO
    write(BLOCKUNIT,*)  (NEIGH(B,I),I=1,6)
    !print *, (NEIGH(B,I),I=1,6)
    !write(3,*)  (LK(K),K=1,NK)
    !write(3,*)  (LI(I),I=1,NI)

    write(BLOCKUNIT,*)  (X(I),I=1,NIJK)
    write(BLOCKUNIT,*)  (Y(I),I=1,NIJK)
    write(BLOCKUNIT,*)  (Z(I),I=1,NIJK)
    write(BLOCKUNIT,*)  (XC(I),I=1,NIJK)
    write(BLOCKUNIT,*)  (YC(I),I=1,NIJK)
    write(BLOCKUNIT,*)  (ZC(I),I=1,NIJK)

    write(BLOCKUNIT,*)  (IJKBDIR(I),I=1,NDIR)
    write(BLOCKUNIT,*)  (IJKPDIR(I),I=1,NDIR)
    write(BLOCKUNIT,*)  (IJKDIR1(I),I=1,NDIR)
    write(BLOCKUNIT,*)  (IJKDIR2(I),I=1,NDIR)
    write(BLOCKUNIT,*)  (IJKDIR3(I),I=1,NDIR)
    write(BLOCKUNIT,*)  (IJKDIR4(I),I=1,NDIR)
    
    write(BLOCKUNIT,*)  (IJKBNEU(I),I=1,NNEU)
    write(BLOCKUNIT,*)  (IJKPNEU(I),I=1,NNEU)
    write(BLOCKUNIT,*)  (IJKNEU1(I),I=1,NNEU)
    write(BLOCKUNIT,*)  (IJKNEU2(I),I=1,NNEU)
    write(BLOCKUNIT,*)  (IJKNEU3(I),I=1,NNEU)
    write(BLOCKUNIT,*)  (IJKNEU4(I),I=1,NNEU)
    
    write(BLOCKUNIT,*)  (IJKBWAL(I),I=1,NWAL)
    write(BLOCKUNIT,*)  (IJKPWAL(I),I=1,NWAL)
    write(BLOCKUNIT,*)  (IJKWAL1(I),I=1,NWAL)
    write(BLOCKUNIT,*)  (IJKWAL2(I),I=1,NWAL)
    write(BLOCKUNIT,*)  (IJKWAL3(I),I=1,NWAL)
    write(BLOCKUNIT,*)  (IJKWAL4(I),I=1,NWAL)

    write(BLOCKUNIT,*)  (IJKBBLO(I),I=1,NBLO)
    write(BLOCKUNIT,*)  (IJKPBLO(I),I=1,NBLO)
    write(BLOCKUNIT,*)  (IJKBLO1(I),I=1,NBLO)
    write(BLOCKUNIT,*)  (IJKBLO2(I),I=1,NBLO)
    write(BLOCKUNIT,*)  (IJKBLO3(I),I=1,NBLO)
    write(BLOCKUNIT,*)  (IJKBLO4(I),I=1,NBLO)
    
    write(BLOCKUNIT,*)  (FX(I), I=1,NIJK)
    write(BLOCKUNIT,*)  (FY(I), I=1,NIJK)
    write(BLOCKUNIT,*)  (FZ(I), I=1,NIJK)

    write(BLOCKUNIT,*)  DX,DY,DZ,VOL
    !write(BLOCKUNIT,*)  (SRDDIR(I),I=1,NDIR)
    !write(BLOCKUNIT,*)  (SRDNEU(I),I=1,NNEU)
    write(BLOCKUNIT,*)  (SRDWAL(I),I=1,NWAL)

    close(unit=BLOCKUNIT)

!
!.....Create .vtk file
!
    write(BLOCK_CH,'(I1)') B
    VTKFILE='grid_'//trim(BLOCK_CH)//'.vtk'
    open (unit=BLOCKUNIT,FILE=VTKFILE)
    rewind BLOCKUNIT

    print *, ' *** GENERATING .VTK *** '
    write(BLOCKUNIT,'(A)') '# vtk DataFile Version 3.0'
    write(BLOCKUNIT,'(A)') 'grid'
    write(BLOCKUNIT,'(A)') 'ASCII'
    write(BLOCKUNIT,'(A)') 'DATASET STRUCTURED_GRID'
    write(BLOCKUNIT,'(A I6 I6 I6)') 'DIMENSIONS', NIM,NJM,NKM
    write(BLOCKUNIT,'(A I10 A)') 'Points ', NIM*NJM*NKM, ' float'
    do K=1,NKM
    do J=1,NJM
    do I=1,NIM
        IJK=(K-1)*NI*NJ+(I-1)*NJ+J
        write(BLOCKUNIT,'(E20.10,1X,E20.10,1X,E20.10)'), X(IJK), Y(IJK),Z(IJK)
    end do
    end do
    end do

    close(unit=BLOCKUNIT)

end subroutine gridExport

!=======================================================
!>  writes the boundary type of each boundary cell into
!>  an array
!########################################################
subroutine setBc
!########################################################

    use boundaryModule
    use indexModule
    implicit none

    do P=1,6
        if(P.LE.2) then
            IK=0
            do K=2,NKM
            do I=2,NIM
                IK=IK+1
                ITB(P,IK)=BTYP(P)
            end do
            end do
        elseif(P.LE.4) then
            JK=0
            do K=2,NKM
            do J=2,NJM
                JK=JK+1
                JTB(P-2,JK)=BTYP(P)
            end do
            end do
        else
            IJ=0
            do I=2,NIM
            do J=2,NJM
                IJ=IJ+1
                KTB(P-4,IJ)=BTYP(P)
            end do
            end do
        end if
    end do

    call defbc(1,NDIR,IJKBDIR,IJKPDIR,IJKDIR1,IJKDIR2,IJKDIR3,IJKDIR4)
    NDIRA=NDIRA+NDIR
    call defbc(2,NNEU,IJKBNEU,IJKPNEU,IJKNEU1,IJKNEU2,IJKNEU3,IJKNEU4)
    NNEUA=NNEUA+NNEU
    call defbc(3,NWAL,IJKBWAL,IJKPWAL,IJKWAL1,IJKWAL2,IJKWAL3,IJKWAL4)
    NWALA=NWALA+NWAL
    call defbc(4,NBLO,IJKBBLO,IJKPBLO,IJKBLO1,IJKBLO2,IJKBLO3,IJKBLO4)
    NBLOA=NBLOA+NBLO

end subroutine setBc

!=========================================================
!>  collects the following information about the boundary
!>    cells of the specified boundary type: 
!>    - IJKBB: index of boundary cell face center (neighbour)
!>    - IJKBP: index of boundary cell center
!>    - IJK1...4: index of boundary the edges of the respective
!>    boundary cell face. They have to be ordered in a 
!>    clockwise manner if the cell face is viewed from 
!>    outside the CV.
!########################################################
subroutine defBc(LT,NBCF,IJKBB,IJKBP,IJK1,IJK2,IJK3,IJK4)
!########################################################

    use boundaryModule
    use indexModule
    implicit none
    integer, intent(in) :: LT
    integer, dimension(N), intent(inout) :: IJKBB,IJKBP,IJK1,IJK2,IJK3,IJK4
    integer, intent(out) :: NBCF
!
!.....COLLECT BOUNDARY CELL FACES OF TYPE 'LT' IN A LIST FOR EACH GRID
!
      NBCF=0
!
!.....SNEUH SIDE
!
    IK=0
    do K=2,NKM
    do I=2,NIM
        IK=IK+1
        IF(ITB(1,IK).EQ.LT) THEN
          NBCF=NBCF+1
          !IJKBB(NBCF)=LK(K)+LI(I)+1
          IJKBB(NBCF)=(K-1)*NI*NJ+(I-1)*NJ+1
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
            !IJKBB(NBCF)=LK(K)+LI(I)+NJ
            IJKBB(NBCF)=(K-1)*NI*NJ+(I-1)*NJ+NJ
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
            !IJKBB(NBCF)=LK(K)+LI(1)+J
            IJKBB(NBCF)=(K-1)*NI*NJ+J
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
            !IJKBB(NBCF)=LK(K)+LI(NI)+J
            IJKBB(NBCF)=(K-1)*NI*NJ+(NI-1)*NJ+J
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
            !IJKBB(NBCF)=LK(1)+LI(I)+J
            IJKBB(NBCF)=(I-1)*NJ+J
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
            !IJKBB(NBCF)=LK(NK)+LI(I)+J
            IJKBB(NBCF)=(NK-1)*NI*NJ+(I-1)*NJ+J
            IJKBP(NBCF)=IJKBB(NBCF)-NIJ
            IJK4(NBCF)=IJKBP(NBCF)
            IJK3(NBCF)=IJK4(NBCF)-NJ
            IJK1(NBCF)=IJK4(NBCF)-1
            IJK2(NBCF)=IJK3(NBCF)-1
        end if
    end do
    end do

!
end subroutine defBc

!========================================================
!>  calculates all missing corner, edge and face
!>  coordinates and the interpolation factors FX,FY,FZ
!########################################################
subroutine calcG
!########################################################

    use boundaryModule
    use geoModule
    use indexModule
    implicit none
    real*8 ::   DLPE,DLEE,DLPN,DLNN,DLPT,DLTT,&
                XE,YE,ZE,XN,YN,ZN,XT,YT,ZT, &
                X1,X2,X4,Y1,Y2,Y4,Z1,Z2,Z4
    real*8,parameter :: SMALL=1.0E-20
    integer :: IJK1, IJK2, IJK3, IJK4
!
!.....CALCULATION OF NODE COORDINATES: CORNER (DUMMY) NODES
!
    !IJK=LK(1)+LI(1)+1
    IJK=1
    XC(IJK)=X(IJK)
    YC(IJK)=Y(IJK)
    ZC(IJK)=Z(IJK)
 
    !IJK=LK(1)+LI(NIM)+1
    IJK=(NIM-1)*NJ+1
    XC(IJK+NJ)=X(IJK)
    YC(IJK+NJ)=Y(IJK)
    ZC(IJK+NJ)=Z(IJK)
 
    !IJK=LK(1)+LI(1)+NJM
    IJK=NJM
    XC(IJK+1)=X(IJK)
    YC(IJK+1)=Y(IJK)
    ZC(IJK+1)=Z(IJK)
 
    !IJK=LK(1)+LI(NIM)+NJM
    IJK=(NIM-1)*NJ+NJM
    XC(IJK+NJ+1)=X(IJK)
    YC(IJK+NJ+1)=Y(IJK)
    ZC(IJK+NJ+1)=Z(IJK)
 
    !IJK=LK(NKM)+LI(1)+1
    IJK=(NKM-1)*NI*NJ+1
    XC(IJK+NIJ)=X(IJK)
    YC(IJK+NIJ)=Y(IJK)
    ZC(IJK+NIJ)=Z(IJK)
 
    !IJK=LK(NKM)+LI(NIM)+1
    IJK=(NKM-1)*NI*NJ+(NIM-1)*NJ+1
    XC(IJK+NIJ+NJ)=X(IJK)
    YC(IJK+NIJ+NJ)=Y(IJK)
    ZC(IJK+NIJ+NJ)=Z(IJK)
 
    !IJK=LK(NKM)+LI(1)+NJM
    IJK=(NKM-1)*NI*NJ+NJM
    XC(IJK+NIJ+1)=X(IJK)
    YC(IJK+NIJ+1)=Y(IJK)
    ZC(IJK+NIJ+1)=Z(IJK)
 
    !IJK=LK(NKM)+LI(NIM)+NJM
    IJK=(NKM-1)*NI*NJ+(NIM-1)*NJ+NJM
    XC(IJK+NIJ+NJ+1)=X(IJK)
    YC(IJK+NIJ+NJ+1)=Y(IJK)
    ZC(IJK+NIJ+NJ+1)=Z(IJK)
 
!
!.....CALCULATION OF NODE COORDINATES: BOUNDARY EDGE NODES
!               
    do I=2,NIM
        !IJK=LK(1)+LI(I)+1
        IJK=(I-1)*NJ+1
        XC(IJK)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK)=0.5*(Y(IJK)+Y(IJK-NJ))
        ZC(IJK)=0.5*(Z(IJK)+Z(IJK-NJ))
        !IJK=LK(1)+LI(I)+NJM
        IJK=(I-1)*NJ+NJM
        XC(IJK+1)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK+1)=0.5*(Y(IJK)+Y(IJK-NJ))
        ZC(IJK+1)=0.5*(Z(IJK)+Z(IJK-NJ))
        !IJK=LK(NKM)+LI(I)+1
        IJK=(NKM-1)*NI*NJ+(I-1)*NJ+1
        XC(IJK+NIJ)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK+NIJ)=0.5*(Y(IJK)+Y(IJK-NJ))
        ZC(IJK+NIJ)=0.5*(Z(IJK)+Z(IJK-NJ))
        !IJK=LK(NKM)+LI(I)+NJM
        IJK=(NKM-1)*NI*NJ+(I-1)*NJ+NJM
        XC(IJK+NIJ+1)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK+NIJ+1)=0.5*(Y(IJK)+Y(IJK-NJ))
        ZC(IJK+NIJ+1)=0.5*(Z(IJK)+Z(IJK-NJ))
    end do

    do J=2,NJM
        !IJK=LK(1)+LI(NIM)+J
        IJK=(NIM-1)*NJ+J
        XC(IJK+NJ)=0.5*(X(IJK)+X(IJK-1))
        YC(IJK+NJ)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK+NJ)=0.5*(Z(IJK)+Z(IJK-1))
        !IJK=LK(1)+LI(1)+J
        IJK=J
        XC(IJK)=0.5*(X(IJK)+X(IJK-1))
        YC(IJK)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK)=0.5*(Z(IJK)+Z(IJK-1))
        !IJK=LK(NKM)+LI(NIM)+J
        IJK=(NKM-1)*NI*NJ+(NIM-1)*NJ+J
        XC(IJK+NIJ+NJ)=0.5*(X(IJK)+X(IJK-1))
        YC(IJK+NIJ+NJ)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK+NIJ+NJ)=0.5*(Z(IJK)+Z(IJK-1))
        !IJK=LK(NKM)+LI(1)+J
        IJK=(NKM-1)*NI*NJ+J
        XC(IJK+NIJ)=0.5*(X(IJK)+X(IJK-1))
        YC(IJK+NIJ)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK+NIJ)=0.5*(Z(IJK)+Z(IJK-1))
    end do

    do K=2,NKM
        !IJK=LK(K)+LI(NIM)+1
        IJK=(K-1)*NI*NJ+(NIM-1)*NJ+1
        XC(IJK+NJ)=0.5*(X(IJK)+X(IJK-NIJ))
        YC(IJK+NJ)=0.5*(Y(IJK)+Y(IJK-NIJ))
        ZC(IJK+NJ)=0.5*(Z(IJK)+Z(IJK-NIJ))
        !IJK=LK(K)+LI(1)+1
        IJK=(K-1)*NI*NJ+1
        XC(IJK)=0.5*(X(IJK)+X(IJK+NIJ))
        YC(IJK)=0.5*(Y(IJK)+Y(IJK+NIJ))
        ZC(IJK)=0.5*(Z(IJK)+Z(IJK-NIJ))
        !IJK=LK(K)+LI(NIM)+NJM
        IJK=(K-1)*NI*NJ+(NIM-1)*NJ+NJM
        XC(IJK+NJ+1)=0.5*(X(IJK)+X(IJK-NIJ))
        YC(IJK+NJ+1)=0.5*(Y(IJK)+Y(IJK-NIJ))
        ZC(IJK+NJ+1)=0.5*(Z(IJK)+Z(IJK-NIJ))
        !IJK=LK(K)+LI(1)+NJM
        IJK=(K-1)*NI*NJ+NJM
        XC(IJK+1)=0.5*(X(IJK)+X(IJK-NIJ))
        YC(IJK+1)=0.5*(Y(IJK)+Y(IJK-NIJ))
        ZC(IJK+1)=0.5*(Z(IJK)+Z(IJK-NIJ))
    end do
!
!.....CALCULATION OF NODE COORDINATES: BOUNDARY FACE NODES
!
    do I=2,NIM
    do J=2,NJM
        !IJK=LK(1)+LI(I)+J
        IJK=(I-1)*NJ+J
        XC(IJK)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK)=0.5*(Z(IJK)+Z(IJK-NJ))
        !IJK=LK(NKM)+LI(I)+J
        IJK=(NKM-1)*NI*NJ+(I-1)*NJ+J
        XC(IJK+NIJ)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK+NIJ)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK+NIJ)=0.5*(Z(IJK)+Z(IJK-NJ))
    end do
    end do
    
    do J=2,NJM
    do K=2,NKM
        !IJK=LK(K)+LI(NIM)+J
        IJK=(K-1)*NI*NJ+(NIM-1)*NJ+J
        XC(IJK+NJ)=0.5*(X(IJK)+X(IJK-1))
        YC(IJK+NJ)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK+NJ)=0.5*(Z(IJK)+Z(IJK-NIJ))
        !IJK=LK(K)+LI(1)+J
        IJK=(K-1)*NI*NJ+J
        XC(IJK)=0.5*(X(IJK)+X(IJK-1))
        YC(IJK)=0.5*(Y(IJK)+Y(IJK-1))
        ZC(IJK)=0.5*(Z(IJK)+Z(IJK-NIJ))
    end do
    end do
            
    do K=2,NKM
    do I=2,NIM
        !IJK=LK(K)+LI(I)+1
        IJK=(K-1)*NI*NJ+(I-1)*NJ+1
        XC(IJK)=0.5*(X(IJK)+X(IJK-NJ))
        YC(IJK)=0.5*(Y(IJK)+Y(IJK-NJ))
        ZC(IJK)=0.5*(Z(IJK)+Z(IJK-NIJ))
        !IJK=LK(K)+LI(I)+NJM
        IJK=(K-1)*NI*NJ+(I-1)*NJ+NJM
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
          !IJK=LK(K)+LI(I)+J
          IJK=(K-1)*NI*NJ+(I-1)*NJ+J
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
        !IJK=LK(K)+LI(I)+J
        IJK=(K-1)*NI*NJ+(I-1)*NJ+J
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
    !do IJKDIR=1,NDIR
    !    IJKB=IJKBDIR(IJKDIR)
    !    IJKP=IJKPDIR(IJKDIR)
    !    !IJK1=IJKDIR1(IJKDIR)
    !    IJK2=IJKDIR2(IJKDIR)
    !    IJK3=IJKDIR3(IJKDIR)
    !    IJK4=IJKDIR4(IJKDIR)
    !    !
    !    call normalArea(IJKP,IJKB,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
    !    !
    !    SRDDIR(IJKDIR)=AR/(DN+SMALL)
    !end do

    !do IJKNEU=1,NNEU
    !    IJKB=IJKBNEU(IJKNEU)
    !    IJKP=IJKPNEU(IJKNEU)
    !    !IJK1=IJKNEU1(IJKNEU)
    !    IJK2=IJKNEU2(IJKNEU)
    !    IJK3=IJKNEU3(IJKNEU)
    !    IJK4=IJKNEU4(IJKNEU)
    !    !
    !    call normalArea(IJKP,IJKB,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
    !    !
    !    SRDNEU(IJKNEU)=AR/(DN+SMALL)
    !end do

    do IJKWAL=1,NWAL
        IJKB=IJKBWAL(IJKWAL)
        IJKP=IJKPWAL(IJKWAL)
        !IJK1=IJKWAL1(IJKWAL)
        IJK2=IJKWAL2(IJKWAL)
        IJK3=IJKWAL3(IJKWAL)
        IJK4=IJKWAL4(IJKWAL)
        !
        call normalArea(IJKP,IJKB,IJK2,IJK3,IJK4,AR,DN,XPN,YPN,ZPN,NX,NY,NZ)
        !
        SRDWAL(IJKWAL)=AR/(DN+SMALL)
    end do

end subroutine calcG

!########################################################
subroutine writeParamMod
!########################################################

    use boundaryModule
    use geoModule
    use indexModule
    implicit none

!...Create solver file
    OPEN(UNIT=9,FILE='parameterModule.f90')
    REWIND 9
    
    write(9,'(A22)') 'module parameterModule'
    write(9,'(4X, A13)') 'implicit none'
    write(9,'(4X, A22, A4, I6)') 'integer, parameter :: ', 'NXA=', NIA
    write(9,'(4X, A22, A4, I6)') 'integer, parameter :: ', 'NYA=', NJA
    write(9,'(4X, A22, A4, I6)') 'integer, parameter :: ', 'NZA=', NKA
    write(9,'(4X, A22, A6, I9)') 'integer, parameter :: ', 'NXYZA=', NIJKA
    write(9,'(4X, A22, A7, I6)') 'integer, parameter :: ', 'NDIRAL=', NDIRA
    write(9,'(4X, A22, A7, I6)') 'integer, parameter :: ', 'NNEUAL=', NNEUA
    write(9,'(4X, A22, A7, I6)') 'integer, parameter :: ', 'NWALAL=', NWALA
    write(9,'(4X, A22, A7, I6)') 'integer, parameter :: ', 'NBLOAL=', NBLOA
    write(9,'(4X, A22, A8, I6)') 'integer, parameter :: ', 'NBLOCKS=',NB
    write(9,'(4X, A22, A5, I1)') 'integer, parameter :: ', 'PREC=',PREC
    write(9,'(4X, A22, A4, I6)') 'integer, parameter :: ', 'NAL=',(NIA-2)*(NJA-2)*(NKA-2)
    write(9,'(4X, A22, A8, I6)') 'integer, parameter :: ', 'NFACEAL=',100000
    write(9,'(A)') 'end module parameterModule'

end subroutine writeParamMod
