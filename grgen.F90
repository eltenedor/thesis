!########################################################
program grgen
!########################################################

    implicit none

    call readData
    call cartesian
    call setBc
    call calcG
    call gridExport

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

    PRINT *, ' ENTER> BOUNDARY TYPE N S W E:  '
    IF(ITYP.EQ.1) THEN
        READ(*,*) BTYP(1:4)
        WRITE(1,*) BTYP(1:4),  '   SBTYP, NBTYP, WBTYP, EBTYP '
    ELSE
        READ(2,*) BTYP(1:4)
    END IF
       
end subroutine readData

!########################################################
subroutine cartesian
!########################################################

    use geo
    use ind
    implicit none

    DX=(XXE-XXS)/NICV
    DY=(YYE-YYS)/NJCV
    VOL=DX*DY

    NI=NICV+2
    NJ=NJCV+2
    NIJ=NI*NJ
    NIM=NI-1
    NJM=NJ-1
    IJST=1
    N=NICV*NJCV

    do i=1,NI
        LI(I)=(I-1)*NJ
    end do

    call createMapping
    
    do I=1,NIM
       do J=1,NJM
            IJ=LI(I)+J
            X(IJ)=XXS+(I-1)*DX
            Y(IJ)=YYS+(J-1)*DY 
        end do
    end do

end subroutine cartesian

!########################################################
subroutine createMapping
!########################################################

    use ind
    implicit none

    IJP=-1
    do I=2,NIM
        do J=2,NJM
            IJP=IJP+1
            IJ=LI(I)+J
            CTD(IJP)=IJ
            DTC(IJ)=IJP
        end do
    end do

end subroutine createMapping

!########################################################
subroutine gridexport
!########################################################

    use bc
    use geo
    use ind
    implicit none

    DX=(XXE-XXS)/NICV
    DY=(YYE-YYS)/NJCV
!...Create solver file
    !OPEN(UNIT=9,FILE='../../../pet_src/paramMod.F90')
    OPEN(UNIT=9,FILE='paramMod.F90')
    REWIND 9
    write(9,'(A16)') 'module param'
    write(9,'(A18)') '  implicit none'
    write(9,'(A22 A4 I9 A5 I9 A6 I9 A5 I9 A6 I9)') 'integer, parameter :: ', &
                    'NXA=', NI, &
                    ',NYA=', NJ, &
                    ',NXYA=', NIJ, &
                    ',NWA=', NWALI, &
                    ',PREC=',PREC
    write(9,'(A)') 'end module param'

    write(3,*) (ITB(1,I),I=1,NI), (ITB(2,I),I=1,NI),&
            (JTB(1,J),J=1,NJ), (JTB(2,J),J=1,NJ),& 
            (LI(I),I=1,NI),(CTD(IJP),IJP=0,NIJ-1),(DTC(I),I=1,NIJ),&
            (IJW(I),I=1,NWALI), (IJPW(I),I=1,NWALI),(IJW1(I),I=1,NWALI),&
            (IJW2(I),I=1,NWALI)
    write(3,*) (X(I),I=1,NIJ), (Y(I),I=1,NIJ), (XC(I),I=1,NIJ),&
            (YC(I),I=1,NIJ), (FX(I), I=1,NIJ), (FY(I), I=1,NIJ),DX,DY, VOL,&
            (SRDW(I),I=1,NWALI),(XTW(I),I=1,NWALI),(YTW(I),I=1,NWALI)

!...Create .vtk file
    print *, ' *** GENERATING .VTK *** '
    open (unit=4,FILE='grid.vtk')
    write(4,'(A)') '# vtk DataFile Version 3.0'
    write(4,'(A)') 'grid'
    write(4,'(A)') 'ASCII'
    write(4,'(A)') 'DATASET STRUCTURED_GRID'
    write(4,'(A I6 I6 I6)') 'DIMENSIONS', NIM, NJM,1
    write(4,'(A I10 A)') 'Points ', NIM*NJM, ' float'
    DO I=1,NIM
        DO J=1,NJM
            IJ=LI(I)+J
            write(4,'(E20.10,1X,E20.10,1X,E20.10)'), X(IJ), Y(IJ),0.0
        END DO
    END DO
    

end subroutine gridexport

!########################################################
subroutine setBc
!########################################################

    use bc
    use ind
    implicit none

    do L=1,4
        if(L.LE.2) then
            do I=2,NIM
                ITB(L,I)=BTYP(L)
            end do
        else
            do J=2,NJM
                JTB(L-2,J)=BTYP(L)
            end do
        end if
    end do

    call defbc(1, NWALI, IJW, IJPW, IJW1, IJW2)

end subroutine setBc

!########################################################
subroutine defBc(LT,NBCF,IJBB,IJBP, IJ1, IJ2)
!########################################################

    use bc
    use ind
    implicit none
    integer, intent(in) :: LT
    integer, dimension(N), intent(inout) :: IJBB, IJBP, IJ1, IJ2
    integer, intent(out) :: NBCF
!
!.....COLLECT BOUNDARY CELL FACES OF TYPE 'LT' IN A LIST FOR EACH GRID
!
      NBCF=0
!
!.....SOUTH SIDE
!
      DO I=2,NIM
        IF(ITB(1,I).EQ.LT) THEN
          NBCF=NBCF+1
          IJBB(NBCF)=LI(I)+1
          IJBP(NBCF)=IJBB(NBCF)+1
          IJ1(NBCF)=IJBB(NBCF)
          IJ2(NBCF)=IJ1(NBCF)-NJ
        ENDIF
      END DO
!
!.....NORTH SIDE
!
      DO I=2,NIM
        IF(ITB(2,I).EQ.LT) THEN
          NBCF=NBCF+1
          IJBB(NBCF)=LI(I)+NJ
          IJBP(NBCF)=IJBB(NBCF)-1
          IJ2(NBCF)=IJBP(NBCF)
          IJ1(NBCF)=IJ2(NBCF)-NJ
        ENDIF
      END DO
!
!.....WEST SIDE
!
      DO J=2,NJM
        IF(JTB(1,J).EQ.LT) THEN
          NBCF=NBCF+1
          IJBB(NBCF)=LI(1)+J
          IJBP(NBCF)=IJBB(NBCF)+NJ
          IJ2(NBCF)=IJBB(NBCF)
          IJ1(NBCF)=IJ2(NBCF)-1
        ENDIF
      END DO
!
!.....EAST SIDE
!
      DO J=2,NJM
        IF(JTB(2,J).EQ.LT) THEN
          NBCF=NBCF+1
          IJBB(NBCF)=LI(NI)+J
          IJBP(NBCF)=IJBB(NBCF)-NJ
          IJ1(NBCF)=IJBP(NBCF)
          IJ2(NBCF)=IJ1(NBCF)-1
        ENDIF
      END DO
!
end subroutine defBc

!########################################################
subroutine calcG
!########################################################

    use bc
    use geo
    use ind
    implicit none
    real*8 :: XE, XN,YE, YN, DLPE, DLEE, DLPN, DLNN, DN, AR,SMALL
    integer :: IJ1, IJ2

    SMALL=1.0E-20


!.....CALCULATION OF NODE COORDINATES: CORNER (DUMMY) NODES
!
        IJ=LI(1)+1
        XC(IJ)=X(IJ)
        YC(IJ)=Y(IJ)
!
        IJ=LI(NIM)+1
        XC(IJ+NJ)=X(IJ)
        YC(IJ+NJ)=Y(IJ)
!
        IJ=LI(1)+NJM
        XC(IJ+1)=X(IJ)
        YC(IJ+1)=Y(IJ)
!
        IJ=LI(NIM)+NJM
        XC(IJ+NJ+1)=X(IJ)
        YC(IJ+NJ+1)=Y(IJ)
!
!.....CALCULATION OF NODE COORDINATES: BOUNDARY NODES
!
        DO J=2,NJM
          IJ=LI(NIM)+J
          XC(IJ+NJ)=0.5d0*(X(IJ)+X(IJ-1))
          YC(IJ+NJ)=0.5d0*(Y(IJ)+Y(IJ-1))
          IJ=LI(1)+J
          XC(IJ)=0.5d0*(X(IJ)+X(IJ-1))
          YC(IJ)=0.5d0*(Y(IJ)+Y(IJ-1))
        END DO
!
        DO I=2,NIM
          IJ=LI(I)+1
          XC(IJ)=0.5d0*(X(IJ)+X(IJ-NJ))
          YC(IJ)=0.5d0*(Y(IJ)+Y(IJ-NJ))
          IJ=LI(I)+NJM
          XC(IJ+1)=0.5d0*(X(IJ)+X(IJ-NJ))
          YC(IJ+1)=0.5d0*(Y(IJ)+Y(IJ-NJ))
        END DO
!
!.....CALCULATION OF NODE COORDINATES: CELL CENTERS
!
        DO I=2,NIM
        DO J=2,NJM
          IJ=LI(I)+J
          XC(IJ)=0.25d0*(X(IJ)+X(IJ-1)+X(IJ-NJ)+X(IJ-NJ-1))
          YC(IJ)=0.25d0*(Y(IJ)+Y(IJ-1)+Y(IJ-NJ)+Y(IJ-NJ-1))
        END DO
        END DO
!
!......CALCULATION OF INTERPOLATION FACTORS
!
        DO I=2,NIM
        DO J=2,NJM
          IJ=LI(I)+J
!
!.....INTERPOLATION IN I-DIRECTION: FX = Pe/PE
!
          XE=0.5d0*(X(IJ)+X(IJ-1))
          YE=0.5d0*(Y(IJ)+Y(IJ-1))
          DLPE=SQRT((XE-XC(IJ))**2+(YE-YC(IJ))**2)
          DLEE=SQRT((XC(IJ+NJ)-XE)**2+(YC(IJ+NJ)-YE)**2)
          FX(IJ)=DLPE/(DLPE+DLEE+SMALL)
!
!.....INTERPOLATION IN J-DIRECTION: FY = Pn/PN
!
          XN=0.5d0*(X(IJ)+X(IJ-NJ))
          YN=0.5d0*(Y(IJ)+Y(IJ-NJ))
          DLPN=SQRT((XN-XC(IJ))**2+(YN-YC(IJ))**2)
          DLNN=SQRT((XC(IJ+1)-XN)**2+(YC(IJ+1)-YN)**2)
          FY(IJ)=DLPN/(DLPN+DLNN+SMALL)
        END DO
        END DO
!
!....Normal distance from cell face center to cell center
!
    do IW=1,NWALI
        IJB=IJW(IW)
        IJP=IJPW(IW)
        IJ1=IJW1(IW)
        IJ2=IJW2(IW)
        DX=X(IJ1)-X(IJ2)
        DY=Y(IJ1)-Y(IJ2)
        AR=sqrt(DX**2+DY**2)
        XTW(IW)=DX/(AR+SMALL)
        YTW(IW)=DY/(AR+SMALL)
        XN=YTW(IW)
        YN=-XTW(IW)
        DN=(XC(IJB)-XC(IJP))*XN+(YC(IJB)-YC(IJP))*YN
        SRDW(IW)=AR/(DN+SMALL)
    end do

end subroutine calcG

