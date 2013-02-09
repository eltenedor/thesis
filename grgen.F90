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

    PRINT *, ' ENTER> ZSTART, ZEND, NUMBER OF CVS:  '
    IF(ITYP.EQ.1) THEN
        READ(*,*) ZZS,ZZE,NKCV
        WRITE(1,*) ZZS,ZZE,NKCV,'   ZS,ZE,NKCV '
    ELSE
        READ(2,*) ZZS,ZZE,NKCV
    END IF

    PRINT *, ' ENTER> BOUNDARY TYPE N S W E T B:  '
    IF(ITYP.EQ.1) THEN
        READ(*,*) BTYP(1:6)
        WRITE(1,*) BTYP(1:6),  '   SBTYP, NBTYP, WBTYP, EBTYP, TBTYP, BBTYP '
    ELSE
        READ(2,*) BTYP(1:6)
    END IF
       
end subroutine readData

!########################################################
subroutine cartesian
!########################################################

    use geo
    use ind
    implicit none

    DX=(XXE-XXS)/dble(NICV)
    DY=(YYE-YYS)/dble(NJCV)
    DZ=(ZZE-ZZS)/dble(NKCV)
    VOL=DX*DY*DZ

    NI=NICV+2
    NJ=NJCV+2
    NK=NKCV+2
    NIJ=NI*NJ
    NIM=NI-1
    NJM=NJ-1
    NKM=NK-1
    IJST=1
    N=NICV*NJCV*NKCV

    do I=1,NI
        LI(I)=(I-1)*NJ
    end do

    do K=1,NK
        LK(K)=(K-1)*NJ*NI
    end do

    call createMapping
    
    !print *, YYS
    do K=1,NKM
        do I=1,NIM
           do J=1,NJM
                IJK=LK(K)+LI(I)+J
                X(IJK)=XXS+dble(I-1)*DX
                Y(IJK)=YYS+dble(J-1)*DY 
                Z(IJK)=ZZS+dble(K-1)*DZ
                !print *, IJ, dble(J-1), Y(IJ), YYS, DY
            end do
        end do
    end do

end subroutine cartesian

!########################################################
subroutine createMapping
!########################################################

    use ind
    use geo
    implicit none

    ! PETSc routines indices start from O!!
    IJKP=-1
    do K=2,NKM
        do I=2,NIM
            do J=2,NJM
                IJKP=IJKP+1
                IJK=LK(K)+LI(I)+J
                CTD(IJKP)=IJK
                DTC(IJK)=IJKP
            end do
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

    DX=(XXE-XXS)/dble(NICV)
    DY=(YYE-YYS)/dble(NJCV)
    DZ=(ZZE-ZZS)/dble(NKCV)
!...Create solver file
    !OPEN(UNIT=9,FILE='../../../pet_src/paramMod.F90')
    OPEN(UNIT=9,FILE='paramMod.F90')
    REWIND 9
    write(9,'(A16)') 'module param'
    write(9,'(A18)') '  implicit none'
    write(9,'(A22 A4 I9 A5 I9 A5 I9 A6 I9 A5 I9 A6 I9)') 'integer, parameter :: ', &
                    'NXA=', NI, &
                    ',NYA=', NJ, &
                    ',NZA=', NK, &
                    ',NXYZA=', NIJK, &
                    ',NWA=', NWALI, &
                    ',PREC=',PREC
    write(9,'(A)') 'end module param'

    write(3,*)  (ITB(1,I),I=1,NI),(ITB(2,I),I=1,NI),&
                (JTB(1,J),J=1,NJ),(JTB(2,J),J=1,NJ),& 
                (KTB(1,K),K=1,NK),(KTB(2,K),K=1,NK),&
                (LK(K),I=1,NK),(LI(I),I=1,NI),&
                (CTD(I),I=0,NIJK-1),(DTC(I),I=1,NIJK),&
                (IJKW(I),I=1,NWALI),(IJKPW(I),I=1,NWALI),&
                (IJKW1(I),I=1,NWALI),(IJKW2(I),I=1,NWALI)
    write(3,*)  (X(I),I=1,NIJK),(Y(I),I=1,NIJK),(Z(I),I=1,NIJK),&
                (XC(I),I=1,NIJK),(YC(I),I=1,NIJK),(ZC(I),I=1,NIJK),&
                (FX(I), I=1,NIJK),(FY(I), I=1,NIJK),(FZ(I), I=1,NIJK),&
                DX,DY,DZ, VOL,&
                (SRDW(I),I=1,NWALI)

!...Create .vtk file
    print *, ' *** GENERATING .VTK *** '
    open (unit=4,FILE='grid.vtk')
    write(4,'(A)') '# vtk DataFile Version 3.0'
    write(4,'(A)') 'grid'
    write(4,'(A)') 'ASCII'
    write(4,'(A)') 'DATASET STRUCTURED_GRID'
    write(4,'(A I6 I6 I6)') 'DIMENSIONS', NIM, NJM,NKM
    write(4,'(A I10 A)') 'Points ', NIM*NJM*NKM, ' float'
    do K=1,NKM
        DO I=1,NIM
            DO J=1,NJM
                IJK=LK(K)+LI(I)+J
                write(4,'(E20.10,1X,E20.10,1X,E20.10)'), X(IJK), Y(IJK),Z(IJK)
            END DO
        END DO
    end do

end subroutine gridexport

!########################################################
subroutine setBc
!########################################################

    use bc
    use ind
    implicit none

    do L=1,6
        if(L.LE.2) then
            do I=2,NIM
                ITB(L,I)=BTYP(L)
            end do
        elseif(L.LE.4) then
            do J=2,NJM
                JTB(L-2,J)=BTYP(L)
            end do
        else
            do K=2,NKM
                KTB(L-4,K)=BTYP(L)
            end do
        end if
    end do

    call defbc(1, NWALI, IJKW, IJKPW, IJKW1, IJKW2, IJKW3, IJKW4)

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
    do K=2,NKM
        do I=2,NIM
            IF(ITB(1,I).EQ.LT) THEN
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
    do K=2,NKM
        do I=2,NIM
            IF(ITB(2,I).EQ.LT) THEN
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
    do K=2,NKM
        do J=2,NJM
            IF(JTB(1,J).EQ.LT) THEN
            NBCF=NBCF+1
            IJKBB(NBCF)=LK(K)+LI(1)+J
            IJKBP(NBCF)=IJKBB(NBCF)+NJ
            IJK3(NBCF)=IJKBB(NBCF)
            IJK4(NBCF)=IJK2(NBCF)-1
            IJK2(NBCF)=IJK3(NBCF)-NIJ
            IJK1(NBCF)=IJK4(NBCF)-NIJ
        ENDIF
      END DO
!
!.....EAST SIDE
!
        do J=2,NJM
            IF(JTB(2,J).EQ.LT) THEN
                NBCF=NBCF+1
                IJKBB(NBCF)=LK(K)+LI(NI)+J
                IJKBP(NBCF)=LK(K)+IJKBB(NBCF)-NJ
                IJK3(NBCF)=IJKBP(NBCF)
                IJK4(NBCF)=IJK3(NBCF)-1
                IJK2(NBCF)=IJK3(NBCF)-NIJ
                IJK1(NBCF)=IJK4(NBCF)-NIJ
            end if
        end do
    end do
!
!......BOTTOM SIDE
!
    do I=2,NIM
        do J=2,NJM
            if (KTB(1,I).EQ.LT) then
                NBCF=NBCF+1
                IJKBB(NBCF)=LK(1)+LI(I)+J
                IJKBP(NBCF)=IJKBB(NBCF)+NIJ
                IJK1(NBCF)=IJKBB(NBCF)-1
                IJK2(NBCF)=IJK1(NBCF)-NJ
                IJK3(NBCF)=IJK2(NBCF)+1
                IJK4(NBCF)=IJK1(NBCF)+1
            end if
        end do
    end do
!
!......TOP SIDE
!
    do I=2,NIM
        do J=2,NJM
            if (KTB(2,I).EQ.LT) then
                NBCF=NBCF+1
                IJKBB(NBCF)=LK(NK)+LI(I)+J
                IJKBP(NBCF)=IJKBB(NBCF)-NIJ
                IJK2(NBCF)=IJKBP(NBCF)-1
                IJK1(NBCF)=IJK2(NBCF)-NJ
                IJK3(NBCF)=IJK2(NBCF)+1
                IJK4(NBCF)=IJK1(NBCF)+1
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
    real*8 ::   DLPE,DLEE,DLPN,DLNN,DLPT,DLTT,DN,AR,SMALL,&
                XE,YE,ZE,XT,YT,ZT,&
                X1,X2,X4,Y1,Y2,Y4,Z1,Z2,Z4
    integer :: IJK1, IJK2, IJK3, IJK4

    SMALL=1.0E-20


!.....CALCULATION OF NODE COORDINATES: CORNER (DUMMY) NODES
!
    do K=1,NK,NK-1
        IJK=LK(K)+LI(1)+1
        XC(IJK)=X(IJK)
        YC(IJK)=Y(IJK)
        ZC(IJK)=Z(IJK)
!
        IJK=LK(K)+LI(NIM)+1
        XC(IJK+NJ)=X(IJK)
        YC(IJK+NJ)=Y(IJK)
        ZC(IJK+NJ)=Z(IJK)
!
        IJK=LK(K)+LI(1)+NJM
        XC(IJK+1)=X(IJK)
        YC(IJK+1)=Y(IJK)
        ZC(IJK+1)=Z(IJK)
!
        IJK=LK(K)+LI(NIM)+NJM
        XC(IJK+NJ+1)=X(IJK)
        YC(IJK+NJ+1)=Y(IJK)
        ZC(IJK+NJ+1)=Z(IJK)
    end do
!
!.....CALCULATION OF NODE COORDINATES: BOUNDARY NODES
!
    do K=2,NKM
        DO J=2,NJM
          !EAST
          IJK=LK(K)+LI(NIM)+J
          XC(IJK+NJ)=0.5d0*(X(IJK)+X(IJK-1))
          YC(IJK+NJ)=0.5d0*(Y(IJK)+Y(IJK-1))
          ZC(IJK+NJ)=0.5d0*(Z(IJK)+Z(IJK-NIJ))
          !WEST
          IJK=LK(K)+LI(1)+J
          XC(IJK)=0.5d0*(X(IJK)+X(IJK-1))
          YC(IJK)=0.5d0*(Y(IJK)+Y(IJK-1))
          ZC(IJK)=0.5d0*(Z(IJK)+Z(IJK-NIJ))
        END DO
    end do
!
    do K=2,NKM
        DO I=2,NIM
          !SOUTH
          IJK=LK(K)+LI(I)+1
          XC(IJK)=0.5d0*(X(IJK)+X(IJK-NJ))
          YC(IJK)=0.5d0*(Y(IJK)+Y(IJK-NJ))
          ZC(IJK)=0.5d0*(Z(IJK)+Z(IJK-NIJ))
          !NORTH
          IJK=LK(K)+LI(I)+NJM
          XC(IJK+1)=0.5d0*(X(IJK)+X(IJK-NJ))
          YC(IJK+1)=0.5d0*(Y(IJK)+Y(IJK-NJ))
          ZC(IJK+1)=0.5d0*(Z(IJK)+Z(IJK-NIJ))
        END DO
    end do

    do I=2,NIM
        do J=2,NJM
            !BOTTOM
            IJK=LK(1)+LI(I)+J
            XC(IJK)=0.5d0*(X(IJK)+X(IJK-NJ))
            YC(IJK)=0.5d0*(Y(IJK)+Y(IJK-NJ))
            ZC(IJK)=Z(IJK)
            !TOP
            IJK=LK(NK)+LI(I)+J
            XC(IJK)=0.5d0*(X(IJK)+X(IJK-NJ))
            YC(IJK)=0.5d0*(Y(IJK)+Y(IJK-NJ))
            ZC(IJK)=Z(IJK)
        end do
    end do

!
!.....CALCULATION OF NODE COORDINATES: CELL CENTERS
!
    do K=2,NKM
        DO I=2,NIM
        DO J=2,NJM
          IJK=LK(K)+LI(I)+J
          XC(IJK)=0.25d0*(X(IJK)+X(IJK-1)+X(IJK-NJ)+X(IJK-NJ-1))
          YC(IJK)=0.25d0*(Y(IJK)+Y(IJK-1)+Y(IJK-NJ)+Y(IJK-NJ-1))
          ZC(IJK)=0.5d0*(Z(IJK)+Z(IJK-NIJ))
        END DO
        END DO
    end do
!
!......CALCULATION OF INTERPOLATION FACTORS
!
    do K=2,NKM
        DO I=2,NIM
        DO J=2,NJM
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
          DLPN=SQRT((XN-XC(IJK))**2+(YN-YC(IJK))**2*(ZN-ZC(IJK)**2))
          DLNN=SQRT((XC(IJK+1)-XN)**2+(YC(IJK+1)-YN)**2+(ZC(IJK+1)-ZN)**2)
          FY(IJK)=DLPN/(DLPN+DLNN+SMALL)
!
!.....INTERPOLATION IN K-DIRECTION: FZ = Pt/PT
!
          XT=0.5d0*(X(IJK)+X(IJK-NJ))
          YT=0.5d0*(Y(IJK)+Y(IJK-1))
          ZT=0.5d0*(Z(IJK))
          DLPT=SQRT((XT-XC(IJK))**2+(YT-YC(IJK))**2*(ZT-ZC(IJK))**2)
          DLTT=SQRT((XC(IJK+1)-XT)**2+(YC(IJK+1)-YT)**2+(ZC(IJK+1)-ZT)**2)
          FZ(IJK)=DLPT/(DLPT+DLTT+SMALL)
        END DO
        END DO
    end do

!
!....Normal distance from cell face center to cell center
!
    do IW=1,NWALI
        IJKB=IJKW(IW)
        IJKP=IJKPW(IW)
        IJK1=IJKW1(IW)
        IJK2=IJKW2(IW)
        IJK3=IJKW3(IW)
        IJK4=IJKW4(IW)
        !
        X1=X(IJK1)
        Y1=Y(IJK1)
        Z1=Z(IJK1)
        !
        X2=X(IJK2)
        Y2=Y(IJK2)
        Z2=Z(IJK2)
        !
        X4=X(IJK4)
        Y4=Y(IJK4)
        Z4=Z(IJK4)
        !
        NX= (Y1-Y2)*(Z1-Z4)-(Y1-Y4)*(Z1-Z2)
        NY= -(X1-X2)*(Z1-Z4)+(X1-X4)*(Z1-Z2)
        NZ= (X1-X2)*(Y1-Y4)-(X1-X4)*(Y1-Y2)
        !
        AR = SQRT(ABS((X1-X2)*(Y1-Y4)-(X1-X4)*(Y1-Y2))**2+&
                  ABS(-(X1-X2)*(Z1-Z4)+(X1-X4)*(Z1-Z2))**2+&
                  ABS((Y1-Y2)*(Z1-Z4)-(Y1-Y4)*(Z1-Z2))**2)
        !
        DN =(XC(IJKB)-XC(IJKP))*NX/AR+&
            (YC(IJKB)-YC(IJKP))*NY/AR+&
            (ZC(IJKB)-ZC(IJKP))*NZ/AR
        SRDW(IW)=AR/(DN+SMALL)
    end do

end subroutine calcG

