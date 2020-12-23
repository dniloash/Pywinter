


function crear_int0(NX,NY,STARTLAT,STARTLON,DELTALAT,DELTALON,int1,fci,int2,int3,flv,NVARS,NS,va)

    integer, parameter :: IUNIT = 10
    integer, parameter :: OUNIT = 11

    integer :: NVARS
    character(len=200) :: NS


    integer :: VERSION = 5
    character(len=24) :: HDATE

    real :: XFCST = 0

    character(len=8) :: STARTLOC = 'SWCORNER'

    character(len=9) :: FIELD
    character(len=10) :: int1
    character(len=9), dimension(NVARS) :: fint

    character(len=24) :: fci


    character(len=25) :: UNITS
    character(len=10) :: int2
    character(len=25), dimension(NVARS) :: funt

    character(len=46) :: DESC
    character(len=10) :: int3
    character(len=46), dimension(NVARS) :: fdes

    character(len=32):: MAP_SOURCE = 'PYWINTER'

    real :: XLVL
    real, dimension(NVARS) :: flv

    real,dimension(NVARS,NX,NY) :: va


!---------------------------------

    integer :: NX
    integer :: NY
    integer :: IPROJ = 0
    real :: STARTLAT
    real :: STARTLON
    real :: DELTALAT
    real :: DELTALON
    real :: EARTH_RADIUS = 6367.470215
    logical :: IS_WIND_EARTH_REL = .FALSE.

    open(1,file=int1)
    open(2,file=int2)
    open(3,file=int3)

    open (IUNIT, CONVERT='BIG_ENDIAN', FILE= NS, form='unformatted', status='replace')

    HDATE = fci


    ciclo : do n = 1,NVARS


        read(1,'(a)') fint(n)
        read(2,'(a)') funt(n)
        read(3,'(a)') fdes(n)


        FIELD = fint(n)
        UNITS = funt(n)
        DESC = fdes(n)

        XLVL = flv(n)


        write (IUNIT, IOSTAT=IERR) VERSION
        write (IUNIT) HDATE, XFCST, MAP_SOURCE, FIELD, UNITS, DESC, XLVL, NX, NY, IPROJ
        write (IUNIT) STARTLOC, STARTLAT, STARTLON, DELTALAT, DELTALON, EARTH_RADIUS
        write (IUNIT) IS_WIND_EARTH_REL
        write (IUNIT) va(n,1:NX,1:NY)


   end do ciclo


   close (IUNIT)

   close(1)
   close(2)
   close(3)

end function crear_int0


!----------------------------------------------------------------------------------------------


function crear_int1(NX,NY,STARTLAT,STARTLON,DX,DY,TRUELAT1,int1,fci,int2,int3,flv,NVARS,NS,va)

    integer, parameter :: IUNIT = 10
    integer, parameter :: OUNIT = 11

    integer :: NVARS
    character(len=200) :: NS


    integer :: VERSION = 5
    character(len=24) :: HDATE

    real :: XFCST = 0

    character(len=8) :: STARTLOC = 'SWCORNER'

    character(len=9) :: FIELD
    character(len=10) :: int1
    character(len=9), dimension(NVARS) :: fint

    character(len=24) :: fci


    character(len=25) :: UNITS
    character(len=10) :: int2
    character(len=25), dimension(NVARS) :: funt

    character(len=46) :: DESC
    character(len=10) :: int3
    character(len=46), dimension(NVARS) :: fdes

    character(len=32):: MAP_SOURCE = 'PYWINTER'

    real :: XLVL
    real, dimension(NVARS) :: flv

    real,dimension(NVARS,NX,NY) :: va

!---------------------------------

    integer :: NX
    integer :: NY
    integer :: IPROJ = 1
    real :: STARTLAT
    real :: STARTLON
    real :: DX
    real :: DY
    real :: XLONC
    real :: TRUELAT1
    real :: TRUELAT2
    real :: EARTH_RADIUS = 6367.470215
    IS_WIND_EARTH_REL = .FALSE.

    open(1,file=int1)
    open(2,file=int2)
    open(3,file=int3)

    open (IUNIT, CONVERT='BIG_ENDIAN', FILE= NS, form='unformatted', status='replace')

    HDATE = fci


    ciclo : do n = 1,NVARS


        read(1,'(a)') fint(n)
        read(2,'(a)') funt(n)
        read(3,'(a)') fdes(n)


        FIELD = fint(n)
        UNITS = funt(n)
        DESC = fdes(n)

        XLVL = flv(n)


        write (IUNIT, IOSTAT=IERR) VERSION
        write (IUNIT) HDATE, XFCST, MAP_SOURCE, FIELD, UNITS, DESC, XLVL, NX, NY, IPROJ
        write (IUNIT) STARTLOC, STARTLAT, STARTLON, DX, DY, TRUELAT1, EARTH_RADIUS
        write (IUNIT) IS_WIND_EARTH_REL
        write (IUNIT) va(n,1:NX,1:NY)


   end do ciclo

   close (IUNIT)

   close(1)
   close(2)
   close(3)

end function crear_int1


!----------------------------------------------------------------------------------------------

function crear_int3(NX,NY,STARTLAT,STARTLON,DX,DY,XLONC,TRUELAT1,TRUELAT2,iswin,int1,fci,int2,int3,flv,NVARS,NS,va)

    integer, parameter :: IUNIT = 10
    integer, parameter :: OUNIT = 11

    integer :: NVARS
    character(len=200) :: NS


    integer :: VERSION = 5
    character(len=24) :: HDATE

    real :: XFCST = 0

    character(len=8) :: STARTLOC = 'SWCORNER'

    character(len=9) :: FIELD
    character(len=10) :: int1
    character(len=9), dimension(NVARS) :: fint

    character(len=24) :: fci


    character(len=25) :: UNITS
    character(len=10) :: int2
    character(len=25), dimension(NVARS) :: funt

    character(len=46) :: DESC
    character(len=10) :: int3
    character(len=46), dimension(NVARS) :: fdes

    character(len=32):: MAP_SOURCE = 'PYWINTER'

    logical :: iswin

    real :: XLVL
    real, dimension(NVARS) :: flv

    real,dimension(NVARS,NX,NY) :: va

!---------------------------------

    integer :: NX
    integer :: NY
    integer :: IPROJ = 3
    real :: STARTLAT
    real :: STARTLON
    real :: DX
    real :: DY
    real :: XLONC
    real :: TRUELAT1
    real :: TRUELAT2
    real :: EARTH_RADIUS = 6367.470215
    IS_WIND_EARTH_REL = iswin

    open(1,file=int1)
    open(2,file=int2)
    open(3,file=int3)

    open (IUNIT, CONVERT='BIG_ENDIAN', FILE= NS, form='unformatted', status='replace')

    HDATE = fci


    ciclo : do n = 1,NVARS


        read(1,'(a)') fint(n)
        read(2,'(a)') funt(n)
        read(3,'(a)') fdes(n)


        FIELD = fint(n)
        UNITS = funt(n)
        DESC = fdes(n)

        XLVL = flv(n)


        write (IUNIT, IOSTAT=IERR) VERSION
        write (IUNIT) HDATE, XFCST, MAP_SOURCE, FIELD, UNITS, DESC, XLVL, NX, NY, IPROJ
        write (IUNIT) STARTLOC, STARTLAT, STARTLON, DX, DY, XLONC, TRUELAT1, TRUELAT2, EARTH_RADIUS
        write (IUNIT) IS_WIND_EARTH_REL
        write (IUNIT) va(n,1:NX,1:NY)


   end do ciclo

   close (IUNIT)

   close(1)
   close(2)
   close(3)

end function crear_int3


!----------------------------------------------------------------------------------------------


function crear_int4(NX,NY,STARTLAT,STARTLON,NLATS,DELTALON,iswin,int1,fci,int2,int3,flv,NVARS,NS,va)

    integer, parameter :: IUNIT = 10
    integer, parameter :: OUNIT = 11

    integer :: NVARS
    character(len=200) :: NS


    integer :: VERSION = 5
    character(len=24) :: HDATE

    real :: XFCST = 0

    character(len=8) :: STARTLOC = 'SWCORNER'

    character(len=9) :: FIELD
    character(len=10) :: int1
    character(len=9), dimension(NVARS) :: fint

    character(len=24) :: fci


    character(len=25) :: UNITS
    character(len=10) :: int2
    character(len=25), dimension(NVARS) :: funt

    character(len=46) :: DESC
    character(len=10) :: int3
    character(len=46), dimension(NVARS) :: fdes

    character(len=32):: MAP_SOURCE = 'PYWINTER'

    logical :: iswin

    real :: XLVL
    real, dimension(NVARS) :: flv

    real,dimension(NVARS,NX,NY) :: va

!---------------------------------

    integer :: NX
    integer :: NY
    integer :: IPROJ = 4
    real :: STARTLAT
    real :: STARTLON
    real :: NLATS
    real :: DELTALON

    real :: EARTH_RADIUS = 6367.470215
    IS_WIND_EARTH_REL = iswin

    open(1,file=int1)
    open(2,file=int2)
    open(3,file=int3)

    open (IUNIT, CONVERT='BIG_ENDIAN', FILE= NS, form='unformatted', status='replace')

    HDATE = fci


    ciclo : do n = 1,NVARS


        read(1,'(a)') fint(n)
        read(2,'(a)') funt(n)
        read(3,'(a)') fdes(n)


        FIELD = fint(n)
        UNITS = funt(n)
        DESC = fdes(n)

        XLVL = flv(n)


        write (IUNIT, IOSTAT=IERR) VERSION
        write (IUNIT) HDATE, XFCST, MAP_SOURCE, FIELD, UNITS, DESC, XLVL, NX, NY, IPROJ
        write (IUNIT) STARTLOC, STARTLAT, STARTLON, NLATS, DELTALON, EARTH_RADIUS
        write (IUNIT) IS_WIND_EARTH_REL
        write (IUNIT) va(n,1:NX,1:NY)


   end do ciclo

   close (IUNIT)

   close(1)
   close(2)
   close(3)

end function crear_int4

!----------------------------------------------------------------------------------------------


function crear_int5(NX,NY,STARTLAT,STARTLON,DX,DY,XLONC,TRUELAT1,iswin,int1,fci,int2,int3,flv,NVARS,NS,va)

    integer, parameter :: IUNIT = 10
    integer, parameter :: OUNIT = 11

    integer :: NVARS
    character(len=200) :: NS


    integer :: VERSION = 5
    character(len=24) :: HDATE

    real :: XFCST = 0

    character(len=8) :: STARTLOC = 'SWCORNER'

    character(len=9) :: FIELD
    character(len=10) :: int1
    character(len=9), dimension(NVARS) :: fint

    character(len=24) :: fci


    character(len=25) :: UNITS
    character(len=10) :: int2
    character(len=25), dimension(NVARS) :: funt

    character(len=46) :: DESC
    character(len=10) :: int3
    character(len=46), dimension(NVARS) :: fdes

    character(len=32):: MAP_SOURCE = 'PYWINTER'

    logical :: iswin

    real :: XLVL
    real, dimension(NVARS) :: flv

    real,dimension(NVARS,NX,NY) :: va

!---------------------------------

    integer :: NX
    integer :: NY
    integer :: IPROJ = 5
    real :: STARTLAT
    real :: STARTLON
    real :: DX
    real :: DY
    real :: XLONC
    real :: TRUELAT1
    real :: EARTH_RADIUS = 6367.470215
    IS_WIND_EARTH_REL = iswin

    open(1,file=int1)
    open(2,file=int2)
    open(3,file=int3)

    open (IUNIT, CONVERT='BIG_ENDIAN', FILE= NS, form='unformatted', status='replace')

    HDATE = fci


    ciclo : do n = 1,NVARS


        read(1,'(a)') fint(n)
        read(2,'(a)') funt(n)
        read(3,'(a)') fdes(n)


        FIELD = fint(n)
        UNITS = funt(n)
        DESC = fdes(n)

        XLVL = flv(n)


        write (IUNIT, IOSTAT=IERR) VERSION
        write (IUNIT) HDATE, XFCST, MAP_SOURCE, FIELD, UNITS, DESC, XLVL, NX, NY, IPROJ
        write (IUNIT) STARTLOC, STARTLAT, STARTLON, DX, DY, XLONC, TRUELAT1, EARTH_RADIUS
        write (IUNIT) IS_WIND_EARTH_REL
        write (IUNIT) va(n,1:NX,1:NY)


   end do ciclo

   close (IUNIT)

   close(1)
   close(2)
   close(3)

end function crear_int5