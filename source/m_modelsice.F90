module m_modelsice

  use netcdf
  use cmor_users_functions
  use m_utilities
  use m_namelists

  implicit none

  ! Netcdf variables
  integer :: ncid, rhid, dimid, status

  ! Grid dimensions and variables
  integer, save                                     :: idm, jdm
  integer, parameter                                :: ncrns = 4
  real(kind=8), allocatable, save, dimension(:)     :: xvec, yvec
  real(kind=8), allocatable, save, dimension(:, :)  :: &
    angle, tlon, tlat, ulon, ulat, vlon, vlat, uvlon, uvlat, tarea, uarea, &
    tlon2, tlat2, ulon2, ulat2, vlon2, vlat2, uvlon2, uvlat2
  real(kind=8), allocatable, save, dimension(:, :, :)   :: &
    tlon_crns, tlat_crns, tlon_crnsp, tlat_crnsp, ulon_crns, ulat_crns, &
    ulon_crnsp, ulat_crnsp, vlon_crns, vlat_crns, vlon_crnsp, vlat_crnsp, &
    uvlon_crns, uvlat_crns, uvlon_crnsp, uvlat_crnsp
  character(len=slenmax), save                          :: tcoord, zcoord

  ! Fram Strait grid info
  integer                       :: seclen
  integer, parameter            :: maxseclen = 100
  integer, dimension(maxseclen) :: iind, jind, iflg, jflg
  logical, save                 :: lsecindex

  ! Dataset related variables
  character(len=slenmax), save  :: ivnm, ovnm, vunits, vpositive, vcomment

  ! Table related variables
  character(len=lenmax)         :: table

  ! String for module special
  character(len=slenmax), save  :: special

  ! Cmor parameters
  integer, save :: iaxid, jaxid, kaxid, taxid, grdid, varid, table_id, &
    table_id_grid, error_flag

  ! Data fields
  real(kind=8), allocatable, save, dimension(:, :)      :: fld, fld2, fld3, fldacc
  real(kind=8), allocatable, save, dimension(:, :, :)   :: fld3d

  ! Auxillary variables for special operations
  character(len=slenmax), save                          :: str1, str2

contains

  ! -----------------------------------------------------------------

  subroutine ice2cmor

    implicit none

    logical :: badrec, last
    integer :: m, n, nrec

    badrec = .false.

    ! Print start information
    if (verbose) then
      write(*, *)
      write(*, *) '------------------------------'
      write(*, *) '--- Process sea ice output ---'
      write(*, *) '------------------------------'
      write(*, *)
    end if

    ! Read grid information from input files
    write(*, *) 'Read grid information from input files'
    itag = tagimon
    call scan_files(reset=.true.)
    if (len_trim(fnm) == 0) return
    call read_gridinfo_ifile

    ! Read Fram Strait grid information from secindex.dat
    write(*, *) 'Read Fram Strait grid information from secindex.dat'
    call read_secindex(trim(griddata)//trim(secindexfile), lsecindex, &
      seclen, iind, jind, iflg, jflg)

    ! Process table fx
    write(*, *) 'Process table fx'
    fnm = pfx
    table = tfx
    do n = 1, nfx
      if (skip_variable(n, nfx, dfx)) cycle

      ! Map namelist variables
      ovnm = vfx(ovnmpos, n)
      ivnm = vfx(ivnmpos, n)
      special = vfx(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if vertical coordinate required
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)

      ! Check if input variable is present
      if (len_trim(pfx) == 0) cycle
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Skip if ouput variable is atmosphere cell area
      if (index(ovnm, 'areacella') > 0) cycle

      ! Prepare output file
      call special_pre
      call open_ofile(fx=.true.)

      ! Read field
      call read_field

      ! Post Processing
      call special_post

      ! Write field
      call write_field

      ! Close output file
      call close_ofile

    end do

    ! Process table omon
    write(*, *) 'Process table omon'
    fnm = pomon
    table = tomon
    do n = 1, nomon
      if (skip_variable(n, nomon, domon)) cycle

      ! Map namelist variables
      ovnm = vomon(ovnmpos, n)
      ivnm = vomon(ivnmpos, n)
      special = vomon(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Do not process ocean current fields uvel,vvel
      if (trim(ivnm) == 'uvel' .or. trim(ivnm) == 'vvel') exit

      ! Check if vertical coordinate required
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)

      ! Choose history file
      select case (trim(special))
      case ('day2mon')
        itag = tagiday
      case default
        itag = tagimon
      end select

      ! Check if input variable is present
      if (len_trim(pomon) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over input files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, romon) == 0) call open_ofile

        ! Choose history file
        select case (trim(special))
        case ('day2mon')
          itag = tagiday
        case default
          itag = tagimon
        end select

        ! Read variable into buffer (average if necessary)
        rec = 0
        nrec = 0
        fldacc = 0.
        last = .false.
        do
          if (len_trim(pomon) == 0) call scan_files(reset=.false.)
          if (rec == 0) then
            last = .true.
            exit
          end if
          nrec = nrec + 1
          call read_tslice(rec, badrec, fnm)
          fldacc = fldacc + fld
          if (tbnd(2) == mbnd(2)) exit
        end do
        if (last) exit
        fld = fldacc / real(nrec)
        tbnds(:, 1) = mbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Post processing
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, romon) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, romon) > 0) call close_ofile

    end do

    ! Process table oimon
    write(*, *) 'Process table oimon'
    fnm = poimon
    table = toimon
    do n = 1, noimon
      if (skip_variable(n, noimon, doimon)) cycle

      ! Map namelist variables
      ovnm = voimon(ovnmpos, n)
      ivnm = voimon(ivnmpos, n)
      special = voimon(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if vertical coordinate required
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)

      ! Choose history file
      select case (trim(special))
      case ('day2mon')
        itag = tagiday
      case default
        itag = tagimon
      end select

      ! Check if input variable is present
      if (len_trim(poimon) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over input files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, roimon) == 0) call open_ofile

        ! Choose history file
        select case (trim(special))
        case ('day2mon')
          itag = tagiday
        case default
          itag = tagimon
        end select

        ! Read variable into buffer (average if necessary)
        rec = 0
        nrec = 0
        fldacc = 0.
        last = .false.
        do
          if (len_trim(poimon) == 0) call scan_files(reset=.false.)
          if (rec == 0) then
            last = .true.
            exit
          end if
          nrec = nrec + 1
          call read_tslice(rec, badrec, fnm)
          fldacc = fldacc + fld
          if (tbnd(2) == mbnd(2)) exit
        end do
        if (last) exit
        fld = fldacc / real(nrec)
        tbnds(:, 1) = mbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Post processing
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, roimon) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, roimon) > 0) call close_ofile

    end do

    ! Process table day
    write(*, *) 'Process table day'
    fnm = pday
    table = tday
    do n = 1, nday
      if (skip_variable(n, nday, dday)) cycle

      ! Map namelist variables
      ovnm = vday(ovnmpos, n)
      ivnm = vday(ivnmpos, n)
      special = vday(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if vertical coordinate required
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)

      ! Choose history file
      itag = tagiday

      ! Check if input variable is present
      if (len_trim(pday) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over input files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rday) == 0) call open_ofile

        ! Read variable into buffer
        rec = 0
        if (len_trim(pday) == 0) call scan_files(reset=.false.)
        if (rec == 0) exit
        tbnds(:, 1) = tbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, rday) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rday) > 0) call close_ofile

    end do

    ! Process table SIday
    write(*, *) 'Process table SIday'
    fnm = pSIday
    table = tSIday
    do n = 1, nSIday
      if (skip_variable(n, nSIday, dSIday)) cycle

      ! Map namelist variables
      ovnm = vSIday(ovnmpos, n)
      ivnm = vSIday(ivnmpos, n)
      special = vSIday(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if vertical coordinate required
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)

      ! Choose history file
      itag = tagiday

      ! Check if input variable is present
      if (len_trim(pSIday) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over input files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rSIday) == 0) call open_ofile

        ! Read variable into buffer
        rec = 0
        if (len_trim(pSIday) == 0) call scan_files(reset=.false.)
        if (rec == 0) exit
        tbnds(:, 1) = tbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, rSIday) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rSIday) > 0) call close_ofile

    end do

  end subroutine ice2cmor

  ! -----------------------------------------------------------------

  subroutine special_pre

    implicit none

    str2 = special
    do
      if (index(str2, ';') > 0) then
        str1 = str2(1:index(str2, ';') - 1)
        str2 = str2(index(str2, ';') + 1:)
      else
        str1 = str2
      end if
      select case (str1)

        ! Fix unitless units
      case ('unitless')
        vunits = '1'

        ! Set correct units for percentage
      case ('percent')
        vunits = '%'

        ! Set correct units for percentage
      case ('fraction')
        vunits = '1'

        ! Unit transformation: kg m-2
      case ('kg m-2')
        vunits = 'kg m-2'

        ! Unit transformation: cm/day -> kg m-2 s-1
      case ('cmFW day-1 -> kg m-2 s-1', 'cmICE day-1 -> kg m-2 s-1', &
          'cmSNOW day-1 -> kg m-2 s-1')
        vunits = 'kg m-2 s-1'

        ! Unit transformation: kg/m^2/day -> kg m-2 s-1
      case ('kg m-2 day-1 -> kg m-2 s-1')
        vunits = 'kg m-2 s-1'

        ! Fix micrometers units
      case ('micrometer')
        vunits = 'micrometers'

        ! Convert degC to K
      case ('Celsius2Kelvin')
        vunits = 'K'

        ! Fix m-2 units
      case ('m-2')
        vunits = 'm-2'

        ! Convert units from radians2 to m2
      case ('rad2m')
        vunits = 'm2'

        ! Convert units from m3 to 1e3 km3
      case ('1e3 km3')
        vunits = '1e3 km3'

        ! Convert units from m2 to 1e6 km2
      case ('1e6 km2')
        vunits = '1e6 km2'

        ! Fix hcice units
      case ('J m-2')
        vunits = 'J m-2'

        ! Set positive attribute
      case ('positiveup')
        vpositive = 'up'
      case ('positivedo')
        vpositive = 'down'

        ! Write comment for hur and hurs
      case ('tsicecomment')
        vcomment = 'This field differs from the CMOR-definition ' &
          // 'because at every time-step, the field has been ' &
          // 'put to 271.314 K in grid-cells without sea ice. ' &
          // 'The time-mean has been calculated from these ' &
          // 'values without taking care of when a grid-cell is ' &
          // 'free of sea ice.'

        ! Write comment for streng and divice
      case ('tavecomment')
        vcomment = 'During the time-averaging there was no weighting ' &
          // 'with respect to the ice concentration in this ' &
          // 'field. A value of zero was used over open water.'

      end select
      if (str1 == str2) exit
    end do

  end subroutine special_pre

  ! -----------------------------------------------------------------

  subroutine special_post

    implicit none

    integer :: i, j

    str2 = special
    do
      if (index(str2, ';') > 0) then
        str1 = str2(1:index(str2, ';') - 1)
        str2 = str2(index(str2, ';') + 1:)
      else
        str1 = str2
      end if
      select case (str1)

        ! Convert units from radians2 to m2
      case ('rad2m')
        fld = fld * 6.37122e6**2

        ! Convert degC to K
      case ('Celsius2Kelvin')
        fld = fld + 273.15

        ! Unit transformation: cm day-1 -> kg m-2 s-1
      case ('cmFW day-1 -> kg m-2 s-1')
        fld = fld / 100. / (24. * 3600.) * 1000.
      case ('cmICE day-1 -> kg m-2 s-1')
        fld = fld / 100. / (24. * 3600.) * 917.
      case ('cmSNOW day-1 -> kg m-2 s-1')
        fld = fld / 100. / (24. * 3600.) * 330.

        ! Unit transformation: kg/m^2/day -> kg m-2 s-1
      case ('kg m-2 day-1 -> kg m-2 s-1')
        fld = fld / (24. * 3600.)

        ! Flip sign
      case ('flipsign')
        do j = 1, jdm
          do i = 1, idm
            if (fld(i, j) < 1e20) fld(i, j) = -fld(i, j)
          end do
        end do

        ! Set ice free points to missing value
      case ('zero2missing')
        do j = 1, jdm
          do i = 1, idm
            if (abs(fld(i, j)) < 1e-6) fld(i, j) = 1e20
          end do
        end do

        ! Divide by cell area
      case ('Xcellarea-1')
        do j = 1, jdm
          do i = 1, idm
            if (abs(fld(i, j)) < 1e20) fld(i, j) = fld(i, j) / tarea(i, j)
          end do
        end do

        ! mask southern hemisphere
      case ('masksh')
        do j = 1, jdm
          do i = 1, idm
            if (tlat(i, j) < 0) fld(i, j) = 0.0
          end do
        end do

        ! mask northern hemisphere
      case ('masknh')
        do j = 1, jdm
          do i = 1, idm
            if (tlat(i, j) > 0) fld(i, j) = 0.0
          end do
        end do

        ! mask sea ice extent less than 15%
      case ('mask15p')
        do j = 1, jdm
          do i = 1, idm
            if (fld(i, j) < 0.15) fld(i, j) = 0.0
          end do
        end do

        ! Multiple by cell area
      case ('Xcellarea')
        do j = 1, jdm
          do i = 1, idm
            if (abs(fld(i, j)) < 1e20) fld(i, j) = fld(i, j) * tarea(i, j)
          end do
        end do

        ! Global sum
      case ('glbsum')
        do j = 1, jdm
          do i = 1, idm
            if (abs(fld(i, j)) > 1e20) fld(i, j) = 0.0
          end do
        end do
        fldacc(1, 1) = sum(fld)
        fld = fldacc(1, 1)

        ! convert to 1e3 km3
      case ('1e3 km3')
        fld = fld / 1e9

        ! convert to 1e6 km2
      case ('1e6 km2')
        fld = fld / 1e12

      end select
      if (str1 == str2) exit
    end do

  end subroutine special_post

  ! -----------------------------------------------------------------

  subroutine read_gridinfo_ifile

    implicit none

    logical         :: check
    integer         :: i, j, n
    real(kind=8)    :: missing, theta, lambda

    ! Open first input file
    call scan_files(reset=.true.)
    status = nf90_open(fnm, nf90_nowrite, ncid)
    call handle_ncerror(status)

    ! Read longitudes and latitudes
    status = nf90_inq_dimid(ncid, 'ni', dimid)
    call handle_ncerror(status)
    status = nf90_inquire_dimension(ncid, dimid, len=idm)
    call handle_ncerror(status)
    status = nf90_inq_dimid(ncid, 'nj', dimid)
    call handle_ncerror(status)
    status = nf90_inquire_dimension(ncid, dimid, len=jdm)
    call handle_ncerror(status)

    allocate(tlon(idm, jdm), tlat(idm, jdm), ulon(idm, jdm), ulat(idm, jdm), &
      vlon(idm, jdm), vlat(idm, jdm), uvlon(idm, jdm), uvlat(idm, jdm), &
      tlon2(0:idm + 2, 0:jdm + 2), tlat2(0:idm + 2, 0:jdm + 2), &
      ulon2(0:idm + 2, 0:jdm + 2), ulat2(0:idm + 2, 0:jdm + 2), &
      vlon2(0:idm + 2, 0:jdm + 2), vlat2(0:idm + 2, 0:jdm + 2), &
      uvlon2(0:idm + 2, 0:jdm + 2), uvlat2(0:idm + 2, 0:jdm + 2), &
      tlon_crns(idm, jdm, ncrns), tlat_crns(idm, jdm, ncrns), &
      ulon_crns(idm, jdm, ncrns), ulat_crns(idm, jdm, ncrns), &
      vlon_crns(idm, jdm, ncrns), vlat_crns(idm, jdm, ncrns), &
      uvlon_crns(idm, jdm, ncrns), uvlat_crns(idm, jdm, ncrns), &
      tlon_crnsp(ncrns, idm, jdm), tlat_crnsp(ncrns, idm, jdm), &
      ulon_crnsp(ncrns, idm, jdm), ulat_crnsp(ncrns, idm, jdm), &
      vlon_crnsp(ncrns, idm, jdm), vlat_crnsp(ncrns, idm, jdm), &
      uvlon_crnsp(ncrns, idm, jdm), uvlat_crnsp(ncrns, idm, jdm), &
      xvec(idm), yvec(jdm), angle(idm, jdm), tarea(idm, jdm), &
      uarea(idm, jdm), stat=status)
    if (status /= 0) stop 'cannot ALLOCATE enough memory (1)'

    do i = 1, idm
      xvec(i) = i
    end do
    do j = 1, jdm
      yvec(j) = j
    end do

    status = nf90_inq_varid(ncid, 'TLON', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, tlon)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'TLAT', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, tlat)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'ULON', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, uvlon)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'ULAT', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, uvlat)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'ANGLE', rhid)

    do j = 1, jdm
      do i = 1, idm
        if (tlat(i, j) > 1e20) tlat(i, j) = 1e20
        if (tlon(i, j) > 1e20) tlon(i, j) = 1e20
        if (uvlat(i, j) > 1e20) uvlat(i, j) = 1e20
        if (uvlon(i, j) > 1e20) uvlon(i, j) = 1e20
      end do
    end do
    if (status == nf90_noerr) then
      status = nf90_get_var(ncid, rhid, angle)
      call handle_ncerror(status)
    else
      write(*, *) 'ANGLE not in sea ice file. Setting to zero.'
      angle = 0
    end if
    status = nf90_inq_varid(ncid, 'tarea', rhid)
    if (status == nf90_noerr) then
      status = nf90_get_var(ncid, rhid, tarea)
      call handle_ncerror(status)
    else
      write(*, *) 'tarea not in sea ice file. Setting to one.'
      tarea = 1
    end if
    status = nf90_inq_varid(ncid, 'uarea', rhid)
    if (status == nf90_noerr) then
      status = nf90_get_var(ncid, rhid, uarea)
      call handle_ncerror(status)
    else
      write(*, *) 'uarea not in sea ice file. Setting to one.'
      uarea = 1
    end if

    ! Read calendar info (override/change units)
    status = nf90_inq_varid(ncid, 'time', rhid)
    call handle_ncerror(status)
    status = nf90_get_att(ncid, rhid, 'calendar', calendar)
    call handle_ncerror(status)
    write(calunits(12:15), '(i4.4)') exprefyear

    ! Close file
    status = nf90_close(ncid)
    call handle_ncerror(status)

    ! Compute extended coordinate fields
    tlat2(1:idm, 1:jdm) = tlat
    tlon2(1:idm, 1:jdm) = tlon
    uvlat2(1:idm, 1:jdm) = uvlat
    uvlon2(1:idm, 1:jdm) = uvlon
    do i = 1, idm
      call sphextpnt(tlat2(i, 2), tlon2(i, 2), &
        tlat2(i, 1), tlon2(i, 1), tlat2(i, 0), tlon2(i, 0))
      call sphextpnt(tlat2(i, jdm - 1), tlon2(i, jdm - 1), &
        tlat2(i, jdm), tlon2(i, jdm), tlat2(i, jdm + 1), tlon2(i, jdm + 1))

      call sphextpnt(uvlat2(i, 2), uvlon2(i, 2), &
        uvlat2(i, 1), uvlon2(i, 1), uvlat2(i, 0), uvlon2(i, 0))
      call sphextpnt(uvlat2(i, jdm - 1), uvlon2(i, jdm - 1), &
        uvlat2(i, jdm), uvlon2(i, jdm), uvlat2(i, jdm + 1), uvlon2(i, jdm + 1))
      call sphextpnt(uvlat2(i, jdm), uvlon2(i, jdm), uvlat2(i, jdm + 1), &
        uvlon2(i, jdm + 1), uvlat2(i, jdm + 2), uvlon2(i, jdm + 2))
    end do
    do j = 0, jdm + 2
      call sphextpnt(tlat2(2, j), tlon2(2, j), &
        tlat2(1, j), tlon2(1, j), tlat2(0, j), tlon2(0, j))
      call sphextpnt(tlat2(idm - 1, j), tlon2(idm - 1, j), &
        tlat2(idm, j), tlon2(idm, j), tlat2(idm + 1, j), tlon2(idm + 1, j))

      call sphextpnt(uvlat2(2, j), uvlon2(2, j), &
        uvlat2(1, j), uvlon2(1, j), uvlat2(0, j), uvlon2(0, j))
      call sphextpnt(uvlat2(idm - 1, j), uvlon2(idm - 1, j), &
        uvlat2(idm, j), uvlon2(idm, j), uvlat2(idm + 1, j), uvlon2(idm + 1, j))
      call sphextpnt(uvlat2(idm, j), uvlon2(idm, j), uvlat2(idm + 1, j), &
        uvlon2(idm + 1, j), uvlat2(idm + 2, j), uvlon2(idm + 2, j))
    end do
    ulon2 = uvlon2
    ulat2 = uvlat2
    vlon2 = uvlon2
    vlat2 = uvlat2

    ! Interpolate u,v points
    do j = 0, jdm + 1
      do i = 0, idm + 1
        call sphmidpnt(ulat2(i, j), ulon2(i, j), ulat2(i + 1, j), &
          ulon2(i + 1, j), theta, lambda)
        ulat2(i, j) = theta
        ulon2(i, j) = lambda
        call sphmidpnt(vlat2(i, j), vlon2(i, j), vlat2(i, j + 1), &
          vlon2(i, j + 1), theta, lambda)
        vlat2(i, j) = theta
        vlon2(i, j) = lambda
      end do
    end do
    ulat = ulat2(1:idm, 1:jdm)
    ulon = ulon2(1:idm, 1:jdm)
    vlat = vlat2(1:idm, 1:jdm)
    vlon = vlon2(1:idm, 1:jdm)

    ! Compute corner points
    do j = 1, jdm
      do i = 1, idm
        call sphmidpnt(tlat2(i - 1, j), tlon2(i - 1, j), tlat2(i, j - 1), &
          tlon2(i, j - 1), tlat_crns(i, j, 1), tlon_crns(i, j, 1))
        call sphmidpnt(tlat2(i - 1, j), tlon2(i - 1, j), tlat2(i, j + 1), &
          tlon2(i, j + 1), tlat_crns(i, j, 2), tlon_crns(i, j, 2))
        call sphmidpnt(tlat2(i + 1, j), tlon2(i + 1, j), tlat2(i, j + 1), &
          tlon2(i, j + 1), tlat_crns(i, j, 3), tlon_crns(i, j, 3))
        call sphmidpnt(tlat2(i + 1, j), tlon2(i + 1, j), tlat2(i, j - 1), &
          tlon2(i, j - 1), tlat_crns(i, j, 4), tlon_crns(i, j, 4))

        call sphmidpnt(ulat2(i - 1, j), ulon2(i - 1, j), ulat2(i, j - 1), &
          ulon2(i, j - 1), ulat_crns(i, j, 1), ulon_crns(i, j, 1))
        call sphmidpnt(ulat2(i - 1, j), ulon2(i - 1, j), ulat2(i, j + 1), &
          ulon2(i, j + 1), ulat_crns(i, j, 2), ulon_crns(i, j, 2))
        call sphmidpnt(ulat2(i + 1, j), ulon2(i + 1, j), ulat2(i, j + 1), &
          ulon2(i, j + 1), ulat_crns(i, j, 3), ulon_crns(i, j, 3))
        call sphmidpnt(ulat2(i + 1, j), ulon2(i + 1, j), ulat2(i, j - 1), &
          ulon2(i, j - 1), ulat_crns(i, j, 4), ulon_crns(i, j, 4))

        call sphmidpnt(vlat2(i - 1, j), vlon2(i - 1, j), vlat2(i, j - 1), &
          vlon2(i, j - 1), vlat_crns(i, j, 1), vlon_crns(i, j, 1))
        call sphmidpnt(vlat2(i - 1, j), vlon2(i - 1, j), vlat2(i, j + 1), &
          vlon2(i, j + 1), vlat_crns(i, j, 2), vlon_crns(i, j, 2))
        call sphmidpnt(vlat2(i + 1, j), vlon2(i + 1, j), vlat2(i, j + 1), &
          vlon2(i, j + 1), vlat_crns(i, j, 3), vlon_crns(i, j, 3))
        call sphmidpnt(vlat2(i + 1, j), vlon2(i + 1, j), vlat2(i, j - 1), &
          vlon2(i, j - 1), vlat_crns(i, j, 4), vlon_crns(i, j, 4))

        call sphmidpnt(uvlat2(i - 1, j), uvlon2(i - 1, j), uvlat2(i, j - 1), &
          uvlon2(i, j - 1), uvlat_crns(i, j, 1), uvlon_crns(i, j, 1))
        call sphmidpnt(uvlat2(i - 1, j), uvlon2(i - 1, j), uvlat2(i, j + 1), &
          uvlon2(i, j + 1), uvlat_crns(i, j, 2), uvlon_crns(i, j, 2))
        call sphmidpnt(uvlat2(i + 1, j), uvlon2(i + 1, j), uvlat2(i, j + 1), &
          uvlon2(i, j + 1), uvlat_crns(i, j, 3), uvlon_crns(i, j, 3))
        call sphmidpnt(uvlat2(i + 1, j), uvlon2(i + 1, j), uvlat2(i, j - 1), &
          uvlon2(i, j - 1), uvlat_crns(i, j, 4), uvlon_crns(i, j, 4))
      end do
    end do

    ! Permute to compensate for dimension bug in CMOR
    do n = 1, ncrns
      do j = 1, jdm
        do i = 1, idm
          tlon_crnsp(n, i, j) = tlon_crns(i, j, n)
          tlat_crnsp(n, i, j) = tlat_crns(i, j, n)
          ulon_crnsp(n, i, j) = ulon_crns(i, j, n)
          ulat_crnsp(n, i, j) = ulat_crns(i, j, n)
          vlon_crnsp(n, i, j) = vlon_crns(i, j, n)
          vlat_crnsp(n, i, j) = vlat_crns(i, j, n)
          uvlon_crnsp(n, i, j) = uvlon_crns(i, j, n)
          uvlat_crnsp(n, i, j) = uvlat_crns(i, j, n)
          if (tlon_crnsp(n, i, j) < 0.) &
            tlon_crnsp(n, i, j) = tlon_crnsp(n, i, j) + 360
          if (ulon_crnsp(n, i, j) < 0.) &
            ulon_crnsp(n, i, j) = ulon_crnsp(n, i, j) + 360
          if (vlon_crnsp(n, i, j) < 0.) &
            vlon_crnsp(n, i, j) = vlon_crnsp(n, i, j) + 360
          if (uvlon_crnsp(n, i, j) < 0.) &
            uvlon_crnsp(n, i, j) = uvlon_crnsp(n, i, j) + 360
          if (ulon(i, j) < 0.) ulon(i, j) = ulon(i, j) + 360
          if (vlon(i, j) < 0.) vlon(i, j) = vlon(i, j) + 360
        end do
      end do
    end do

  end subroutine read_gridinfo_ifile

  ! -----------------------------------------------------------------

  subroutine open_ofile(fx)

    implicit none

    logical, optional, intent(in)   :: fx
    logical                         :: fxflag

    real :: fac1, fac2
    integer, parameter :: ndimmax = 10
    integer :: n, ndims, dimids(ndimmax), dimlens(ndimmax)
    integer :: physics_version = 1, initialization_method = 1
    character(len=slenmax) :: coord, ivnm1a, ivnm2a, ivnm1b, ivnm2b
    real(kind=8), dimension(:), allocatable :: tmp1d, tmp1d_2
    real(kind=8), dimension(:, :), allocatable :: tmp2d

    ! Check if output variable should have time coordinate
    fxflag = .false.
    if (present(fx)) then
      if (fx) fxflag = .true.
    end if

    ! Inquire variable units and dimensions in input file
    status = nf90_open(fnm, nf90_nowrite, ncid)
    call handle_ncerror(status)

    call resolve_vnm2(slenmax, ivnm, ivnm1a, ivnm2a, ivnm1b, ivnm2b, fac1, fac2)
    if (verbose) then
      if (index(special, 'Dfield2') > 0) then
        write(*, *) 'Resolve variable term: ', trim(ivnm1a), '/', &
          trim(ivnm1b), '*', fac1, '+', trim(ivnm2a), '/', trim(ivnm2b), &
          '*', fac2
      else
        write(*, *) 'Resolve variable term: ', trim(ivnm1a), '*', &
          trim(ivnm1b), '*', fac1, '+', trim(ivnm2a), '*', trim(ivnm2b), &
          '*', fac2
      end if
      status = nf90_inq_varid(ncid, trim(ivnm1a), rhid)
    end if
    if (status /= nf90_noerr) then
      write(*, *) 'cannot find input variable ', trim(ivnm1a)
      stop
    end if
    status = nf90_inquire_variable(ncid, rhid, ndims=ndims)
    call handle_ncerror(status)
    status = nf90_inquire_variable(ncid, rhid, dimids=dimids(1:ndims))
    call handle_ncerror(status)
    if (ndims < 3) then
      write(*, *) 'Variable ', trim(ivnm1a), ' has too few dimensions'
    end if
    dimlens = 1
    do n = 1, ndims
      status = nf90_inquire_dimension(ncid, dimids(n), len=dimlens(n))
      call handle_ncerror(status)
    end do
    if (dimlens(1) /= idm) then
      write(*, *) 'unexpected first dimension of variable ', &
        trim(ivnm1a), ': ', dimlens(1), ' versus idm=', idm
      stop
    end if
    if (dimlens(2) /= jdm) then
      write(*, *) 'unexpected second dimension of variable ', &
        trim(ivnm1a), ': ', dimlens(2), ' versus jdm=', idm
      stop
    end if
    if (allocated(fld)) deallocate(fld, fld2, fld3, fldacc)
    allocate(fld(idm, jdm), fld2(idm, jdm), fld3(idm, jdm), &
      fldacc(idm, jdm), stat=status)
    if (status /= 0) stop 'cannot ALLOCATE enough memory (4)'
    if (allocated(fld3d)) deallocate(fld3d)
    allocate(fld3d(idm, jdm, 5), stat=status)
    if (status /= 0) stop 'cannot ALLOCATE enough memory (5)'

    if (len_trim(vunits) == 0) then
      status = nf90_get_att(ncid, rhid, 'units', vunits)
      call handle_ncerror(status)
      if (trim(vunits) == 'mm/s') vunits = 'kg m-2 s-1'
    end if

    coord = ' '
    status = nf90_get_att(ncid, rhid, 'coordinates', coord)

    status = nf90_close(ncid)
    call handle_ncerror(status)

    ! Inquire time dimension of output variable
    write(*, *) 'get tcoord'
    if (.not. fxflag) then
      call get_timecoord(trim(tabledir)//trim(table), ovnm, tcoord)
    end if
    write(*, *) 'get tcoord (after)'

    ! Call CMOR setup
    if (verbose) then
      if (createsubdirs) then
        error_flag = cmor_setup(inpath=trim(ibasedir), &
          netcdf_file_action=CMOR_REPLACE_4, set_verbosity=CMOR_NORMAL, &
          create_subdirectories=1)
      else
        error_flag = cmor_setup(inpath=trim(ibasedir), &
          netcdf_file_action=CMOR_REPLACE_4, set_verbosity=CMOR_NORMAL, &
          create_subdirectories=0)
      end if
    else
      if (createsubdirs) then
        error_flag = cmor_setup(inpath=trim(ibasedir), &
          netcdf_file_action=CMOR_REPLACE_4, set_verbosity=CMOR_QUIET, &
          create_subdirectories=1)
      else
        error_flag = cmor_setup(inpath=trim(ibasedir), &
          netcdf_file_action=CMOR_REPLACE_4, set_verbosity=CMOR_QUIET, &
          create_subdirectories=0)
      end if
    end if
    if (error_flag /= 0) stop 'Problem setting up CMOR'

    ! Derive physics_version and initialization_method from
    ! parent_experiment_rip
    if (trim(parent_experiment_rip) /= 'r1i1p1' .and. &
        trim(parent_experiment_rip) /= 'N/A' .and. &
        trim(parent_experiment_rip) /= 'no parent') then
      read(parent_experiment_rip(index(parent_experiment_rip, 'i') + 1: &
        index(parent_experiment_rip, 'p') - 1), *) initialization_method
      read(parent_experiment_rip(index(parent_experiment_rip, 'p') + 1:), &
        *) physics_version
    end if

    ! Define output dataset
    call write_namelist_json(icegrid, icegrid_label, icegrid_resolution, ovnm)
    error_flag = cmor_dataset_json(namelist_file_json)
    call system('rm '//trim(namelist_file_json))

    ! Define horizontal axes
    iaxid = cmor_axis( &
      table=trim(tabledir)//trim(tgrids), &
      table_entry='i_index', &
      units='1', &
      length=idm, &
      coord_vals=xvec)
    jaxid = cmor_axis( &
      table=trim(tabledir)//trim(tgrids), &
      table_entry='j_index', &
      units='1', &
      length=jdm, &
      coord_vals=yvec)

    if (coord(1:1) == 'T' .and. ovnm(1:6) /= 'sidmasstran') then
      grdid = cmor_grid( &
        axis_ids=(/iaxid, jaxid/), &
        latitude=tlat, &
        longitude=tlon, &
        latitude_vertices=tlat_crnsp, &
        longitude_vertices=tlon_crnsp)
    else if (trim(ovnm) == 'sidmasstranx') then
      grdid = cmor_grid( &
        axis_ids=(/iaxid, jaxid/), &
        latitude=ulat, &
        longitude=ulon, &
        latitude_vertices=ulat_crnsp, &
        longitude_vertices=ulon_crnsp)
    else if (trim(ovnm) == 'sidmasstrany') then
      grdid = cmor_grid( &
        axis_ids=(/iaxid, jaxid/), &
        latitude=vlat, &
        longitude=vlon, &
        latitude_vertices=vlat_crnsp, &
        longitude_vertices=vlon_crnsp)
    else
      grdid = cmor_grid( &
        axis_ids=(/iaxid, jaxid/), &
        latitude=uvlat, &
        longitude=uvlon, &
        latitude_vertices=uvlat_crnsp, &
        longitude_vertices=uvlon_crnsp)
    end if

    ! Define vertical dummy coordinate
    if (trim(zcoord) == 'olevel') then
      allocate(tmp1d(1), tmp2d(2, 1))
      tmp1d(:) = (/5.d0/)
      tmp2d(:, 1) = (/0.d0, 10.d0/)
      kaxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry='depth_coord', &
        units='m', &
        length=1, &
        coord_vals=tmp1d, &
        cell_bounds=tmp2d)
      deallocate(tmp1d, tmp2d)
    else if (trim(zcoord) == 'iceband') then
      allocate(tmp1d(5), tmp2d(2, 5))
      tmp1d(:) = (/0.6445072d0, 1.391433d0, 2.470179d0, &
                   4.567288d0, 1.d8/)
      tmp2d(1, :) = (/0.d0, 0.6445072d0, 1.391433d0, &
                      2.470179d0, 4.567288d0/)
      tmp2d(2, :) = (/0.6445072d0, 1.391433d0, 2.470179d0, &
                      4.567288d0, 1.d8/)
      kaxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry='iceband', &
        units='m', &
        length=5, &
        coord_vals=tmp1d, &
        cell_bounds=tmp2d)
      deallocate(tmp1d, tmp2d)
    end if

    ! Define time axis
    if (.not. fxflag) then
      taxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry=trim(tcoord), &
        units=trim(calunits), &
        length=1)
    end if

    ! Define output variable
    if (fxflag) then
      varid = cmor_variable( &
        table=trim(tabledir)//trim(table), &
        table_entry=trim(ovnm), &
        units=trim(vunits), &
        axis_ids=(/grdid/), &
        missing_value=1e20, &
        original_name=trim(ivnm), &
        comment=trim(vcomment))
    else
      if (trim(ovnm) == 'transifs') then
        varid = cmor_variable( &
          table=trim(tabledir)//trim(table), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/taxid/), &
          missing_value=1e20, &
          comment=trim(vcomment))
      else
        if (trim(zcoord) == 'olevel' .or. trim(zcoord) == 'iceband') then
          varid = cmor_variable( &
            table=trim(tabledir)//trim(table), &
            table_entry=trim(ovnm), &
            units=trim(vunits), &
            axis_ids=(/grdid, kaxid, taxid/), &
            original_name=trim(ivnm), &
            missing_value=1e20, &
            positive=trim(vpositive), &
            comment=trim(vcomment))
        else if (index(special, 'glbsum') > 0) then
          varid = cmor_variable( &
            table=trim(tabledir)//trim(table), &
            table_entry=trim(ovnm), &
            units=trim(vunits), &
            axis_ids=(/taxid/), &
            original_name=trim(ivnm), &
            missing_value=1e20, &
            positive=trim(vpositive), &
            comment=trim(vcomment))
        else
          varid = cmor_variable( &
            table=trim(tabledir)//trim(table), &
            table_entry=trim(ovnm), &
            units=trim(vunits), &
            axis_ids=(/grdid, taxid/), &
            original_name=trim(ivnm), &
            missing_value=1e20, &
            positive=trim(vpositive), &
            comment=trim(vcomment))
        end if
      end if
    end if
#ifdef DEFLATE
    error_flag = cmor_set_deflate(varid, 1, 1, 5)
#endif

  end subroutine open_ofile

  ! -----------------------------------------------------------------

  subroutine close_ofile

    implicit none

    status = cmor_close()
    if (status /= 0) stop 'problem closing CMOR output file'

  end subroutine close_ofile

  ! -----------------------------------------------------------------

  subroutine read_field

    implicit none

    real :: fac1, fac2
    integer :: ind
    character(len=slenmax) :: ivnm1a, ivnm2a, ivnm1b, ivnm2b

    ! Open input file
    status = nf90_open(fnm, nf90_nowrite, ncid)
    call handle_ncerror(status)

    ! Read data
    call resolve_vnm2(slenmax, ivnm, ivnm1a, ivnm2a, ivnm1b, ivnm2b, fac1, fac2)
    status = nf90_inq_varid(ncid, trim(ivnm1a), rhid)
    if (status /= nf90_noerr) then
      write(*, *) 'cannot find input variable ', trim(ivnm1a)
      stop
    end if
    status = nf90_get_var(ncid, rhid, fld)
    call handle_ncerror(status)
    if (fac1 /= 1) then
      fld = fld * fac1
    end if

    if (len_trim(ivnm1b) > 0) then
      status = nf90_inq_varid(ncid, trim(ivnm1b), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm1b)
        stop
      end if
      status = nf90_get_var(ncid, rhid, fld2)
      call handle_ncerror(status)
      if (index(special, 'Dfield2') > 0) then
        fld = fld / fld2
      else
        fld = fld * fld2
      end if
    end if

    if (len_trim(ivnm2a) > 0) then
      status = nf90_inq_varid(ncid, trim(ivnm2a), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm2a)
        stop
      end if
      status = nf90_get_var(ncid, rhid, fld2)
      call handle_ncerror(status)
      if (len_trim(ivnm2b) > 0) then
        status = nf90_inq_varid(ncid, trim(ivnm2b), rhid)
        if (status /= nf90_noerr) then
          write(*, *) 'cannot find input variable ', trim(ivnm2b)
          stop
        end if
        status = nf90_get_var(ncid, rhid, fld3)
        call handle_ncerror(status)
        if (index(special, 'Dfield2') > 0) then
          fld2 = fld2 / fld3
        else
          fld2 = fld2 * fld3
        end if
      end if
      fld = fld + fld2 * fac2
    end if

    status = nf90_close(ncid)
    call handle_ncerror(status)

  end subroutine read_field

  ! -----------------------------------------------------------------

  subroutine read_tslice(rec, badrec, fname)

    implicit none

    real                                    :: fac1, fac2
    integer, intent(in)                     :: rec
    logical, intent(out)                    :: badrec
    character(len=*), intent(in), optional  :: fname
    integer, save                           :: fid
    integer                                 :: i, j, i1, j1
    integer                                 :: nc
    character(len=slenmax)                  :: ivnm1a, ivnm2a, ivnm1b, ivnm2b

    ! Open input file
    if (present(fname)) then
      status = nf90_open(fname, nf90_nowrite, fid)
      call handle_ncerror(status)
    else
      status = nf90_open(fnm, nf90_nowrite, fid)
      call handle_ncerror(status)
    end if

    if (.false.) then
      ! Read time information
      badrec = .false.
      status = nf90_inq_varid(fid, 'time', rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find time variable'
        stop
      end if
      status = nf90_get_var(fid, rhid, tval, (/rec/), (/1/))
      call handle_ncerror(status)
      status = nf90_inq_varid(fid, 'time_bounds', rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find time_bounds variable'
        stop
      end if
      status = nf90_get_var(fid, rhid, tbnds, (/1, rec/), (/2, 1/))
      call handle_ncerror(status)
      if (linstant) then
        ! Exception for instantaneous 6+3 hourly data
        if (rec == 1) then
          status = nf90_inq_varid(fid, 'time', rhid)
          call handle_ncerror(status)
          status = nf90_get_var(fid, rhid, tval2, (/2/), (/2/))
          call handle_ncerror(status)
          if (tval(1) == tval2(1)) then
            tbnds(2, 1) = tval(1) + tval2(1) - tval2(2)
            badrec = .true.
          end if
        end if
        tbnds(1, 1) = tbnds(2, 1)
      end if
      tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

      ! correct erroneous intial time bound
      tbnds(1, 1) = max(0., tbnds(1, 1))
      tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))
    end if

    ! Read data
    call resolve_vnm2(slenmax, ivnm, ivnm1a, ivnm2a, ivnm1b, ivnm2b, fac1, fac2)
    if (trim(ovnm) == 'transifs') then
      status = nf90_inq_varid(fid, 'transix', rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable transix'
        stop
      end if
      status = nf90_get_var(fid, rhid, fld, (/1, 1, rec/), (/idm, jdm, 1/))
      call handle_ncerror(status)
      status = nf90_inq_varid(fid, 'transiy', rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable transiy'
        stop
      end if
      status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idm, jdm, 1/))
      call handle_ncerror(status)
      ! Compute Fram Strait transport
      fld(1, 1) = transifs(seclen, iind, jind, iflg, jflg, fld, fld2)
      write(*, *) 'transifs=', fld(1, 1)
      vunits = 'kg/s'
    else
      status = nf90_inq_varid(fid, trim(ivnm1a), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm1a), status
        stop
      end if
      if (trim(zcoord) == 'iceband') then
        status = nf90_inq_dimid(fid, 'nc', dimid)
        call handle_ncerror(status)
        status = nf90_inquire_dimension(fid, dimid, len=nc)
        call handle_ncerror(status)
        status = nf90_get_var(fid, rhid, fld3d, (/1, 1, 1, rec/), &
          (/idm, jdm, nc, 1/))
      else
        status = nf90_get_var(fid, rhid, fld, (/1, 1, rec/), (/idm, jdm, 1/))
      end if
      call handle_ncerror(status)
    end if

    ! Rotate to east/north alignment if variable is a velocity (BYPASSED)
    if (.false.) then
      if (ivnm1a(1:4) == 'uvel' .or. ivnm1a(1:4) == 'vvel') then
        if (ivnm1a(1:4) == 'uvel') then
          status = nf90_inq_varid(fid, 'v'//trim(ivnm1a(2:)), rhid)
          if (status /= nf90_noerr) then
            write(*, *) 'cannot find input variable ', 'v'//trim(ivnm1a(2:))
            stop
          end if
        else
          status = nf90_inq_varid(fid, 'u'//trim(ivnm1a(2:)), rhid)
          if (status /= nf90_noerr) then
            write(*, *) 'cannot find input variable ', 'u'//trim(ivnm1a(2:))
            stop
          end if
        end if

        status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idm, jdm, 1/))
        call handle_ncerror(status)

        if (ivnm1a(1:4) == 'uvel') then
          call rotate_uv(idm, jdm, angle, fld, fld2)
        else
          call rotate_uv(idm, jdm, angle, fld2, fld)
        end if
      end if
    end if

    ! Apply user defined factors and linear combinations
    if (fac1 /= 1) then
      if (trim(zcoord) == 'iceband') then
        fld3d = fld3d * fac1
      else
        fld = fld * fac1
      end if
    end if

    if (len_trim(ivnm1b) > 0) then
      status = nf90_inq_varid(fid, trim(ivnm1b), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm1b)
        stop
      end if
      status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idm, jdm, 1/))
      call handle_ncerror(status)
      fld = fld * fld2
    end if

    if (len_trim(ivnm2a) > 0) then
      status = nf90_inq_varid(fid, trim(ivnm2a), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm2a)
        stop
      end if
      status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idm, jdm, 1/))
      call handle_ncerror(status)
      if (len_trim(ivnm2b) > 0) then
        status = nf90_inq_varid(fid, trim(ivnm2b), rhid)
        if (status /= nf90_noerr) then
          write(*, *) 'cannot find input variable ', trim(ivnm2b)
          stop
        end if
        status = nf90_get_var(fid, rhid, fld3, (/1, 1, rec/), (/idm, jdm, 1/))
        call handle_ncerror(status)
        fld2 = fld2 * fld3
      end if
      fld = fld + fld2 * fac2
    end if

    ! Do sea ice fraction weighting if required
    if (index(special, 'Xaiu-1') > 0) then
      ivnm1a = 'aice'
      status = nf90_inq_varid(fid, trim(ivnm1a), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm1a)
        stop
      end if
      status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idm, jdm, 1/))
      call handle_ncerror(status)
      do j = 1, jdm
        do i = 1, idm
          if (fld2(i, j) > 1e20) fld2(i, j) = 0
        end do
      end do
      do j = 1, jdm
        j1 = min(j + 1, jdm)
        do i = 1, idm
          i1 = mod(i, idm) + 1
          fac1 = 0.01 * 0.25 * (tarea(i, j) * fld2(i, j) + tarea(i1, j) * fld2(i1, j) + &
            tarea(i, j1) * fld2(i, j1) + tarea(i1, j1) * fld2(i1, j1)) / uarea(i, j)
          if (fac1 > 0.001 .and. fld(i, j) < 1e20) then
            fld(i, j) = fld(i, j) / fac1
          else
            fld(i, j) = 1e20
          end if
        end do
      end do
    end if

    status = nf90_close(fid)
    call handle_ncerror(status)

  end subroutine read_tslice

  ! -----------------------------------------------------------------

  subroutine write_field

    implicit none

    integer :: i, j

    ! Set zero on ocean grid cells
    do j = 1, jdm
      do i = 1, idm
        if (abs(fld(i, j)) > 2e20) fld(i, j) = 0.
      end do
    end do

    ! Store variable
    error_flag = cmor_write( &
      var_id=varid, &
      data=fld)

  end subroutine write_field

  ! -----------------------------------------------------------------

  subroutine write_tslice

    implicit none

    integer :: i, j, k

    ! Set missing on land grid cells
    do j = 1, jdm
      do i = 1, idm
        if (abs(fld(i, j)) > 1e20) fld(i, j) = 1e20
      end do
    end do

    ! Set missing on land grid cells
    do k = 1, 5
      do j = 1, jdm
        do i = 1, idm
          if (abs(fld3d(i, j, k)) > 1e20) fld3d(i, j, k) = 1e20
        end do
      end do
    end do

    ! Store variable
    if (trim(ovnm) == 'transifs') then
      error_flag = cmor_write( &
        var_id=varid, &
        data=reshape(fld(1:1, 1:1), (/1/)), &
        ntimes_passed=1, &
        time_vals=tval, &
        time_bnds=tbnds)
    else if (trim(tcoord) == 'time1') then
      error_flag = cmor_write( &
        var_id=varid, &
        data=fld, &
        ntimes_passed=1, &
        time_vals=tval)
    else
      if (trim(zcoord) == 'olevel') then
        error_flag = cmor_write( &
          var_id=varid, &
          data=reshape(fld, (/idm, jdm, 1/)), &
          ntimes_passed=1, &
          time_vals=tval, &
          time_bnds=tbnds)
      else if (trim(zcoord) == 'iceband') then
        error_flag = cmor_write( &
          var_id=varid, &
          data=reshape(fld3d, (/idm, jdm, 5/)), &
          ntimes_passed=1, &
          time_vals=tval, &
          time_bnds=tbnds)
      else if (index(special, 'glbsum') > 0) then
        error_flag = cmor_write( &
          var_id=varid, &
          data=reshape(fld(1:1, 1:1), (/1/)), &
          ntimes_passed=1, &
          time_vals=tval, &
          time_bnds=tbnds)
      else
        error_flag = cmor_write( &
          var_id=varid, &
          data=fld, &
          ntimes_passed=1, &
          time_vals=tval, &
          time_bnds=tbnds)
      end if
    end if

  end subroutine write_tslice

end module m_modelsice
