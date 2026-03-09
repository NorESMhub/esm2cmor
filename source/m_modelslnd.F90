module m_modelslnd

  use netcdf
  use cmor_users_functions
  use m_utilities
  use m_namelists

  implicit none

  ! Netcdf variables
  integer               :: ncid, rhid, dimid, status

  ! Grid dimensions and variables
  integer, save         :: ii, jj, idm, jdm, kdm, idmrof, jdmrof, ldm
  real(kind=8), allocatable, save, dimension(:) :: lon, lat, lev, lonrof, latrof
  real(kind=8), allocatable, save, dimension(:, :)  :: lon_bnds, lat_bnds,&
    lev_bnds, lonrof_bnds, latrof_bnds
  character(len=slenmax), allocatable, save, dimension(:)   :: ltype
  character(len=slenmax), save                              :: zcoord, tcoord

  ! Dataset related variables
  character(len=slenmax), save  :: ivnm, ovnm, vunits, vpositive

  ! Table related variables
  character(len=lenmax)         :: table

  ! Cmor parameters
  integer, save :: iaxid, jaxid, kaxid, taxid, varid, table_id, error_flag

  ! Data fields
  real(kind=8), allocatable, save, dimension(:, :, :)   :: fld, fld2, fldacc

  ! Auxillary variables for special operations
  character(len=slenmax), save                          :: special, str1, str2

contains

  ! -----------------------------------------------------------------

  subroutine lnd2cmor

    implicit none

    logical :: badrec, last
    integer :: m, n, nrec

    badrec = .false.

    ! Print start information
    if (verbose) then
      write(*, *)
      write(*, *) '---------------------------'
      write(*, *) '--- Process land output ---'
      write(*, *) '---------------------------'
      write(*, *)
    end if

    ! Read grid information from input files
    write(*, *) 'Read grid information from input files'
    itag = taglmon
    call scan_files(reset=.true.)
    if (len_trim(fnm) == 0) return
    call read_gridinfo_ifile

    ! Process table fx
    fnm = pfx
    table = tfx
    linstant = .false.
    do n = 1, nfx
      if (skip_variable(n, nfx, dfx)) cycle

      ! Map namelist variables
      ovnm = vfx(ovnmpos, n)
      ivnm = vfx(ivnmpos, n)
      special = vfx(3, n)
      vunits = ' '
      vpositive = ' '

      ! Skip variable?
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (trim(zcoord) == 'sdepth') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Check if input variable is present
      if (len_trim(pfx) == 0) cycle
      if (.not. var_in_file(fnm, ivnm)) cycle

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

    ! Process table Eyr
    write(*, *) 'Process table Eyr'
    fnm = pEyr
    table = tEyr
    linstant = .false.
    do n = 1, nEyr
      if (skip_variable(n, nEyr, dEyr)) cycle

      ! Map namelist variables
      ovnm = vEyr(ovnmpos, n)
      ivnm = vEyr(ivnmpos, n)
      special = vEyr(3, n)
      vunits = ' '
      vpositive = ' '

      ! Skip variable?
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)

      ! Choose history file
      itag = taglyr

      ! Check if input variable is present
      if (len_trim(pEyr) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over input files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rEyr) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        nrec = 0
        fldacc = 0.
        last = .false.
        do
          if (len_trim(pEyr) == 0) call scan_files(reset=.false.)
          if (rec == 0) then
            last = .true.
            exit
          end if
          nrec = nrec + 1
          call read_tslice(rec, badrec, fnm)
          fldacc = fldacc + fld
          write(*, *) "WARNING: NO TIME_BOUNDS CHECK FOR Eyr VARIABLE!"
          exit
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
        if (mod(m, rEyr) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rEyr) > 0) call close_ofile

    end do

    ! Process table Lmon
    write(*, *) 'Process table Lmon'
    fnm = plmon
    table = tlmon
    linstant = .false.
    do n = 1, nlmon
      if (skip_variable(n, nlmon, dlmon)) cycle

      ! Map namelist variables
      ovnm = vlmon(ovnmpos, n)
      ivnm = vlmon(ivnmpos, n)
      special = vlmon(3, n)
      vunits = ' '
      vpositive = ' '

      ! Skip variable?
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (trim(zcoord) == 'sdepth') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      select case (trim(special))
      case ('day2mon')
        itag = taglday
      case ('3hr2mon')
        itag = tagl3hr
      case default
        itag = taglmon
      end select

      ! Check if input variable is present
      if (len_trim(plmon) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over input files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rlmon) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        nrec = 0
        fldacc = 0.
        last = .false.
        do
          if (len_trim(plmon) == 0) call scan_files(reset=.false.)
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
        if (mod(m, rlmon) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rlmon) > 0) call close_ofile

    end do

    ! Process table Emon
    write(*, *) 'Process table Emon'
    fnm = pEmon
    table = tEmon
    linstant = .false.
    do n = 1, nEmon
      if (skip_variable(n, nEmon, dEmon)) cycle

      ! Map namelist variables
      ovnm = vEmon(ovnmpos, n)
      ivnm = vEmon(ivnmpos, n)
      special = vEmon(3, n)
      vunits = ' '
      vpositive = ' '

      ! Skip variable?
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (trim(zcoord) == 'sdepth') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      select case (trim(special))
      case ('day2mon')
        itag = taglday
      case ('3hr2mon')
        itag = tagl3hr
      case default
        itag = taglmon
      end select

      ! Check if input variable is present
      if (len_trim(pEmon) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over input files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rEmon) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        nrec = 0
        fldacc = 0.
        last = .false.
        do
          if (len_trim(pEmon) == 0) call scan_files(reset=.false.)
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
        if (mod(m, rEmon) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rEmon) > 0) call close_ofile

    end do

    ! Process table limon
    write(*, *) 'Process table limon'
    fnm = plimon
    table = tlimon
    linstant = .false.
    do n = 1, nlimon
      if (skip_variable(n, nlimon, dlimon)) cycle

      ! Map namelist variables
      ovnm = vlimon(ovnmpos, n)
      ivnm = vlimon(ivnmpos, n)
      special = vlimon(3, n)
      vunits = ' '
      vpositive = ' '

      ! Skip variable?
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (trim(zcoord) == 'sdepth') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      select case (trim(special))
      case ('day2mon')
        itag = taglday
      case ('3hr2mon')
        itag = tagl3hr
      case default
        itag = taglmon
      end select

      ! Check if input variable is present
      if (len_trim(plimon) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over input files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rlimon) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        nrec = 0
        fldacc = 0.
        last = .false.
        do
          if (len_trim(plimon) == 0) call scan_files(reset=.false.)
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
        if (mod(m, rlimon) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rlimon) > 0) call close_ofile

    end do

    ! Process table day
    write(*, *) 'Process table day'
    fnm = pday
    table = tday
    linstant = .false.
    do n = 1, nday
      if (skip_variable(n, nday, dday)) cycle

      ! Map namelist variables
      ovnm = vday(ovnmpos, n)
      ivnm = vday(ivnmpos, n)
      special = vday(3, n)
      vunits = ' '
      vpositive = ' '

      ! Skip variable?
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (trim(zcoord) == 'sdepth') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = taglday

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

    ! Process table Eday
    write(*, *) 'Process table Eday'
    fnm = pEday
    table = tEday
    linstant = .false.
    do n = 1, nEday
      if (skip_variable(n, nEday, dEday)) cycle

      ! Map namelist variables
      ovnm = vEday(ovnmpos, n)
      ivnm = vEday(ivnmpos, n)
      special = vEday(3, n)
      vunits = ' '
      vpositive = ' '

      ! Skip variable?
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (trim(zcoord) == 'sdepth') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = taglday

      ! Check if input variable is present
      if (len_trim(pEday) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over input files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rEday) == 0) call open_ofile

        ! Read variable into buffer
        rec = 0
        if (len_trim(pEday) == 0) call scan_files(reset=.false.)
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
        if (mod(m, rEday) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rEday) > 0) call close_ofile

    end do

    ! Process table 3hr (averaged fields)
    write(*, *) 'Process table 3hr (averaged fields)'
    fnm = p3hr
    table = t3hr
    linstant = .false.
    do n = 1, n3hr
      if (skip_variable(n, n3hr, d3hr)) cycle

      ! Map namelist variables
      ovnm = v3hr(ovnmpos, n)
      ivnm = v3hr(ivnmpos, n)
      special = v3hr(3, n)
      vunits = ' '
      vpositive = ' '

      ! Skip variable?
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (trim(zcoord) == 'sdepth') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = tagl3hr

      ! Check if input variable is present
      if (len_trim(p3hr) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over input files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, r3hr) == 0) call open_ofile

        ! Read variable into buffer
        rec = 0
        if (len_trim(p3hr) == 0) call scan_files(reset=.false.)
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
        if (mod(m, r3hr) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, r3hr) > 0) call close_ofile

    end do

    ! Process table 3hr (instantaneous fields)
    write(*, *) 'Process table 3hr (instantaneous fields)'
    fnm = p3hri
    table = t3hri
    linstant = .true.
    do n = 1, n3hri
      if (skip_variable(n, n3hri, d3hri)) cycle

      ! Map namelist variables
      ovnm = v3hri(ovnmpos, n)
      ivnm = v3hri(ivnmpos, n)
      special = v3hri(3, n)
      vunits = ' '
      vpositive = ' '

      ! Skip variable?
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (trim(zcoord) == 'sdepth') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = tagl3hri

      ! Check if input variable is present
      if (len_trim(p3hri) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over input files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, r3hri) == 0) call open_ofile

        ! Read variable into buffer
        rec = 0
        if (len_trim(p3hri) == 0) call scan_files(reset=.false.)
        if (rec == 0) exit
        tbnds(:, 1) = tbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        call special_post
        if (badrec) fld = 1e20

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, r3hri) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, r3hri) > 0) call close_ofile

    end do

    ! Process table E3hr (averaged fields)
    write(*, *) 'Process table E3hr (averaged fields)'
    fnm = pE3hr
    table = tE3hr
    linstant = .false.
    do n = 1, nE3hr
      if (skip_variable(n, nE3hr, dE3hr)) cycle

      ! Map namelist variables
      ovnm = vE3hr(ovnmpos, n)
      ivnm = vE3hr(ivnmpos, n)
      special = vE3hr(3, n)
      vunits = ' '
      vpositive = ' '

      ! Skip variable?
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (trim(zcoord) == 'sdepth') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = tagl3hr

      ! Check if input variable is present
      if (len_trim(pE3hr) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over input files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rE3hr) == 0) call open_ofile

        ! Read variable into buffer
        rec = 0
        if (len_trim(pE3hr) == 0) call scan_files(reset=.false.)
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
        if (mod(m, rE3hr) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rE3hr) > 0) call close_ofile

    end do

  end subroutine lnd2cmor

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

        ! Unit transformation: g m-2 -> kg m-2
      case ('kg m-2')
        vunits = 'kg m-2'

        ! Unit transformation: mm s-1 -> kg m-2 s-1
      case ('kg m-2 s-1')
        vunits = 'kg m-2 s-1'

        ! Fix micrometers units
      case ('micrometer')
        vunits = 'micrometers'

        ! Fix m-2 units
      case ('m-2')
        vunits = 'm-2'

        ! Convert units from radians2 to m2
      case ('rad2m')
        vunits = 'm2'

        ! Set positive attribute
      case ('positiveup')
        vpositive = 'up'
      case ('positivedo')
        vpositive = 'down'

      end select
      if (str1 == str2) exit
    end do

  end subroutine special_pre

  ! -----------------------------------------------------------------

  subroutine special_post

    implicit none

    integer                 :: i, j, k
    integer, dimension(12)  :: ndays

    data ndays /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

    str2 = special
    do
      if (index(str2, ';') > 0) then
        str1 = str2(1:index(str2, ';') - 1)
        str2 = str2(index(str2, ';') + 1:)
      else
        str1 = str2
      end if
      select case (str1)

        ! Fraction to percent
      case ('percent')
        fld = fld * 1e2

        ! Convert units from radians2 to m2
      case ('rad2m')
        fld = fld * 6.37122e6**2

        ! Unit transformation: g m-2 -> kg m-2
      case ('kg m-2')
        fld = fld * 1e-3

        ! Unit transformation: g m-2 s-1 -> kg m-2 s-1
      case ('kg m-2 s-1')
        fld = fld * 1e-3

        ! Unit transformation: mol to kg CH4
      case ('mol to kg CH4')
        fld = fld * 16.04246 * 1e-3

        ! Limit solid soil moisture (mask ice sheets)
      case ('limitmoist')
        do k = 1, kdm
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) < 1e20 .and. fld(i, j, k) > 5000) &
                fld(i, j, k) = 5000
            end do
          end do
        end do

        ! Set ocean points to missing value
      case ('missingval')
        do k = 1, kdm
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) > 1e20) fld(i, j, k) = 1e20
            end do
          end do
        end do

        ! Set ocean points to zero
      case ('miss2zero')
        do k = 1, kdm
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) > 1e20) fld(i, j, k) = 0.
            end do
          end do
        end do

        ! Compute vertical sum
      case ('vertsum')
        fld(:, :, 1) = sum(fld, 3)

        ! Unit transformation: s-1 -> month-1
      case ('sec2mon')
        fld = fld * ndays(month) * 24 * 3600

      end select
      if (str1 == str2) exit
    end do

  end subroutine special_post

  ! -----------------------------------------------------------------

  subroutine read_gridinfo_ifile

    implicit none

    logical :: check
    integer :: i, j, k
    real    :: missing

    ! Open first input file
    call scan_files(reset=.true.)
    status = nf90_open(fnm, nf90_nowrite, ncid)
    call handle_ncerror(status)

    ! Read longitudes
    status = nf90_inq_dimid(ncid, 'lon', dimid)
    call handle_ncerror(status)
    status = nf90_inquire_dimension(ncid, dimid, len=idm)
    call handle_ncerror(status)
    allocate(lon(idm), lon_bnds(2, idm), stat=status)
    if (status /= 0) stop 'cannot ALLOCATE enough memory (1)'
    status = nf90_inq_varid(ncid, 'lon', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, lon)
    call handle_ncerror(status)
    lon_bnds(1, 1) = lon(1) - 0.5 * (lon(2) - lon(1))
    lon_bnds(2, 1) = lon(1) + 0.5 * (lon(2) - lon(1))
    do i = 2, idm
      lon_bnds(1, i) = lon_bnds(2, i - 1)
      lon_bnds(2, i) = lon(i) + 0.5 * (lon(2) - lon(1))
    end do

    status = nf90_inq_varid(ncid, 'lonrof', rhid)
    if (status == 0) then
      status = nf90_inq_dimid(ncid, 'lonrof', dimid)
      call handle_ncerror(status)
      status = nf90_inquire_dimension(ncid, dimid, len=idmrof)
      call handle_ncerror(status)
      allocate(lonrof(idmrof), lonrof_bnds(2, idmrof), stat=status)
      if (status /= 0) stop 'cannot ALLOCATE enough memory (1b)'
      status = nf90_inq_varid(ncid, 'lonrof', rhid)
      call handle_ncerror(status)
      status = nf90_get_var(ncid, rhid, lonrof)
      call handle_ncerror(status)
      lonrof_bnds(1, 1) = lonrof(1) - 0.5 * (lonrof(2) - lonrof(1))
      lonrof_bnds(2, 1) = lonrof(1) + 0.5 * (lonrof(2) - lonrof(1))
      do i = 2, idmrof
        lonrof_bnds(1, i) = lonrof_bnds(2, i - 1)
        lonrof_bnds(2, i) = lonrof(i) + 0.5 * (lonrof(2) - lonrof(1))
      end do
    else
      idmrof = idm
    end if

    ! Read latitudes
    status = nf90_inq_dimid(ncid, 'lat', dimid)
    call handle_ncerror(status)
    status = nf90_inquire_dimension(ncid, dimid, len=jdm)
    call handle_ncerror(status)
    allocate(lat(jdm), lat_bnds(2, jdm), stat=status)
    if (status /= 0) stop 'cannot ALLOCATE enough memory (2)'
    status = nf90_inq_varid(ncid, 'lat', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, lat)
    call handle_ncerror(status)
    lat_bnds(2, 1) = lat(1) + 90. / real(jdm - 1)
    lat_bnds(1, 1) = max(-90., lat(1) - 90. / real(jdm - 1))
    do j = 2, jdm
      lat_bnds(2, j) = min(90., lat(j) + 90. / real(jdm - 1))
      lat_bnds(1, j) = lat_bnds(2, j - 1)
    end do

    status = nf90_inq_varid(ncid, 'latrof', rhid)
    if (status == 0) then
      status = nf90_inq_dimid(ncid, 'latrof', dimid)
      call handle_ncerror(status)
      status = nf90_inquire_dimension(ncid, dimid, len=jdmrof)
      call handle_ncerror(status)
      allocate(latrof(jdmrof), latrof_bnds(2, jdmrof), stat=status)
      if (status /= 0) stop 'cannot ALLOCATE enough memory (2b)'
      status = nf90_inq_varid(ncid, 'latrof', rhid)
      call handle_ncerror(status)
      status = nf90_get_var(ncid, rhid, latrof)
      call handle_ncerror(status)
      latrof_bnds(2, 1) = latrof(1) + 90. / real(jdmrof - 1)
      latrof_bnds(1, 1) = max(-90., latrof(1) - 90. / real(jdmrof - 1))
      do j = 2, jdmrof
        latrof_bnds(2, j) = min(90., latrof(j) + 90. / real(jdmrof - 1))
        latrof_bnds(1, j) = latrof_bnds(2, j - 1)
      end do
    else
      jdmrof = jdm
    end if

    ! Read soil depths
    status = nf90_inq_dimid(ncid, 'levgrnd', dimid)
    if (status == 0) then
      status = nf90_inquire_dimension(ncid, dimid, len=kdm)
      call handle_ncerror(status)
      allocate(lev(kdm), lev_bnds(2, kdm), stat=status)
      if (status /= 0) stop 'cannot allocate enough memory'
      status = nf90_inq_varid(ncid, 'levgrnd', rhid)
      call handle_ncerror(status)
      status = nf90_get_var(ncid, rhid, lev)
      call handle_ncerror(status)
    else
      write(*, *) "WARNING: dimension 'levgrnd' is not defined"
    end if

    ! Read landuse types
    status = nf90_inq_dimid(ncid, 'landUse', dimid)
    if (status == 0) then
      status = nf90_inquire_dimension(ncid, dimid, len=ldm)
      call handle_ncerror(status)
      allocate(ltype(ldm), stat=status)
      if (status /= 0) stop 'cannot allocate enough memory'
      if (ldm == 4) then
        ltype(1) = "primary_and_secondary_land"
        ltype(2) = "pastures"
        ltype(3) = "crops"
        ltype(4) = "urban"
      else
        write(*, *) "The landUse type is hard-coded with four &
          &elements. Need update the defined variable ltype"
        stop
      end if
    end if

    ! Compute level bounds from soil thickness
    if (allocated(fld)) deallocate(fld)
    allocate(fld(idm, jdm, kdm), stat=status)
    if (status /= 0) stop 'cannot allocate enough memory'
    status = nf90_inq_varid(ncid, 'DZSOI', rhid)
    if (status == nf90_noerr) then
      call handle_ncerror(status)
      status = nf90_get_att(ncid, rhid, '_FillValue', missing)
      call handle_ncerror(status)
      status = nf90_get_var(ncid, rhid, fld)
      call handle_ncerror(status)
      check = .false.
      do i = 1, idm
        do j = 1, jdm
          if (fld(i, j, 1) /= missing) then
            lev_bnds(1, 1) = 0.
            lev_bnds(2, 1) = fld(i, j, 1)
            do k = 2, kdm
              lev_bnds(1, k) = lev_bnds(1, k - 1) + fld(i, j, k - 1)
              lev_bnds(2, k) = lev_bnds(1, k) + fld(i, j, k)
            end do
            check = .true.
            exit
          end if
          if (check) exit
        end do
      end do
    else
      if (allocated(lev)) then
        write(*, *) 'WARNING: cannot find varible DZSOI. will try to &
          &guess soil depth bounds.'
        lev_bnds(1, 1) = 0.
        lev_bnds(2, 1) = 0.5 * (lev(1) + lev(2))
        do k = 2, kdm - 1
          lev_bnds(1, k) = lev_bnds(2, k - 1)
          lev_bnds(2, k) = 0.5 * (lev(k) + lev(k + 1))
        end do
        lev_bnds(1, kdm) = lev_bnds(2, kdm - 1)
        lev_bnds(2, kdm) = lev_bnds(1, kdm) + lev(kdm) - lev(kdm - 1)
      end if
    end if
    deallocate(fld)

    ! Read calendar info
    status = nf90_inq_varid(ncid, 'time', rhid)
    call handle_ncerror(status)
    status = nf90_get_att(ncid, rhid, 'calendar', calendar)
    call handle_ncerror(status)
    write(calunits(12:15), '(i4.4)') exprefyear

    ! Close file
    status = nf90_close(ncid)
    call handle_ncerror(status)

  end subroutine read_gridinfo_ifile

  ! -----------------------------------------------------------------

  subroutine open_ofile(fx)

    implicit none

    logical, optional, intent(in)   :: fx
    logical                         :: fxflag

    real                            :: fac1, fac2, fac3, fac4, fac5, fac6
    integer, parameter              :: ndimmax = 10
    integer                         :: n, ndims, dimids(ndimmax), dimlens(ndimmax)
    integer                         :: physics_version = 1, initialization_method = 1
    character(len=slenmax)          :: ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6

    ! Check if output variable should have time coordinate
    fxflag = .false.
    if (present(fx)) then
      if (fx) fxflag = .true.
    end if

    ! Inquire variable units and dimensions in input file
    status = nf90_open(fnm, nf90_nowrite, ncid)
    call handle_ncerror(status)

    call resolve_vnm(slenmax, ivnm, ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6, &
      fac1, fac2, fac3, fac4, fac5, fac6)
    status = nf90_inq_varid(ncid, trim(ivnm1), rhid)
    if (status /= nf90_noerr) then
      write(*, *) 'cannot find input variable ', trim(ivnm1)
      stop
    end if
    status = nf90_inquire_variable(ncid, rhid, ndims=ndims)
    call handle_ncerror(status)
    if (.not. fxflag .and. ndims < 3) then
      write(*, *) 'Variable ', trim(ivnm1), ' has too few dimensions', ndims
    end if
    status = nf90_inquire_variable(ncid, rhid, dimids=dimids(1:ndims))
    call handle_ncerror(status)
    dimlens = 1
    do n = 1, ndims
      status = nf90_inquire_dimension(ncid, dimids(n), len=dimlens(n))
      call handle_ncerror(status)
    end do
    if (dimlens(1) /= idm .and. dimlens(1) /= idmrof) then
      write(*, *) 'unexpected first dimension of variable ', &
        trim(ivnm1), ': ', dimlens(1), ' versus idm=', idm
      stop
    end if
    if (dimlens(2) /= jdm .and. dimlens(2) /= jdmrof) then
      write(*, *) 'unexpected second dimension of variable ', &
        trim(ivnm1), ': ', dimlens(2), ' versus jdm=', idm
      stop
    end if
    if (fxflag .and. ndims > 2 .or. ndims > 3) then
      kdm = dimlens(3)
    else
      kdm = 1
    end if
    if (allocated(fld)) deallocate(fld, fld2, fldacc)
    if (index(ivnm, 'QCHOCNR') > 0) then
      ii = idmrof
      jj = jdmrof
    else
      ii = idm
      jj = jdm
    end if
    allocate(fld(ii, jj, kdm), fld2(ii, jj, kdm), fldacc(ii, jj, kdm), stat=status)
    if (status /= 0) stop 'cannot ALLOCATE enough memory (4)'

    if (len_trim(vunits) == 0) then
      status = nf90_get_att(ncid, rhid, 'units', vunits)
      call handle_ncerror(status)
      if (trim(vunits) == 'mm/s') vunits = 'kg m-2 s-1'
    end if

    status = nf90_close(ncid)
    call handle_ncerror(status)

    ! Inquire time and vertical dimension of output variable
    if (.not. fxflag) then
      call get_timecoord(trim(tabledir)//trim(table), ovnm, tcoord)
    end if
    call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)

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
    call write_namelist_json(lndgrid, lndgrid_label, lndgrid_resolution, ovnm)
    error_flag = cmor_dataset_json(namelist_file_json)
    call system('rm '//trim(namelist_file_json))

    ! Define horizontal axes
    if (index(ivnm, 'QCHOCNR') > 0) then
      iaxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry='longitude', &
        units='degrees_east', &
        length=idmrof, &
        coord_vals=lonrof, &
        cell_bounds=lonrof_bnds)
      jaxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry='latitude', &
        units='degrees_north', &
        length=jdmrof, &
        coord_vals=latrof, &
        cell_bounds=latrof_bnds)
    else
      iaxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry='longitude', &
        units='degrees_east', &
        length=idm, &
        coord_vals=lon, &
        cell_bounds=lon_bnds)
      jaxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry='latitude', &
        units='degrees_north', &
        length=jdm, &
        coord_vals=lat, &
        cell_bounds=lat_bnds)
    end if

    ! Define time axis
    if (.not. fxflag) then
      taxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry=trim(tcoord), &
        units=trim(calunits), &
        length=1)
    end if

    ! Define vertical axis
    if (trim(zcoord) == 'sdepth') then
      kaxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry=trim(zcoord), &
        units='m', &
        length=kdm, &
        coord_vals=lev, &
        cell_bounds=lev_bnds)
    end if
    if (trim(zcoord) == 'landUse') then
      kaxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry=trim(zcoord), &
        units='non', &
        length=ldm, &
        coord_vals=ltype)
    end if

    ! Define output variable
    if (fxflag) then
      varid = cmor_variable( &
        table_entry=trim(ovnm), &
        units=trim(vunits), &
        axis_ids=(/iaxid, jaxid/), &
        missing_value=1e20, &
        original_name=trim(ivnm))
    else
      if (trim(zcoord) == 'sdepth' .or. trim(zcoord) == 'vegtype' &
        .or. trim(zcoord) == 'landUse') then
        varid = cmor_variable( &
          table=trim(tabledir)//trim(table), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/iaxid, jaxid, kaxid, taxid/), &
          original_name=trim(ivnm), &
          missing_value=1e20, &
          positive=trim(vpositive))
      else
        varid = cmor_variable( &
          table=trim(tabledir)//trim(table), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/iaxid, jaxid, taxid/), &
          original_name=trim(ivnm), &
          missing_value=1e20, &
          positive=trim(vpositive))
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

    real                    :: fac1, fac2, fac3, fac4, fac5, fac6
    integer                 :: ind
    character(len=slenmax)  :: ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6

    ! Open input file
    status = nf90_open(fnm, nf90_nowrite, ncid)
    call handle_ncerror(status)

    ! Read data
    call resolve_vnm(slenmax, ivnm, ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6, &
      fac1, fac2, fac3, fac4, fac5, fac6)
    if (verbose .and. len_trim(ivnm2) /= 0) &
      write(*, *) 'Compound variable: ', trim(ivnm1), '*', fac1, '+', &
      trim(ivnm2), '*', fac2, '+', trim(ivnm3), '*', fac3, &
      trim(ivnm4), '*', fac4, '+', trim(ivnm5), '*', fac5, &
      trim(ivnm6), '*', fac6
    status = nf90_inq_varid(ncid, trim(ivnm1), rhid)
    if (status /= nf90_noerr) then
      write(*, *) 'cannot find input variable ', trim(ivnm1)
      stop
    end if
    status = nf90_get_var(ncid, rhid, fld)
    call handle_ncerror(status)
    if (fac1 /= 1) then
      fld = fld * fac1
    end if

    if (len_trim(ivnm2) > 0) then
      status = nf90_inq_varid(ncid, trim(ivnm2), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm2)
        stop
      end if
      status = nf90_get_var(ncid, rhid, fld2)
      call handle_ncerror(status)
      fld = fld + fld2 * fac2
    end if

    if (len_trim(ivnm3) > 0) then
      status = nf90_inq_varid(ncid, trim(ivnm3), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm3)
        stop
      end if
      status = nf90_get_var(ncid, rhid, fld2)
      call handle_ncerror(status)
      fld = fld + fld2 * fac3
    end if

    if (len_trim(ivnm4) > 0) then
      status = nf90_inq_varid(ncid, trim(ivnm4), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm4)
        stop
      end if
      status = nf90_get_var(ncid, rhid, fld2)
      call handle_ncerror(status)
      fld = fld + fld2 * fac4
    end if

    if (len_trim(ivnm5) > 0) then
      status = nf90_inq_varid(ncid, trim(ivnm5), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm5)
        stop
      end if
      status = nf90_get_var(ncid, rhid, fld2)
      call handle_ncerror(status)
      fld = fld + fld2 * fac5
    end if

    if (len_trim(ivnm6) > 0) then
      status = nf90_inq_varid(ncid, trim(ivnm6), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm6)
        stop
      end if
      status = nf90_get_var(ncid, rhid, fld2)
      call handle_ncerror(status)
      fld = fld + fld2 * fac6
    end if

    status = nf90_close(ncid)
    call handle_ncerror(status)

  end subroutine read_field

  ! -----------------------------------------------------------------

  subroutine read_tslice(rec, badrec, fname)

    implicit none

    real                    :: fac1, fac2, fac3, fac4, fac5, fac6
    integer, intent(in)     :: rec
    logical, intent(out)    :: badrec
    character(len=*), intent(in), optional :: fname
    integer, save           :: fid
    character(len=slenmax)  :: ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6

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
      ! correct erroneous intial time bound
      tbnds(1, 1) = max(0., tbnds(1, 1))
      tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))
    end if

    ! Read data
    call resolve_vnm(slenmax, ivnm, ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6, &
      fac1, fac2, fac3, fac4, fac5, fac6)
    if (verbose .and. rec == 1 .and. len_trim(ivnm2) /= 0) &
      write(*, *) 'Compound variable: ', trim(ivnm1), '*', fac1, '+', &
      trim(ivnm2), '*', fac2, '+', trim(ivnm3), '*', fac3, &
      trim(ivnm4), '*', fac4, '+', trim(ivnm5), '*', fac5, &
      trim(ivnm6), '*', fac6
    status = nf90_inq_varid(fid, trim(ivnm1), rhid)
    if (status /= nf90_noerr) then
      write(*, *) 'cannot find input variable ', trim(ivnm1)
      stop
    end if
    if (kdm == 1) then
      if (index(ivnm, 'QCHOCNR') > 0) then
        status = nf90_get_var(fid, rhid, fld, (/1, 1, rec/), (/idmrof, jdmrof, 1/))
      else
        status = nf90_get_var(fid, rhid, fld, (/1, 1, rec/), (/idm, jdm, 1/))
      end if
    else
      status = nf90_get_var(fid, rhid, fld, (/1, 1, 1, rec/), (/idm, jdm, kdm, 1/))
    end if
    call handle_ncerror(status)
    if (fac1 /= 1) then
      fld = fld * fac1
    end if

    if (len_trim(ivnm2) > 0) then
      status = nf90_inq_varid(fid, trim(ivnm2), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm2)
        stop
      end if
      if (kdm == 1) then
        if (index(ivnm, 'QCHOCNR') > 0) then
          status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idmrof, jdmrof, 1/))
        else
          status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idm, jdm, 1/))
        end if
      else
        status = nf90_get_var(fid, rhid, fld2, (/1, 1, 1, rec/), (/idm, jdm, kdm, 1/))
      end if
      call handle_ncerror(status)
      fld = fld + fld2 * fac2
    end if

    if (len_trim(ivnm3) > 0) then
      status = nf90_inq_varid(fid, trim(ivnm3), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm3)
        stop
      end if
      if (kdm == 1) then
        if (index(ivnm, 'QCHOCNR') > 0) then
          status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idmrof, jdmrof, 1/))
        else
          status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idm, jdm, 1/))
        end if
      else
        status = nf90_get_var(fid, rhid, fld2, (/1, 1, 1, rec/), (/idm, jdm, kdm, 1/))
      end if
      call handle_ncerror(status)
      fld = fld + fld2 * fac3
    end if

    if (len_trim(ivnm4) > 0) then
      status = nf90_inq_varid(fid, trim(ivnm4), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm4)
        stop
      end if
      if (kdm == 1) then
        if (index(ivnm, 'QCHOCNR') > 0) then
          status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idmrof, jdmrof, 1/))
        else
          status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idm, jdm, 1/))
        end if
      else
        status = nf90_get_var(fid, rhid, fld2, (/1, 1, 1, rec/), (/idm, jdm, kdm, 1/))
      end if
      call handle_ncerror(status)
      fld = fld + fld2 * fac4
    end if

    if (len_trim(ivnm5) > 0) then
      status = nf90_inq_varid(fid, trim(ivnm5), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm5)
        stop
      end if
      if (kdm == 1) then
        if (index(ivnm, 'QCHOCNR') > 0) then
          status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idmrof, jdmrof, 1/))
        else
          status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idm, jdm, 1/))
        end if
      else
        status = nf90_get_var(fid, rhid, fld2, (/1, 1, 1, rec/), (/idm, jdm, kdm, 1/))
      end if
      call handle_ncerror(status)
      fld = fld + fld2 * fac5
    end if

    if (len_trim(ivnm6) > 0) then
      status = nf90_inq_varid(fid, trim(ivnm6), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm6)
        stop
      end if
      if (kdm == 1) then
        if (index(ivnm, 'QCHOCNR') > 0) then
          status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idmrof, jdmrof, 1/))
        else
          status = nf90_get_var(fid, rhid, fld2, (/1, 1, rec/), (/idm, jdm, 1/))
        end if
      else
        status = nf90_get_var(fid, rhid, fld2, (/1, 1, 1, rec/), (/idm, jdm, kdm, 1/))
      end if
      call handle_ncerror(status)
      fld = fld + fld2 * fac6
    end if

    status = nf90_close(fid)
    call handle_ncerror(status)

  end subroutine read_tslice

  ! -----------------------------------------------------------------

  subroutine write_field

    implicit none

    integer :: i, j, k

    ! Set zero on ocean grid cells
    do k = 1, kdm
      do j = 1, jj
        do i = 1, ii
          if (abs(fld(i, j, k)) > 2e20) fld(i, j, k) = 0.
        end do
      end do
    end do

    ! Store variable
    error_flag = cmor_write( &
      var_id=varid, &
      data=reshape(fld, (/idm, jdm/)))

  end subroutine write_field

  ! -----------------------------------------------------------------

  subroutine write_tslice

    implicit none

    integer :: i, j, k

    ! Set zero on ocean grid cells
    do k = 1, kdm
      do j = 1, jj
        do i = 1, ii
          if (abs(fld(i, j, k)) > 2e20) fld(i, j, k) = 0.
        end do
      end do
    end do

    ! Store variable
    if (len_trim(zcoord) > 0 .or. kdm == 1) then
      if (trim(tcoord) /= 'time1') then
        if (trim(zcoord) == 'typecrop') then
          error_flag = cmor_write( &
            var_id=varid, &
            data=fld(1:idm, 1:jdm, 2), &
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
      else
        error_flag = cmor_write( &
          var_id=varid, &
          data=fld, &
          ntimes_passed=1, &
          time_vals=tval)
      end if
    else
      if (trim(tcoord) /= 'time1') then
        error_flag = cmor_write( &
          var_id=varid, &
          data=fld(1:idm, 1:jdm, 1), &
          ntimes_passed=1, &
          time_vals=tval, &
          time_bnds=tbnds)
      else
        error_flag = cmor_write( &
          var_id=varid, &
          data=fld(1:idm, 1:jdm, 1), &
          ntimes_passed=1, &
          time_vals=tval)
      end if
    end if

  end subroutine write_tslice

end module m_modelslnd
