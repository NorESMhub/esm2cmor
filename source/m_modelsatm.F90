module m_modelsatm

  use netcdf
  use cmor_users_functions
  use m_utilities
  use m_namelists

  implicit none

  ! Netcdf variables
  integer :: ncid, rhid, dimid, status

  ! Grid dimensions and variables
  integer, save :: idm, jdm, kdm, ldm, pdm
  real(kind=8), save :: p0 = 1000.d0
  real(kind=8), allocatable, save, dimension(:) :: lon, slon, lat, slat, lev, &
    ilev, hyam, hybm, hyai, hybi, plevs, iilev, hyaii, hybii
  real(kind=8), allocatable, save, dimension(:,:) :: lon_bnds, lat_bnds
  character(len=slenmax), save :: zcoord, tcoord

  ! Dataset related variables
  logical, save :: lreqphis, lreqsst, lreadplev
  character(len=slenmax), save :: ivnm, ovnm, vunits
  character(len=slenmax), save :: vcomment
  character(len=slenmax), save :: vpositive

  ! String for module special
  character(len=slenmax), save :: special

  ! Table related variables
  character(len=lenmax) :: table

  ! Cmor parameters
  integer, save :: iaxid, jaxid, kaxid, taxid, varid, zfacid, table_id, error_flag

  ! Data fields
  real(kind=8), allocatable, save, dimension(:,:)   :: ps, phis, sst, tbot
  real(kind=8), allocatable, save, dimension(:,:,:) :: ifld, ifld2, ofld, ifldacc

  ! Gravity (same as in nco_lifetimes_A2.csh by Alf Kirkevaag)
  real(kind=8), parameter :: g = 9.80665, ginv = 1.d0 / g

  ! Auxillary variables for special operations
  character(len=slenmax), save :: str1, str2

contains

  ! -----------------------------------------------------------------

  subroutine atm2cmor

    implicit none

    logical :: badrec, last, skipini
    integer :: m, n, nrec

    badrec = .false.

    ! Print start information
    if (verbose) then
      write(*, *)
      write(*, *) '----------------------------------'
      write(*, *) '--- Process atmospheric output ---'
      write(*, *) '----------------------------------'
      write(*, *)
    end if

    ! Read grid information from input files
    itag = tagamon
    call scan_files(reset=.true.)
    if (len_trim(fnm) == 0) then
      write(*, *) 'ERROR: fnm is empty'
      return
    end if
    call read_gridinfo_ifile

    ! Initialise pressure level dimension
    pdm = 1
    allocate(plevs(pdm))

    ! Process table fx
    write(*, *) 'Process table fx'
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
      vcomment = ' '

      write(*, *) 'tabledir:', tabledir

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      write(*, *) 'zcoord:', trim(zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Check if input variable is present
      if (len_trim(pfx) == 0) then
        if (index(trim(ovnm), 'sftlf') > 0) then
          call scan_files(reset=.true.)
        else
          fnm = trim(griddata)//trim(atmgridfile)
        end if
        if (index(trim(ovnm), 'sftgif') > 0) then
          call scan_files(reset=.true.)
          fnm2 = fnm
          fnm = trim(griddata)//trim(atmgridfile)
        end if
      end if
      write(*, *) 'ivnm, ovnm:'//trim(ivnm)//','//trim(ovnm)
      write(*, *) 'fnm:'//trim(fnm)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre
      call open_ofile(fx=.true.)

      ! Read field
      call read_field

      ! Post Processing
      ofld = ifld
      call special_post

      ! Write field
      call write_field

      ! Close output file
      call close_ofile

    end do

    ! Process table amon
    write(*, *) 'Process table amon'
    fnm = pamon
    table = tamon
    linstant = .false.
    lreadplev = .true.
    do n = 1, namon
      if (skip_variable(n, namon, damon)) cycle

      ! Map namelist variables
      ovnm = vamon(ovnmpos, n)
      ivnm = vamon(ivnmpos, n)
      special = vamon(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      write(*, *) 'zcoord(1:4):', zcoord(1:4)
      !write(*, *) trim(ovnm), trim(zcoord), do_3d
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
        ! Read pressure levels from table
        if (lreadplev .and. zcoord(1:4) == 'plev') then
          call read_gridinfo_plev(trim(zcoord))
        end if
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      select case (trim(special))
      case ('day2mon')
        itag = tagaday
      case ('6hr2mon')
        itag = taga6hr
      case ('3hr2mon')
        itag = taga3hr
      case default
        itag = tagamon
      end select

      ! Check if input variable is present
      if (len_trim(pamon) == 0) call scan_files(reset=.true.)
      if (index(special, 'catplev') < 1) then
        if (.not. var_in_file(fnm, ivnm)) cycle
      end if

      ! Prepare special operations
      call special_pre

      ! Loop over input files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, ramon) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        nrec = 0
        ifldacc = 0.
        last = .false.
        do
          if (len_trim(pamon) == 0) call scan_files
          if (rec == 0) then
            last = .true.
            exit
          end if
          nrec = nrec + 1
          call read_tslice(rec, badrec, fnm)
          ifldacc = ifldacc + ifld
          if (tbnd(2) == mbnd(2)) exit
        end do
        if (last) exit
        ifld = ifldacc / real(nrec)
        tbnds(:, 1) = mbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Post processing
        if (zcoord(1:4) == 'plev' .and. index(special, 'catplev') < 1) then
          call interp_z
        else
          if (.not. (index(special, 'ps2pfull') > 0 .or. &
            index(special, 'ps2phalf') > 0)) then
            ofld = ifld
          end if
        end if
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, ramon) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, ramon) > 0) call close_ofile

    end do

    ! Process table Emon
    write(*, *) 'Process table Emon'
    fnm = pEmon
    table = tEmon
    linstant = .false.
    lreadplev = .true.
    do n = 1, nEmon
      if (skip_variable(n, nEmon, dEmon)) cycle

      ! Map namelist variables
      ovnm = vEmon(ovnmpos, n)
      ivnm = vEmon(ivnmpos, n)
      special = vEmon(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        ! Read pressure levels from table
        if (lreadplev .and. zcoord(1:4) == 'plev') then
          call read_gridinfo_plev(trim(zcoord))
        end if
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      select case (trim(special))
      case ('day2mon')
        itag = tagaday
      case ('6hr2mon')
        itag = taga6hr
      case ('3hr2mon')
        itag = taga3hr
      case default
        itag = tagamon
      end select

      ! Check if input variable is present
      if (len_trim(pEmon) == 0) call scan_files(reset=.true.)
      if (index(special, 'catplev') < 1) then
        if (.not. var_in_file(fnm, ivnm)) cycle
      end if

      ! Prepare special operations
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
        ifldacc = 0.
        last = .false.
        do
          if (len_trim(pEmon) == 0) call scan_files
          if (rec == 0) then
            last = .true.
            exit
          end if
          nrec = nrec + 1
          call read_tslice(rec, badrec, fnm)
          ifldacc = ifldacc + ifld
          if (tbnd(2) == mbnd(2)) exit
        end do
        if (last) exit
        ifld = ifldacc / real(nrec)
        tbnds(:, 1) = mbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Post processing
        if (zcoord(1:4) == 'plev' .and. index(special, 'catplev') < 1) then
          call interp_z
        else
          ofld = ifld
        end if
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, rEmon) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rEmon) > 0) call close_ofile

    end do

    ! Process table EmonZ
    write(*, *) 'Process table EmonZ'
    fnm = pEmonZ
    table = tEmonZ
    linstant = .false.
    lreadplev = .true.
    do n = 1, nEmonZ
      if (skip_variable(n, nEmonZ, dEmonZ)) cycle

      ! Map namelist variables
      ovnm = vEmonZ(ovnmpos, n)
      ivnm = vEmonZ(ivnmpos, n)
      special = vEmonZ(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
        ! Read pressure levels from table
        if (lreadplev .and. zcoord(1:4) == 'plev') then
          call read_gridinfo_plev(trim(zcoord))
        end if
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      select case (trim(special))
      case ('day2mon')
        itag = tagaday
      case ('6hr2mon')
        itag = taga6hr
      case ('3hr2mon')
        itag = taga3hr
      case default
        itag = tagamon
      end select

      ! Check if input variable is present
      if (len_trim(pEmonZ) == 0) call scan_files(reset=.true.)
      if (index(special, 'catplev') < 1) then
        if (.not. var_in_file(fnm, ivnm)) cycle
      end if

      ! Prepare special operations
      call special_pre

      ! Loop over input files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rEmonZ) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        nrec = 0
        ifldacc = 0.
        last = .false.
        do
          if (len_trim(pEmonZ) == 0) call scan_files
          if (rec == 0) then
            last = .true.
            exit
          end if
          nrec = nrec + 1
          call read_tslice(rec, badrec, fnm)
          ifldacc = ifldacc + ifld
          if (tbnd(2) == mbnd(2)) exit
        end do
        if (last) exit
        ifld = ifldacc / real(nrec)
        tbnds(:, 1) = mbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Post processing
        if (zcoord(1:4) == 'plev' .and. index(special, 'catplev') < 1) then
          call interp_z
        else
          ofld = ifld
        end if
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, rEmonZ) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rEmonZ) > 0) call close_ofile

    end do

    ! Process table aero
    write(*, *) 'Process table aero'
    fnm = paero
    table = taero
    linstant = .false.
    lreadplev = .true.
    do n = 1, naero
      if (skip_variable(n, naero, daero)) cycle

      ! Map namelist variables
      ovnm = vaero(ovnmpos, n)
      ivnm = vaero(ivnmpos, n)
      special = vaero(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
        ! Read pressure levels from table
        if (lreadplev .and. zcoord(1:4) == 'plev') then
          call read_gridinfo_plev(trim(zcoord))
        end if
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      select case (trim(special))
      case ('day2mon')
        itag = tagaday
      case ('6hr2mon')
        itag = taga6hr
      case ('3hr2mon')
        itag = taga3hr
      case default
        itag = tagamon
      end select

      ! Check if input variable is present
      if (len_trim(paero) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over files and records
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, raero) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        nrec = 0
        ifldacc = 0.
        last = .false.
        do
          if (len_trim(paero) == 0) call scan_files
          if (rec == 0) then
            last = .true.
            exit
          end if
          nrec = nrec + 1
          call read_tslice(rec, badrec, fnm)
          ifldacc = ifldacc + ifld
          if (tbnd(2) == mbnd(2)) exit
        end do
        if (last) exit

        ifld = ifldacc / real(nrec)
        tbnds(:, 1) = mbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Post processing
        if (.not. (index(special, 'ps2pfull') > 0 .or. &
          index(special, 'ps2phalf') > 0)) then
          ofld = ifld
        end if
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, raero) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, raero) > 0) call close_ofile

    end do

    ! Process table CFmon
    write(*, *) 'Process table CFmon'
    fnm = pcfmon
    table = tcfmon
    linstant = .false.
    lreadplev = .true.
    do n = 1, ncfmon
      if (skip_variable(n, ncfmon, dcfmon)) cycle

      ! Map namelist variables
      ovnm = vcfmon(ovnmpos, n)
      ivnm = vcfmon(ivnmpos, n)
      special = vcfmon(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
        ! Read pressure levels from table
        if (lreadplev .and. zcoord(1:4) == 'plev') then
          call read_gridinfo_plev(trim(zcoord))
        end if
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      select case (trim(special))
      case ('day2mon')
        itag = tagaday
      case ('6hr2mon')
        itag = taga6hr
      case ('3hr2mon')
        itag = taga3hr
      case default
        itag = tagamon
      end select

      ! Check if input variable is present
      if (len_trim(pcfmon) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over files and records
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rcfmon) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        nrec = 0
        ifldacc = 0.
        last = .false.
        do
          if (len_trim(pcfmon) == 0) call scan_files
          if (rec == 0) then
            last = .true.
            exit
          end if
          nrec = nrec + 1
          call read_tslice(rec, badrec, fnm)
          ifldacc = ifldacc + ifld
          if (tbnd(2) == mbnd(2)) exit
        end do
        if (last) exit

        ifld = ifldacc / real(nrec)
        tbnds(:, 1) = mbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Post processing
        ofld = ifld
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, rcfmon) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rcfmon) > 0) call close_ofile

    end do

    ! Process table day
    write(*, *) 'Process table day'
    fnm = pday
    table = tday
    linstant = .false.
    lreadplev = .true.
    do n = 1, nday
      if (skip_variable(n, nday, dday)) cycle

      ! Map namelist variables
      ovnm = vday(ovnmpos, n)
      ivnm = vday(ivnmpos, n)
      special = vday(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
        ! Read pressure levels from table
        if (lreadplev .and. zcoord(1:4) == 'plev') then
          call read_gridinfo_plev(trim(zcoord))
        end if
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = tagaday

      ! Check if input variable is present
      if (len_trim(pday) == 0) call scan_files(reset=.true.)
      if (index(special, 'catplev') < 1) then
        if (.not. var_in_file(fnm, ivnm)) cycle
      end if

      ! Prepare output file
      call special_pre

      ! Loop over files
      skipini = .false.
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rday) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        if (len_trim(pday) == 0) call scan_files
        if (rec == 0) exit
        if (tbnd(1) == tbnd(2)) then
          if (verbose) write(*, *) 'Skip initialisation record'
          skipini = .true.
          cycle
        end if
        tbnds(:, 1) = tbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        if (zcoord(1:4) == 'plev' .and. index(special, 'catplev') < 1) then
          call interp_z
        else
          ofld = ifld
        end if
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (skipini .and. mod(m - 1, rday) == 0 .or. &
          .not. skipini .and. mod(m, rday) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rday) > 0) call close_ofile

    end do

    ! Process table AERday
    write(*, *) 'Process table AERday'
    fnm = pAERday
    table = tAERday
    linstant = .false.
    lreadplev = .true.
    do n = 1, nAERday
      if (skip_variable(n, nAERday, dAERday)) cycle

      ! Map namelist variables
      ovnm = vAERday(ovnmpos, n)
      ivnm = vAERday(ivnmpos, n)
      special = vAERday(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
        ! Read pressure levels from table
        if (lreadplev .and. zcoord(1:4) == 'plev') then
          call read_gridinfo_plev(trim(zcoord))
        end if
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = tagaday

      ! Check if input variable is present
      if (len_trim(pAERday) == 0) call scan_files(reset=.true.)
      if (index(special, 'catplev') < 1) then
        if (.not. var_in_file(fnm, ivnm)) cycle
      end if

      ! Prepare output file
      call special_pre

      ! Loop over files
      skipini = .false.
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rAERday) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        if (len_trim(pAERday) == 0) call scan_files
        if (rec == 0) exit
        if (tbnd(1) == tbnd(2)) then
          if (verbose) write(*, *) 'Skip initialisation record'
          skipini = .true.
          cycle
        end if
        tbnds(:, 1) = tbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        if (zcoord(1:4) == 'plev' .and. index(special, 'catplev') < 1) then
          call interp_z
        else
          ofld = ifld
        end if
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (skipini .and. mod(m - 1, rAERday) == 0 .or. &
          .not. skipini .and. mod(m, rAERday) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rAERday) > 0) call close_ofile

    end do

    ! Process table Eday
    write(*, *) 'Process table Eday'
    fnm = pEday
    table = tEday
    linstant = .false.
    lreadplev = .true.
    do n = 1, nEday
      if (skip_variable(n, nEday, dEday)) cycle

      ! Map namelist variables
      ovnm = vEday(ovnmpos, n)
      ivnm = vEday(ivnmpos, n)
      special = vEday(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
        ! Read pressure levels from table
        if (lreadplev .and. zcoord(1:4) == 'plev') then
          call read_gridinfo_plev(trim(zcoord))
        end if
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = tagaday

      ! Check if input variable is present
      if (len_trim(pEday) == 0) call scan_files(reset=.true.)
      if (index(special, 'catplev') < 1) then
        if (.not. var_in_file(fnm, ivnm)) cycle
      end if

      ! Prepare output file
      call special_pre

      ! Loop over files
      skipini = .false.
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rEday) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        if (len_trim(pEday) == 0) call scan_files
        if (rec == 0) exit
        if (tbnd(1) == tbnd(2)) then
          if (verbose) write(*, *) 'Skip initialisation record'
          skipini = .true.
          cycle
        end if
        tbnds(:, 1) = tbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        if (zcoord(1:4) == 'plev' .and. index(special, 'catplev') < 1) then
          call interp_z
        else
          ofld = ifld
        end if
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (skipini .and. mod(m - 1, rEday) == 0 .or. &
          .not. skipini .and. mod(m, rEday) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rEday) > 0) call close_ofile

    end do

    ! Process table CFday
    write(*, *) 'Process table CFday'
    fnm = pCFday
    table = tCFday
    linstant = .false.
    lreadplev = .true.
    do n = 1, nCFday
      if (skip_variable(n, nCFday, dCFday)) cycle

      ! Map namelist variables
      ovnm = vCFday(ovnmpos, n)
      ivnm = vCFday(ivnmpos, n)
      special = vCFday(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
        ! Read pressure levels from table
        if (lreadplev .and. zcoord(1:4) == 'plev') then
          call read_gridinfo_plev(trim(zcoord))
        end if
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = tagaday

      ! Check if input variable is present
      if (len_trim(pCFday) == 0) call scan_files(reset=.true.)
      if (index(special, 'catplev') < 1) then
        if (.not. var_in_file(fnm, ivnm)) cycle
      end if

      ! Prepare output file
      call special_pre

      ! Loop over files
      skipini = .false.
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rCFday) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        if (len_trim(pCFday) == 0) call scan_files
        if (rec == 0) exit
        if (tbnd(1) == tbnd(2)) then
          if (verbose) write(*, *) 'Skip initialisation record'
          skipini = .true.
          cycle
        end if
        tbnds(:, 1) = tbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        if (zcoord(1:4) == 'plev' .and. index(special, 'catplev') < 1) then
          call interp_z
        else
          if (.not. (index(special, 'ps2pfull') > 0 .or. &
            index(special, 'ps2phalf') > 0)) then
            ofld = ifld
          end if
        end if
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (skipini .and. mod(m - 1, rCFday) == 0 .or. &
          .not. skipini .and. mod(m, rCFday) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rCFday) > 0) call close_ofile

    end do

    ! Process table 6hrlev
    write(*, *) 'Process table 6hrlev'
    fnm = p6hrlev
    table = t6hrlev
    linstant = .false.
    do n = 1, n6hrlev
      if (skip_variable(n, n6hrlev, d6hrlev)) cycle

      ! Map namelist variables
      ovnm = v6hrlev(ovnmpos, n)
      ivnm = v6hrlev(ivnmpos, n)
      special = v6hrlev(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = taga6hr

      ! Check if input variable is present
      if (len_trim(p6hrlev) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, r6hrlev) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        if (len_trim(p6hrlev) == 0) call scan_files
        if (rec == 0) exit
        tbnds(:, 1) = tbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        if (.not. (index(special, 'ps2pfull') > 0 .or. &
          index(special, 'ps2phalf') > 0)) then
          ofld = ifld
        end if
        call special_post
        if (badrec) ofld = 1e20

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, r6hrlev) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, r6hrlev) > 0) call close_ofile

    end do

    ! Process table 6hrlevi
    write(*, *) 'Process table 6hrlevi'
    fnm = p6hrlevi
    table = t6hrlevi
    linstant = .true.
    do n = 1, n6hrlevi
      if (skip_variable(n, n6hrlevi, d6hrlevi)) cycle

      ! Map namelist variables
      ovnm = v6hrlevi(ovnmpos, n)
      ivnm = v6hrlevi(ivnmpos, n)
      special = v6hrlevi(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = taga6hri

      ! Check if input variable is present
      if (len_trim(p6hrlevi) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, r6hrlevi) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        if (len_trim(p6hrlevi) == 0) call scan_files
        if (rec == 0) exit
        tbnds(:, 1) = tbnd
        tval = tbnds(2, 1)

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        ofld = ifld
        call special_post
        if (badrec) ofld = 1e20

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, r6hrlevi) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, r6hrlevi) > 0) call close_ofile

    end do

    ! Process table 6hrplev
    write(*, *) 'Process table 6hrplev'
    fnm = p6hrplev
    table = t6hrplev
    linstant = .false.
    lreadplev = .true.
    do n = 1, n6hrplev
      if (skip_variable(n, n6hrplev, d6hrplev)) cycle

      ! Map namelist variables
      ovnm = v6hrplev(ovnmpos, n)
      ivnm = v6hrplev(ivnmpos, n)
      special = v6hrplev(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
        ! Read pressure levels from table
        if (lreadplev .and. zcoord(1:4) == 'plev') then
          call read_gridinfo_plev(trim(zcoord))
        end if
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = taga6hr

      ! Check if input variable is present
      if (len_trim(p6hrplev) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, r6hrplev) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        if (len_trim(p6hrplev) == 0) call scan_files
        if (rec == 0) exit
        tbnds(:, 1) = tbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        if (zcoord(1:4) == 'plev' .and. index(special, 'catplev') < 1) then
          call interp_z
        else
          ofld = ifld
        end if
        call special_post
        if (badrec) ofld = 1e20

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, r6hrplev) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, r6hrplev) > 0) call close_ofile

    end do

    ! Process table 6hrPlevPt
    write(*, *) 'Process table 6hrPlevPt'
    fnm = p6hrPlevPt
    table = t6hrPlevPt
    linstant = .true.
    lreadplev = .true.
    do n = 1, n6hrPlevPt
      if (skip_variable(n, n6hrPlevPt, d6hrPlevPt)) cycle

      ! Map namelist variables
      ovnm = v6hrPlevPt(ovnmpos, n)
      ivnm = v6hrPlevPt(ivnmpos, n)
      special = v6hrPlevPt(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
        ! Read pressure levels from table
        if (lreadplev .and. zcoord(1:4) == 'plev') then
          call read_gridinfo_plev(trim(zcoord))
        end if
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = taga6hri

      ! Check if input variable is present
      if (len_trim(p6hrPlevPt) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, r6hrPlevPt) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        if (len_trim(p6hrPlevPt) == 0) call scan_files
        if (rec == 0) exit
        tbnds(:, 1) = tbnd
        tval = tbnds(2, 1)

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        if (zcoord(1:4) == 'plev' .and. index(special, 'catplev') < 1) then
          call interp_z
        else
          ofld = ifld
        end if
        call special_post
        if (badrec) ofld = 1e20

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, r6hrPlevPt) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, r6hrPlevPt) > 0) call close_ofile

    end do

    ! Process table 3hr
    write(*, *) 'Process table 3hr'
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
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = taga3hr

      ! Check if input variable is present
      if (len_trim(p3hr) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, r3hr) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        if (len_trim(p3hr) == 0) call scan_files
        if (rec == 0) exit
        tbnds(:, 1) = tbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        ofld = ifld
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, r3hr) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, r3hr) > 0) call close_ofile

    end do

    ! Process table 3hr (instantaneous)
    write(*, *) 'Process table 3hri'
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
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = taga3hri

      ! Check if input variable is present
      if (len_trim(p3hri) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, r3hri) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        if (len_trim(p3hri) == 0) call scan_files
        if (rec == 0) exit
        tbnds(:, 1) = tbnd
        tval = tbnds(2, 1)

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        ofld = ifld
        call special_post
        if (badrec) ofld = 1e20

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, r3hri) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, r3hri) > 0) call close_ofile

    end do

    ! Process table E3hr
    write(*, *) 'Process table E3hr'
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
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = taga3hr

      ! Check if input variable is present
      if (len_trim(pE3hr) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rE3hr) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        if (len_trim(pE3hr) == 0) call scan_files
        if (rec == 0) exit
        tbnds(:, 1) = tbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        ofld = ifld
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, rE3hr) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rE3hr) > 0) call close_ofile

    end do

    ! Process table E3hrPt
    write(*, *) 'Process table E3hrPt'
    fnm = pE3hrPt
    table = tE3hrPt
    linstant = .true.
    do n = 1, nE3hrPt
      if (skip_variable(n, nE3hrPt, dE3hrPt)) cycle

      ! Map namelist variables
      ovnm = vE3hrPt(ovnmpos, n)
      ivnm = vE3hrPt(ivnmpos, n)
      special = vE3hrPt(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = taga3hri

      ! Check if input variable is present
      if (len_trim(pE3hrPt) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rE3hrPt) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        if (len_trim(pE3hrPt) == 0) call scan_files
        if (rec == 0) exit
        tbnds(:, 1) = tbnd
        tval = tbnds(2, 1)

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        ofld = ifld
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, rE3hrPt) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rE3hrPt) > 0) call close_ofile

    end do

    ! Process table CF3hr
    write(*, *) 'Process table CF3hr'
    fnm = pCF3hr
    table = tCF3hr
    linstant = .false.
    do n = 1, nCF3hr
      if (skip_variable(n, nCF3hr, dCF3hr)) cycle

      ! Map namelist variables
      ovnm = vCF3hr(ovnmpos, n)
      ivnm = vCF3hr(ivnmpos, n)
      special = vCF3hr(3, n)
      vunits = ' '
      vpositive = ' '
      vcomment = ' '

      ! Check if 2d and/or 3d fields are to be written
      call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
      !if (zcoord(1:4) == 'plev' .or. zcoord(2:5) == 'leve' .or. &
        !zcoord(2:5) == 'levh') then
        !if (.not. do_3d) cycle
      !else
        !if (.not. do_2d) cycle
      !end if

      ! Choose history file
      itag = taga3hri

      ! Check if input variable is present
      if (len_trim(pCF3hr) == 0) call scan_files(reset=.true.)
      if (.not. var_in_file(fnm, ivnm)) cycle

      ! Prepare output file
      call special_pre

      ! Loop over files
      m = 0
      do
        m = m + 1

        ! Open output file
        if (mod(m - 1, rCF3hr) == 0) call open_ofile

        ! Read variable into buffer (average if necessary)
        rec = 0
        if (len_trim(pCF3hr) == 0) call scan_files
        if (rec == 0) exit
        tbnds(:, 1) = tbnd
        tval = tbnds(2, 1)

        ! Read data
        call read_tslice(rec, badrec, fnm)

        ! Post processing
        if (.not. (index(special, 'ps2pfull') > 0 .or. &
          index(special, 'ps2phalf') > 0)) then
          ofld = ifld
        end if
        call special_post

        ! Write time slice to output file
        call write_tslice

        ! Close output file if max rec has been reached
        if (mod(m, rCF3hr) == 0) call close_ofile

      end do

      ! Close output file if still open
      if (mod(m, rCF3hr) > 0) call close_ofile

    end do

  end subroutine atm2cmor

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

        ! CO2 units
      case ('co2units')
        vunits = '1e-6'

        ! Set correct units for percentage
      case ('percent')
        vunits = '%'

        ! Unit transformation: m -> kg m-2
      case ('kg m-2', 'calcload')
        vunits = 'kg m-2'

        ! Unit transformation: m s-1 -> kg m-2 s-1
      case ('kg m-2 s-1')
        vunits = 'kg m-2 s-1'

        ! Fix micrometers units
      case ('micrometer')
        vunits = 'micrometers'

        ! Fix m units
      case ('m')
        vunits = 'm'

        ! Fix m/s units
      case ('m s-1')
        vunits = 'm s-1'

        ! Fix m-2 units
      case ('m-2')
        vunits = 'm-2'

        ! Set N_AER units
      case ('cm-3')
        vunits = 'cm-3'

        ! Set OH_aft units
      case ('mol mol-1')
        vunits = 'mol mol-1'

        ! Convert degress Celsius to K
      case ('Celsius2Kelvin')
        vunits = 'K'

        ! Convert units from radians2 to m2
      case ('rad2m')
        vunits = 'm2'

        ! Set positive attribute
      case ('positiveup')
        vpositive = 'up'
        write(*, *) 'set positive attribute to ', trim(vpositive)
      case ('positivedo')
        vpositive = 'down'
        write(*, *) 'set positive attribute to ', trim(vpositive)

        ! Write comment for hur and hurs
      case ('hurcomment')
        vcomment = 'field is weighted towards water-ice fraction'

        ! Write comment for evspsbl
      case ('evscomment')
        vcomment = 'field includes dew fall and therefore may ' &
          // 'have negative values'

        ! Times dry mass
      case ('timesmass')
        vunits = 'kg'

        ! Weight with CLDFOC^-1 or FOCHANA^-1
      case ('cldfoc', 'fochana')
        vcomment = &
          'Variable definition deviates from the CMIP5 table ' &
          // 'definition. The variable is weighted toward ' &
          // 'frequency of occurence of warm clouds.'

        ! Weight with CLDFOC^-1 or FOCHANA^-1
      case ('blayer')
        vcomment = &
          'Variable definition deviates from the CMIP5 table ' &
          // 'definition. Values are taken from the lowest model layer.'

      end select
      if (str1 == str2) exit
    end do

  end subroutine special_pre

  ! -----------------------------------------------------------------

  subroutine special_post

    implicit none

    integer :: i, j, k

    str2 = special
    do
      if (index(str2, ';') > 0) then
        str1 = str2(1:index(str2, ';') - 1)
        str2 = str2(index(str2, ';') + 1:)
      else
        str1 = str2
      end if
      select case (str1)

        ! Convert degress Celsius to K
      case ('Celsius2Kelvin')
        ofld = ofld + 273.15

        ! CO2 units
      case ('co2units')
        ofld = ofld * 1e6 * 29. / (12. + 2. * 16.)

        ! Density weighted vertical integration
      case ('calcload')
        call calcload

        ! Convert units from radians2 to m2
      case ('rad2m')
        ofld = ofld * 6.37122e6**2

        ! Unit transformation: m -> kg m-2
      case ('kg m-2')
        ofld = ofld * 1e3

        ! Unit transformation: m s-1 -> kg m-2 s-1
      case ('kg m-2 s-1')
        ofld = ofld * 1e3

        ! Weight with DAYFOC^-1, CLDFOC^-1 or FOCHANA^-1
      case ('dayfoc', 'cldfoc', 'fochana', 'Dfield2')
        do k = 1, kdm
          do j = 1, jdm
            do i = 1, idm
              if (abs(ifld2(i, j, k)) > 1e-15) then
                ofld(i, j, k) = ofld(i, j, k) / ifld2(i, j, k)
              else
                ofld(i, j, k) = 1e20
              end if
            end do
          end do
        end do

        ! Weight with land fraction
      case ('landfrac')
        ofld = ofld * ifld2

        ! Times dry mass (computed offline from PS and Q)
      case ('timesmass')
        ofld = ofld * 5.11253535805483d+18

        ! Invert fraction field
      case ('invert')
        ofld = 1 - ofld

        ! Divide by g
      case ('xginv')
        ofld = ofld / 9.80665

        ! Derive pressure at level mid-points
      case ('ps2pfull')
        do j = 1, jdm
          do i = 1, idm
            do k = 1, ldm
              ofld(i, j, k) = p0 * hyam(k) + ps(i, j) * hybm(k)
            end do
          end do
        end do

        ! Derive pressure at level interfaces
      case ('ps2phalf')
        do j = 1, jdm
          do i = 1, idm
            do k = 1, ldm + 1
              ofld(i, j, k) = p0 * hyai(k) + ps(i, j) * hybi(k)
            end do
          end do
        end do

        ! Set as missing if abs(ofld) is larger then 1.0e+20
      case ('resetfill')
        do k = 1, kdm
          do j = 1, jdm
            do i = 1, idm
              if (abs(ofld(i, j, k)) > 1e20) then
                ofld(i, j, k) = 1e20
              end if
            end do
          end do
        end do

        ! Reverse sign
      case ('Xminus')
        ofld = -ofld

      end select
      if (str1 == str2) exit
    end do

  end subroutine special_post

  ! -----------------------------------------------------------------

  subroutine calcload

    implicit none

    real(kind=8)    :: pt, pb, ps1, fldint
    integer         :: i, j, k

    ! Compute dp and dz
    do j = 1, jdm
      do i = 1, idm
        fldint = 0.
        do k = 1, kdm
          ps1 = ps(i, j)
          pt = p0 * hyai(k) + ps1 * hybi(k)
          pb = p0 * hyai(k + 1) + ps1 * hybi(k + 1)
          fldint = fldint + ofld(i, j, k) * (pb - pt) * ginv
        end do
        ofld(i, j, kdm) = fldint
      end do
    end do

  end subroutine calcload

  ! -----------------------------------------------------------------

  subroutine read_gridinfo_ifile

    implicit none

    integer :: i, j

    ! Open first input file
    call scan_files(reset=.true.)
    status = nf90_open(fnm, nf90_nowrite, ncid)
    call handle_ncerror(status)

    ! Read longitudes
    status = nf90_inq_dimid(ncid, 'lon', dimid)
    call handle_ncerror(status)
    status = nf90_inquire_dimension(ncid, dimid, len=idm)
    call handle_ncerror(status)

    allocate(lon(idm), slon(idm), lon_bnds(2, idm), stat=status)
    if (status /= 0) stop 'cannot ALLOCATE enough memory (1)'

    status = nf90_inq_varid(ncid, 'lon', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, lon)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'slon', rhid)
    if (status == 0) then
      status = nf90_get_var(ncid, rhid, slon)
      call handle_ncerror(status)
    else
      slon(1) = 0.5 * (lon(1) + lon(idm) - 360)
      do i = 2, idm
        slon(i) = 0.5 * (lon(i) + lon(i - 1))
      end do
    end if

    do i = 1, idm - 1
      lon_bnds(1, i) = slon(i)
      lon_bnds(2, i) = slon(i + 1)
    end do
    lon_bnds(1, idm) = slon(idm)
    lon_bnds(2, idm) = slon(1) + 360.

    ! Read latitudes
    status = nf90_inq_dimid(ncid, 'lat', dimid)
    call handle_ncerror(status)
    status = nf90_inquire_dimension(ncid, dimid, len=jdm)
    call handle_ncerror(status)

    allocate(lat(jdm), slat(jdm - 1), lat_bnds(2, jdm), stat=status)
    if (status /= 0) stop 'cannot ALLOCATE enough memory (2)'

    status = nf90_inq_varid(ncid, 'lat', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, lat)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'slat', rhid)
    if (status == 0) then
      status = nf90_get_var(ncid, rhid, slat)
      call handle_ncerror(status)
    else
      do j = 1, jdm - 1
        slat(j) = 0.5 * (lat(j + 1) + lat(j))
      end do
    end if

    lat_bnds(1, 1) = -90.
    lat_bnds(2, 1) = slat(1)
    do j = 2, jdm - 1
      lat_bnds(1, j) = slat(j - 1)
      lat_bnds(2, j) = slat(j)
    end do
    lat_bnds(1, jdm) = slat(jdm - 1)
    lat_bnds(2, jdm) = 90.

    ! Read vertical hybrid coefficients
    status = nf90_inq_dimid(ncid, 'lev', dimid)
    if (status == 0) then
      status = nf90_inquire_dimension(ncid, dimid, len=ldm)
      call handle_ncerror(status)
    else
      write(*, *) 'WARNING: no lev coord found, search for ilev...'
      status = nf90_inq_dimid(ncid, 'ilev', dimid)
      if (status == 0) then
        status = nf90_inquire_dimension(ncid, dimid, len=ldm)
        call handle_ncerror(status)
        ldm = ldm - 1
      else
        call handle_ncerror(status)
      end if
    end if

    allocate(lev(ldm), hyam(ldm), hybm(ldm), ilev(ldm + 1), hyai(ldm + 1), &
      hybi(ldm + 1), iilev(ldm + 2), hyaii(ldm + 2), hybii(ldm + 2), stat=status)
    if (status /= 0) stop 'cannot allocate enough memory'

    status = nf90_inq_varid(ncid, 'lev', rhid)
    if (status == 0) then
      status = nf90_get_var(ncid, rhid, lev)
    end if
    status = nf90_inq_varid(ncid, 'ilev', rhid)
    if (status == 0) then
      status = nf90_get_var(ncid, rhid, ilev)
    end if
    status = nf90_inq_varid(ncid, 'hyam', rhid)
    if (status == 0) then
      status = nf90_get_var(ncid, rhid, hyam)
    else
      write(*, *) 'WARNING: no hyam found...'
    end if
    status = nf90_inq_varid(ncid, 'hybm', rhid)
    if (status == 0) then
      status = nf90_get_var(ncid, rhid, hybm)
    else
      write(*, *) 'WARNING: no hybm found...'
    end if
    status = nf90_inq_varid(ncid, 'ilev', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, ilev)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'hyai', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, hyai)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'hybi', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, hybi)
    call handle_ncerror(status)

    iilev(1) = ilev(1)
    iilev(2:ldm + 1) = lev
    iilev(ldm + 2) = ilev(ldm + 1)
    hyaii(1) = hyai(1)
    hyaii(2:ldm + 1) = hyam
    hyaii(ldm + 2) = hyai(ldm + 1)
    hybii(1) = hybi(1)
    hybii(2:ldm + 1) = hybm
    hybii(ldm + 2) = hybi(ldm + 1)

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

  subroutine read_gridinfo_plev(str)

    use json_module

    implicit none

    character(len=*), intent(in)                    :: str
    integer, parameter                              :: kmax = 200
    real(kind=8)                                    :: kvec(kmax)
    character(len=200)                              :: c200
    logical                                         :: fexists
    type(json_file)                                 :: json
    logical                                         :: found
    integer                                         :: k
    character(len=10), dimension(:), allocatable    :: cval

    inquire(file=trim(tabledir)//trim(coordtable), exist=fexists)
    if (.not. fexists) then
      write(*, *) 'Table ', trim(tabledir)//trim(coordtable), &
        ' does not exist'
      stop
    end if
    call json%initialize()
    call json%load_file(filename=trim(tabledir)//trim(coordtable))
    call json%get('axis_entry.'//trim(str)//'.requested', cval, found)
    call json%destroy()
    pdm = size(cval)
    do k = 1, pdm
      read(cval(k)(1:10), *) kvec(k)
    end do
    deallocate(cval)

    if (allocated(plevs)) deallocate(plevs)
    allocate(plevs(pdm))
    if (status /= 0) stop 'cannot ALLOCATE enough memory (3)'
    plevs = kvec(1:pdm)
    write(*, *) trim(str), ': ', plevs

  end subroutine read_gridinfo_plev

  ! -----------------------------------------------------------------

  subroutine open_ofile(fx)

    implicit none

    logical, optional, intent(in) :: fx
    logical :: fxflag

    real :: fac1, fac2, fac3, fac4, fac5, fac6
    real :: fac7, fac8, fac9, fac10, fac11, fac12
    integer, parameter  :: ndimmax = 10
    integer             :: n, ndims, dimids(ndimmax), dimlens(ndimmax)
    character(len=slenmax) :: ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6, str
    character(len=slenmax) :: ivnm7, ivnm8, ivnm9, ivnm10, ivnm11, ivnm12
    real(kind=8), dimension(:), allocatable :: tmp1d, tmp1d_2
    real(kind=8), dimension(:), allocatable :: tmp2d

    ! Check if output variable should have time coordinate
    fxflag = .false.
    if (present(fx)) then
      if (fx) fxflag = .true.
    end if

    ! Inquire variable units and dimensions in input file
    status = nf90_open(fnm, nf90_nowrite, ncid)
    call handle_ncerror(status)

    call resolve_vnm(slenmax, ivnm, ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6, &
      fac1, fac2, fac3, fac4, fac5, fac6, &
      ivnm7, ivnm8, ivnm9, ivnm10, ivnm11, ivnm12, &
      fac7, fac8, fac9, fac10, fac11, fac12)
    if (verbose) then
      write(*, *) 'Resolve input variable term ', trim(ivnm), ' ='
      if (len_trim(ivnm1) > 0) write(*, *) ' ', trim(ivnm1), '*', fac1
      if (len_trim(ivnm2) > 0) write(*, *) ' + ', trim(ivnm2), '*', fac2
      if (len_trim(ivnm3) > 0) write(*, *) ' + ', trim(ivnm3), '*', fac3
      if (len_trim(ivnm4) > 0) write(*, *) ' + ', trim(ivnm4), '*', fac4
      if (len_trim(ivnm5) > 0) write(*, *) ' + ', trim(ivnm5), '*', fac5
      if (len_trim(ivnm6) > 0) write(*, *) ' + ', trim(ivnm6), '*', fac6
      if (len_trim(ivnm7) > 0) write(*, *) ' + ', trim(ivnm7), '*', fac7
      if (len_trim(ivnm8) > 0) write(*, *) ' + ', trim(ivnm8), '*', fac8
      if (len_trim(ivnm9) > 0) write(*, *) ' + ', trim(ivnm9), '*', fac9
      if (len_trim(ivnm10) > 0) write(*, *) ' + ', trim(ivnm10), '*', fac10
      if (len_trim(ivnm11) > 0) write(*, *) ' + ', trim(ivnm11), '*', fac11
      if (len_trim(ivnm12) > 0) write(*, *) ' + ', trim(ivnm12), '*', fac12
    end if
    if (index(special, 'catplev') >= 1) then
      str = ' '
      if (nint(plevs(1) * 1e-2) < 1000) then
        write(str, '(i3.3)') nint(plevs(1) * 1e-2)
      else
        write(str, '(i4.4)') nint(plevs(1) * 1e-2)
      end if
      if (ivnm1(1:2) == 'Z3') then
        ivnm1 = 'Z'//trim(str)
      else
        ivnm1 = trim(ivnm1)//trim(str)
      end if
    end if
    status = nf90_inq_varid(ncid, trim(ivnm1), rhid)
    if (status /= nf90_noerr) then
      write(*, *) 'cannot find input variable ', trim(ivnm1)
      stop
    end if
    status = nf90_inquire_variable(ncid, rhid, ndims=ndims)
    call handle_ncerror(status)
    status = nf90_inquire_variable(ncid, rhid, dimids=dimids(1:ndims))
    call handle_ncerror(status)
    if (ndims < 3) then
      write(*, *) 'Variable ', trim(ivnm1), ' has too few dimensions'
    end if
    dimlens = 1
    do n = 1, ndims
      status = nf90_inquire_dimension(ncid, dimids(n), len=dimlens(n))
      call handle_ncerror(status)
    end do
    if (dimlens(1) /= idm) then
      write(*, *) 'unexpected first dimension of variable ', &
        trim(ivnm1), ': ', dimlens(1), ' versus idm=', idm
    end if
    if (dimlens(2) /= jdm) then
      write(*, *) 'unexpected second dimension of variable ', &
        trim(ivnm1), ': ', dimlens(2), ' versus jdm=', idm
    end if
    if (ndims > 3) then
      kdm = dimlens(3)
    else
      kdm = 1
    end if
    if (index(special, 'catplev') >= 1) then
      kdm = pdm
    end if
    if (allocated(ifld)) deallocate(ifld, ifld2, ifldacc)
    if (dimlens(1) == 1) then
      allocate(ifld(1, jdm, kdm), ifld2(1, jdm, kdm), &
        ifldacc(1, jdm, kdm), stat=status)
    else
      allocate(ifld(idm, jdm, kdm), ifld2(idm, jdm, kdm), &
        ifldacc(idm, jdm, kdm), stat=status)
    end if
    if (status /= 0) stop 'cannot ALLOCATE enough memory (4)'

    if (len_trim(vunits) == 0) then
      status = nf90_get_att(ncid, rhid, 'units', vunits)
      call handle_ncerror(status)
    end if

    status = nf90_close(ncid)
    call handle_ncerror(status)

    ! Inquire time and vertical dimension of output variable
    if (.not. fxflag) then
      call get_timecoord(trim(tabledir)//trim(table), ovnm, tcoord)
      if (len_trim(tcoord) == 0) tcoord = 'time'
    end if
    call get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)
    if (zcoord(2:4) == 'lev' .or. index(special, 'calcload') > 0) then
      if (.not. allocated(ps)) then
        allocate(ps(idm, jdm), stat=status)
        if (status /= 0) stop 'cannot allocate enough memory'
      end if
    end if

    lreqphis = .false.
    if (trim(ivnm1) == 'Z3' .or. trim(ivnm1) == 'T') then
      lreqphis = .true.
      if (.not. allocated(phis)) then
        allocate(phis(idm, jdm), tbot(idm, jdm), stat=status)
        if (status /= 0) stop 'cannot allocate enough memory'
      end if
    end if

    lreqsst = .false.
    if (trim(ovnm) == 'tslsi') then
      lreqsst = .true.
      if (.not. allocated(sst)) then
        allocate(sst(idm, jdm), stat=status)
        if (status /= 0) stop 'cannot allocate enough memory'
      end if
    end if

    ! Call CMOR setup
    if (verbose) then
      if (createsubdirs) then
        error_flag = cmor_setup(inpath=trim(tabledir), &
          netcdf_file_action=CMOR_REPLACE_4, &
          set_verbosity=CMOR_NORMAL, &
          create_subdirectories=1)
      else
        error_flag = cmor_setup(inpath=trim(tabledir), &
          netcdf_file_action=CMOR_REPLACE_4, &
          set_verbosity=CMOR_NORMAL, &
          create_subdirectories=0)
      end if
    else
      if (createsubdirs) then
        error_flag = cmor_setup(inpath=trim(tabledir), &
          netcdf_file_action=CMOR_REPLACE_4, &
          set_verbosity=CMOR_QUIET, &
          create_subdirectories=1)
      else
        error_flag = cmor_setup(inpath=trim(tabledir), &
          netcdf_file_action=CMOR_REPLACE_4, &
          set_verbosity=CMOR_QUIET, &
          create_subdirectories=0)
      end if
    end if
    if (error_flag /= 0) stop 'Problem setting up CMOR'

    ! Define output dataset
    if (ndims == 1) then
      call json_write_attributes('global mean or integral', 'gm', &
        atmgrid_resolution, ovnm)
    else
      call json_write_attributes(atmgrid, atmgrid_label, &
        atmgrid_resolution, ovnm)
    end if
    error_flag = cmor_dataset_json(json_file_attributes)

    ! Define horizontal axes
    write(*, *) 'define horizontal axes'
    write(*, *) 'tabledir:', trim(tabledir)
    write(*, *) 'table:', trim(table)

    iaxid = cmor_axis( &
      table=trim(tabledir)//trim(table), &
      table_entry='longitude', &
      units='degrees_east', &
      length=idm, &
      coord_vals=lon, &
      cell_bounds=lon_bnds)

    write(*, *) 'line 2458'
    jaxid = cmor_axis( &
      table=trim(tabledir)//trim(table), &
      table_entry='latitude', &
      units='degrees_north', &
      length=jdm, &
      coord_vals=lat, &
      cell_bounds=lat_bnds)

    ! Define time axis
    if (.not. fxflag) then
      taxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry=trim(tcoord), &
        units=trim(calunits), &
        length=1)
    end if

    ! Define vertical axis
    write(*, *) 'define vertical axes'
    str = 'ps'
    if (tcoord(1:5) == 'time1') str = 'ps1'
    if (zcoord(1:4) == 'plev') then
      kaxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry=trim(zcoord), &
        units='Pa', &
        length=pdm, &
        coord_vals=plevs)
    else if (trim(zcoord) == 'alev1') then
      allocate(tmp1d(kdm), tmp1d_2(kdm + 1))
      tmp1d(:) = lev(kdm) * 1d-3
      tmp1d_2(1) = dble(ilev(kdm + 1)) * 1.d-3
      tmp1d_2(2) = dble(ilev(kdm)) * 1.d-3
      kaxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry='standard_hybrid_sigma', &
        units='1', &
        length=1, &
        coord_vals=tmp1d, &
        cell_bounds=tmp1d_2)
      deallocate(tmp1d, tmp1d_2)
      error_flag = cmor_zfactor( &
        zaxis_id=kaxid, &
        zfactor_name='p0', &
        units='Pa', &
        zfactor_values=p0)
      error_flag = cmor_zfactor( &
        zaxis_id=kaxid, &
        zfactor_name='b', &
        axis_ids=(/kaxid/), &
        zfactor_values=(/hybm(kdm)/), &
        zfactor_bounds=(/hybi(kdm + 1), hybi(kdm)/))
      error_flag = cmor_zfactor( &
        zaxis_id=kaxid, &
        zfactor_name='a', &
        axis_ids=(/kaxid/), &
        zfactor_values=(/hyam(kdm)/), &
        zfactor_bounds=(/hyai(kdm + 1), hyai(kdm)/))
      zfacid = cmor_zfactor( &
        zaxis_id=kaxid, &
        zfactor_name=trim(str), &
        axis_ids=(/iaxid, jaxid, taxid/), &
        units='Pa')
    else if (trim(zcoord) == 'alevel') then
      write(*, *) 'l 2532'
      allocate(tmp1d(size(lev)), tmp1d_2(size(ilev)))
      tmp1d(:) = 1.d-3 * lev(:)
      tmp1d_2(:) = 1.d-3 * ilev(:)
      write(*, *) 'tmp1d:', tmp1d
      write(*, *) 'tmp1d_2:', tmp1d_2
      kaxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry='standard_hybrid_sigma', &
        units='1', &
        length=ldm, &
        coord_vals=tmp1d, &
        cell_bounds=tmp1d_2)
      deallocate(tmp1d, tmp1d_2)
      error_flag = cmor_zfactor( &
        zaxis_id=kaxid, &
        zfactor_name='p0', &
        units='Pa', &
        zfactor_values=p0)
      error_flag = cmor_zfactor( &
        zaxis_id=kaxid, &
        zfactor_name='b', &
        axis_ids=(/kaxid/), &
        zfactor_values=hybm, &
        zfactor_bounds=hybi)
      error_flag = cmor_zfactor( &
        zaxis_id=kaxid, &
        zfactor_name='a', &
        axis_ids=(/kaxid/), &
        zfactor_values=hyam, &
        zfactor_bounds=hyai)
      zfacid = cmor_zfactor( &
        zaxis_id=kaxid, &
        zfactor_name=trim(str), &
        axis_ids=(/iaxid, jaxid, taxid/), &
        units='Pa')
    else if (trim(zcoord) == 'alevhalf') then
      allocate(tmp1d(size(ilev)), tmp1d_2(size(iilev)))
      tmp1d(:) = ilev(:) * 1.d-3
      tmp1d_2(:) = iilev(:) * 1.d-3
      kaxid = cmor_axis( &
        table=trim(tabledir)//trim(table), &
        table_entry='standard_hybrid_sigma_half', &
        units='1', &
        length=ldm + 1, &
        coord_vals=tmp1d)
      deallocate(tmp1d, tmp1d_2)
      error_flag = cmor_zfactor( &
        zaxis_id=kaxid, &
        zfactor_name='p0', &
        units='Pa', &
        zfactor_values=p0)
      error_flag = cmor_zfactor( &
        zaxis_id=kaxid, &
        zfactor_name='b_half', &
        axis_ids=(/kaxid/), &
        zfactor_values=hybi)
      error_flag = cmor_zfactor( &
        zaxis_id=kaxid, &
        zfactor_name='a_half', &
        axis_ids=(/kaxid/), &
        zfactor_values=hyai)
      zfacid = cmor_zfactor( &
        zaxis_id=kaxid, &
        zfactor_name=trim(str), &
        axis_ids=(/iaxid, jaxid, taxid/), &
        units='Pa')
    end if

    ! Define output variable
    write(*, *) 'define outpout variable'
    write(*, *) 'zcoord(1:4):', zcoord(1:4)
    if (fxflag) then
      varid = cmor_variable( &
        table=trim(tabledir)//trim(table), &
        table_entry=trim(ovnm), &
        units=trim(vunits), &
        axis_ids=(/iaxid, jaxid/), &
        original_name=trim(ivnm), &
        missing_value=1e20, &
        comment=trim(vcomment))
    else
      if (zcoord(1:4) == 'plev' .or. zcoord(1:4) == 'alev') then
        write(*, *) 'dimlens(1):', dimlens(1)
        if (dimlens(1) == 1) then
          varid = cmor_variable( &
            table=trim(tabledir)//trim(table), &
            table_entry=trim(ovnm), &
            units=trim(vunits), &
            axis_ids=(/jaxid, kaxid, taxid/), &
            original_name=trim(ivnm), &
            missing_value=1e20, &
            comment=trim(vcomment), &
            positive=trim(vpositive))
        else
          write(*, *) 'l2646'
          varid = cmor_variable( &
            table=trim(tabledir)//trim(table), &
            table_entry=trim(ovnm), &
            units=trim(vunits), &
            axis_ids=(/iaxid, jaxid, kaxid, taxid/), &
            original_name=trim(ivnm), &
            missing_value=1e20, &
            comment=trim(vcomment), &
            positive=trim(vpositive))
          write(*, *) 'l2656'
        end if
      else if (ndims == 1) then
        varid = cmor_variable( &
          table=trim(tabledir)//trim(table), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/taxid/), &
          original_name=trim(ivnm), &
          missing_value=1e20, &
          comment=trim(vcomment), &
          positive=trim(vpositive))
      else
        if (lreqsst) then
          varid = cmor_variable( &
            table=trim(tabledir)//trim(table), &
            table_entry=trim(ovnm), &
            units=trim(vunits), &
            axis_ids=(/iaxid, jaxid, taxid/), &
            original_name=trim(ivnm), &
            missing_value=1e20, &
            comment=trim(vcomment), &
            positive=trim(vpositive))
        else
          varid = cmor_variable( &
            table=trim(tabledir)//trim(table), &
            table_entry=trim(ovnm), &
            units=trim(vunits), &
            axis_ids=(/iaxid, jaxid, taxid/), &
            original_name=trim(ivnm), &
            missing_value=1e20, &
            comment=trim(vcomment), &
            positive=trim(vpositive))
        end if
      end if
    end if

#ifdef DEFLATE
    error_flag = cmor_set_deflate(varid, 1, 1, 5)
#endif

    ! Allocate memory for output variable
    if (allocated(ofld)) deallocate(ofld)
    if (zcoord(1:4) == 'plev') then
      if (dimlens(1) == 1) then
        allocate(ofld(1, jdm, pdm), stat=status)
      else
        allocate(ofld(idm, jdm, pdm), stat=status)
      end if
    else if (index(special, 'ps2pfull') > 0) then
      allocate(ofld(idm, jdm, ldm), stat=status)
    else if (index(special, 'ps2phalf') > 0) then
      allocate(ofld(idm, jdm, ldm + 1), stat=status)
    else
      allocate(ofld(idm, jdm, kdm), stat=status)
    end if
    if (status /= 0) stop 'cannot allocate enough memory (0)'

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
    status = nf90_inq_varid(ncid, trim(ivnm1), rhid)
    if (status /= nf90_noerr) then
      write(*, *) 'cannot find input variable ', trim(ivnm1)
      stop
    end if
    status = nf90_get_var(ncid, rhid, ifld)
    call handle_ncerror(status)
    if (fac1 /= 1) then
      ifld = ifld * fac1
    end if

    if (len_trim(ivnm2) > 0) then
      status = nf90_inq_varid(ncid, trim(ivnm2), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm2)
        stop
      end if
      status = nf90_get_var(ncid, rhid, ifld2)
      call handle_ncerror(status)
      ifld = ifld + ifld2 * fac2
    end if

    if (len_trim(ivnm3) > 0) then
      status = nf90_inq_varid(ncid, trim(ivnm3), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm3)
        stop
      end if
      status = nf90_get_var(ncid, rhid, ifld2)
      call handle_ncerror(status)
      ifld = ifld + ifld2 * fac3
    end if

    if (len_trim(ivnm4) > 0) then
      status = nf90_inq_varid(ncid, trim(ivnm4), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm4)
        stop
      end if
      status = nf90_get_var(ncid, rhid, ifld2)
      call handle_ncerror(status)
      ifld = ifld + ifld2 * fac4
    end if

    if (len_trim(ivnm5) > 0) then
      status = nf90_inq_varid(ncid, trim(ivnm5), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm5)
        stop
      end if
      status = nf90_get_var(ncid, rhid, ifld2)
      call handle_ncerror(status)
      ifld = ifld + ifld2 * fac5
    end if

    if (len_trim(ivnm6) > 0) then
      status = nf90_inq_varid(ncid, trim(ivnm6), rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable ', trim(ivnm6)
        stop
      end if
      status = nf90_get_var(ncid, rhid, ifld2)
      call handle_ncerror(status)
      ifld = ifld + ifld2 * fac6
    end if

    status = nf90_close(ncid)
    call handle_ncerror(status)

    if (index(special, 'landfrac') > 0) then
      ifld2 = 1.
      status = nf90_open(fnm2, nf90_nowrite, ncid)
      status = nf90_inq_varid(ncid, 'LANDFRAC', rhid)
      status = nf90_get_var(ncid, rhid, ifld2)
      status = nf90_close(ncid)
      call handle_ncerror(status)
    end if

  end subroutine read_field

  ! -----------------------------------------------------------------

  subroutine read_tslice(rec, badrec, fname)

    implicit none

    real                :: fac1, fac2, fac3, fac4, fac5, fac6
    real                :: fac7, fac8, fac9, fac10, fac11, fac12
    integer, intent(in) :: rec
    logical, intent(out):: badrec
    character(len=*), intent(in), optional :: fname
    integer, save           :: fid, fidphis
    character(len=slenmax)  :: ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6, str
    character(len=slenmax)  :: ivnm7, ivnm8, ivnm9, ivnm10, ivnm11, ivnm12
    integer :: i, j

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
      status = nf90_inq_varid(fid, 'time_bnds', rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find time_bnds variable'
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
    end if

    ! Read data
    call resolve_vnm(slenmax, ivnm, ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6, &
      fac1, fac2, fac3, fac4, fac5, fac6, &
      ivnm7, ivnm8, ivnm9, ivnm10, ivnm11, ivnm12, &
      fac7, fac8, fac9, fac10, fac11, fac12)
    if (index(special, 'catplev') < 1) then
      status = nf90_inq_varid(fid, trim(ivnm1), rhid)
    end if
    if (index(ivnm1, 'vmr') > 0) then
      status = nf90_get_var(fid, rhid, ifld, (/rec/), (/1/))
    else
      if (index(special, 'catplev') < 1) then
        if (status /= nf90_noerr) then
          write(*, *) 'cannot find input variable ', trim(ivnm1)
          stop
        end if
        if (.not. readdummy) then
          if (kdm == 1) then
            status = nf90_get_var(fid, rhid, ifld, (/1, 1, rec/), &
              (/idm, jdm, 1/))
          else
            if (index(table, 'EmonZ') > 0) then
              status = nf90_get_var(fid, rhid, ifld(1,:,:), &
                (/1, 1, 1, rec/), &
                (/1, jdm, kdm, 1/))
            else
              status = nf90_get_var(fid, rhid, ifld, (/1, 1, 1, rec/), &
                (/idm, jdm, kdm, 1/))
            end if
          end if
          call handle_ncerror(status)
        else
          ifld = 0.
        end if
        if (fac1 /= 1) then
          ifld = ifld * fac1
        end if
      else
        do i = 1, kdm
          str = ' '
          if (nint(plevs(i) * 1e-2) < 1000) then
            write(str, '(i3.3)') nint(plevs(i) * 1e-2)
          else
            write(str, '(i4.4)') nint(plevs(i) * 1e-2)
          end if
          if (ivnm1(1:2) == 'Z3') then
            status = nf90_inq_varid(fid, 'Z'//trim(str), rhid)
          else
            status = nf90_inq_varid(fid, trim(ivnm1)//trim(str), rhid)
          end if
          if (status /= nf90_noerr) then
            write(*, *) 'cannot find input variable ', &
              trim(ivnm1)//trim(str)
            stop
          end if
          if (.not. readdummy) then
            status = nf90_get_var(fid, rhid, ps, (/1, 1, rec/), &
              (/idm, jdm, 1/))
            ifld(:,:, i) = ps
          else
            ifld = 0.
          end if
        end do
      end if

      if (len_trim(ivnm2) > 0) then
        status = nf90_inq_varid(fid, trim(ivnm2), rhid)
        if (status /= nf90_noerr) then
          write(*, *) 'cannot find input variable ', trim(ivnm2)
          stop
        end if
        if (kdm == 1) then
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, rec/), &
            (/idm, jdm, 1/))
        else
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, 1, rec/), &
            (/idm, jdm, kdm, 1/))
        end if
        call handle_ncerror(status)
        if (index(special, 'Dfield2') > 0) then
          return
        else
          ifld = ifld + ifld2 * fac2
        end if
      end if

      if (len_trim(ivnm3) > 0) then
        status = nf90_inq_varid(fid, trim(ivnm3), rhid)
        if (status /= nf90_noerr) then
          write(*, *) 'cannot find input variable ', trim(ivnm3)
          stop
        end if
        if (kdm == 1) then
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, rec/), &
            (/idm, jdm, 1/))
        else
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, 1, rec/), &
            (/idm, jdm, kdm, 1/))
        end if
        call handle_ncerror(status)
        ifld = ifld + ifld2 * fac3
      end if

      if (len_trim(ivnm4) > 0) then
        status = nf90_inq_varid(fid, trim(ivnm4), rhid)
        if (status /= nf90_noerr) then
          write(*, *) 'cannot find input variable ', trim(ivnm4)
          stop
        end if
        if (kdm == 1) then
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, rec/), &
            (/idm, jdm, 1/))
        else
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, 1, rec/), &
            (/idm, jdm, kdm, 1/))
        end if
        call handle_ncerror(status)
        ifld = ifld + ifld2 * fac4
      end if

      if (len_trim(ivnm5) > 0) then
        status = nf90_inq_varid(fid, trim(ivnm5), rhid)
        if (status /= nf90_noerr) then
          write(*, *) 'cannot find input variable ', trim(ivnm5)
          stop
        end if
        if (kdm == 1) then
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, rec/), &
            (/idm, jdm, 1/))
        else
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, 1, rec/), &
            (/idm, jdm, kdm, 1/))
        end if
        call handle_ncerror(status)
        ifld = ifld + ifld2 * fac5
      end if

      if (len_trim(ivnm6) > 0) then
        status = nf90_inq_varid(fid, trim(ivnm6), rhid)
        if (status /= nf90_noerr) then
          write(*, *) 'cannot find input variable ', trim(ivnm6)
          stop
        end if
        if (kdm == 1) then
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, rec/), &
            (/idm, jdm, 1/))
        else
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, 1, rec/), &
            (/idm, jdm, kdm, 1/))
        end if
        call handle_ncerror(status)
        ifld = ifld + ifld2 * fac6
      end if

      if (len_trim(ivnm7) > 0) then
        status = nf90_inq_varid(fid, trim(ivnm7), rhid)
        if (status /= nf90_noerr) then
          write(*, *) 'cannot find input variable ', trim(ivnm7)
          stop
        end if
        if (kdm == 1) then
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, rec/), &
            (/idm, jdm, 1/))
        else
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, 1, rec/), &
            (/idm, jdm, kdm, 1/))
        end if
        call handle_ncerror(status)
        ifld = ifld + ifld2 * fac7
      end if

      if (len_trim(ivnm8) > 0) then
        status = nf90_inq_varid(fid, trim(ivnm8), rhid)
        if (status /= nf90_noerr) then
          write(*, *) 'cannot find input variable ', trim(ivnm8)
          stop
        end if
        if (kdm == 1) then
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, rec/), &
            (/idm, jdm, 1/))
        else
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, 1, rec/), &
            (/idm, jdm, kdm, 1/))
        end if
        call handle_ncerror(status)
        ifld = ifld + ifld2 * fac8
      end if

      if (len_trim(ivnm9) > 0) then
        status = nf90_inq_varid(fid, trim(ivnm9), rhid)
        if (status /= nf90_noerr) then
          write(*, *) 'cannot find input variable ', trim(ivnm9)
          stop
        end if
        if (kdm == 1) then
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, rec/), &
            (/idm, jdm, 1/))
        else
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, 1, rec/), &
            (/idm, jdm, kdm, 1/))
        end if
        call handle_ncerror(status)
        ifld = ifld + ifld2 * fac9
      end if

      if (len_trim(ivnm10) > 0) then
        status = nf90_inq_varid(fid, trim(ivnm10), rhid)
        if (status /= nf90_noerr) then
          write(*, *) 'cannot find input variable ', trim(ivnm10)
          stop
        end if
        if (kdm == 1) then
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, rec/), &
            (/idm, jdm, 1/))
        else
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, 1, rec/), &
            (/idm, jdm, kdm, 1/))
        end if
        call handle_ncerror(status)
        ifld = ifld + ifld2 * fac10
      end if

      if (len_trim(ivnm11) > 0) then
        status = nf90_inq_varid(fid, trim(ivnm11), rhid)
        if (status /= nf90_noerr) then
          write(*, *) 'cannot find input variable ', trim(ivnm11)
          stop
        end if
        if (kdm == 1) then
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, rec/), &
            (/idm, jdm, 1/))
        else
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, 1, rec/), &
            (/idm, jdm, kdm, 1/))
        end if
        call handle_ncerror(status)
        ifld = ifld + ifld2 * fac11
      end if

      if (len_trim(ivnm12) > 0) then
        status = nf90_inq_varid(fid, trim(ivnm12), rhid)
        if (status /= nf90_noerr) then
          write(*, *) 'cannot find input variable ', trim(ivnm12)
          stop
        end if
        if (kdm == 1) then
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, rec/), &
            (/idm, jdm, 1/))
        else
          status = nf90_get_var(fid, rhid, ifld2, (/1, 1, 1, rec/), &
            (/idm, jdm, kdm, 1/))
        end if
        call handle_ncerror(status)
        ifld = ifld + ifld2 * fac12
      end if
    end if

    ! Read auxillary data
    if (zcoord(2:4) == 'lev' .or. index(special, 'calcload') > 0) then
      if (.not. readdummy) then
        status = nf90_inq_varid(fid, 'PS', rhid)
        if (status /= 0) then
          ps = 0. ! for 6hr (cam.h2), PS are not available
        else
          call handle_ncerror(status)
          status = nf90_get_var(fid, rhid, ps, (/1, 1, rec/), (/idm, jdm, 1/))
          call handle_ncerror(status)
        end if
      end if
    end if

    if (index(special, 'dayfoc') > 0) then
      ifld2 = 1.
      status = nf90_inq_varid(fid, 'DAYFOC', rhid)
      call handle_ncerror(status)
      status = nf90_get_var(fid, rhid, ifld2, (/1, 1, rec/), (/idm, jdm, 1/))
      call handle_ncerror(status)
    else if (index(special, 'fochana') > 0) then
      ifld2 = 1.
      status = nf90_inq_varid(fid, 'FOCHANA', rhid)
      call handle_ncerror(status)
      status = nf90_get_var(fid, rhid, ifld2, (/1, 1, rec/), (/idm, jdm, 1/))
      call handle_ncerror(status)
    else if (index(special, 'cldfoc') > 0) then
      ifld2 = 1.
      status = nf90_inq_varid(fid, 'CLDFOC', rhid)
      call handle_ncerror(status)
      status = nf90_get_var(fid, rhid, ifld2, (/1, 1, 1, rec/), &
        (/idm, jdm, kdm, 1/))
      call handle_ncerror(status)
    end if

    if (index(special, 'landfrac') > 0) then
      ifld2 = 1.
      status = nf90_inq_varid(fid, 'LANDFRAC', rhid)
      call handle_ncerror(status)
      status = nf90_get_var(fid, rhid, ifld2)
      call handle_ncerror(status)
    end if

    if (lreqphis) then
      if (.not. readdummy) then
        status = nf90_inq_varid(fid, 'PHIS', rhid)
        if (status /= nf90_noerr) then
          status = nf90_open(trim(griddata)//trim(atmgridfile), &
            nf90_nowrite, fidphis)
          call handle_ncerror(status)
          status = nf90_inq_varid(fidphis, 'PHIS', rhid)
          call handle_ncerror(status)
          status = nf90_get_var(fidphis, rhid, phis, (/1, 1/), (/idm, jdm/))
          call handle_ncerror(status)
          status = nf90_close(fidphis)
        else
          status = nf90_get_var(fid, rhid, phis, (/1, 1, rec/), &
            (/idm, jdm, 1/))
        end if
        call handle_ncerror(status)
        status = nf90_inq_varid(fid, 'T', rhid)
        call handle_ncerror(status)
        status = nf90_get_var(fid, rhid, tbot, (/1, 1, ldm, rec/), &
          (/idm, jdm, 1, 1/))
        call handle_ncerror(status)
      else
        phis = 0.
      end if
    end if

    ! Mask SSTs
    if (lreqsst) then
      if (.not. readdummy) then
        status = nf90_inq_varid(fid, 'SST', rhid)
        call handle_ncerror(status)
        status = nf90_get_var(fid, rhid, sst, (/1, 1, rec/), (/idm, jdm, 1/))
        call handle_ncerror(status)
      else
        sst = 0.
      end if
      do j = 1, jdm
        do i = 1, idm
          if (ifld(i, j, 1) == sst(i, j)) ifld(i, j, 1) = 1e20
        end do
      end do
    end if

    status = nf90_close(fid)
    call handle_ncerror(status)

  end subroutine read_tslice

  ! -----------------------------------------------------------------

  subroutine write_field

    implicit none

    ! Store variable
    write(*, *) 'l3245'
    if (zcoord(1:4) == 'plev' .or. trim(zcoord) == 'alevel') then
      error_flag = cmor_write( &
        var_id=varid, &
        data=ofld)
      write(*, *) 'l3250'
    else
      error_flag = cmor_write( &
        var_id=varid, &
        data=reshape(ofld, (/idm, jdm/)))
    end if

    ! Store auxillary data
    if (zcoord(1:4) == 'alev') then
      if (trim(tcoord) /= 'time1') then
        error_flag = cmor_write( &
          var_id=zfacid, &
          data=ps, &
          ntimes_passed=1, &
          time_vals=tval, &
          time_bnds=tbnds, &
          store_with=varid)
      else
        error_flag = cmor_write( &
          var_id=zfacid, &
          data=ps, &
          ntimes_passed=1, &
          time_vals=tval, &
          store_with=varid)
      end if
    end if

  end subroutine write_field

  ! -----------------------------------------------------------------

  subroutine write_tslice

    implicit none

    if (dry_run) return

    ! Store variable
    if (trim(tcoord) /= 'time1') then
      if (trim(zcoord) == 'alev1' .or. index(special, 'blayer') > 0 .or. &
        index(special, 'calcload') > 0) then
        error_flag = cmor_write( &
          var_id=varid, &
          data=reshape(ofld(:,:, kdm:kdm), (/idm, jdm/)), &
          ntimes_passed=1, &
          time_vals=tval, &
          time_bnds=tbnds)
      else if (index(ivnm, 'vmr') > 0) then
        error_flag = cmor_write( &
          var_id=varid, &
          data=reshape(ofld(1:1, 1:1, 1:1), (/1/)), &
          ntimes_passed=1, &
          time_vals=tval, &
          time_bnds=tbnds)
      else if (index(table, 'EmonZ') > 0) then
        error_flag = cmor_write( &
          var_id=varid, &
          data=reshape(ofld(1,:,:), (/jdm, pdm/)), &
          ntimes_passed=1, &
          time_vals=tval, &
          time_bnds=tbnds)
      else
        error_flag = cmor_write( &
          var_id=varid, &
          data=ofld, &
          ntimes_passed=1, &
          time_vals=tval, &
          time_bnds=tbnds)
      end if
    else
      if (trim(zcoord) == 'alev1' .or. index(special, 'blayer') > 0 .or. &
        index(special, 'calcload') > 0) then
        error_flag = cmor_write( &
          var_id=varid, &
          data=reshape(ofld(:,:, kdm:kdm), (/idm, jdm/)), &
          ntimes_passed=1, &
          time_vals=tval)
      else
        error_flag = cmor_write( &
          var_id=varid, &
          data=ofld, &
          ntimes_passed=1, &
          time_vals=tval)
      end if
    end if

    ! Store auxillary data
    if (zcoord(1:4) == 'alev') then
      if (trim(tcoord) /= 'time1') then
        error_flag = cmor_write( &
          var_id=zfacid, &
          data=ps, &
          ntimes_passed=1, &
          time_vals=tval, &
          time_bnds=tbnds, &
          store_with=varid)
      else
        error_flag = cmor_write( &
          var_id=zfacid, &
          data=ps, &
          ntimes_passed=1, &
          time_vals=tval, &
          store_with=varid)
      end if
    end if

  end subroutine write_tslice

  ! -----------------------------------------------------------------

  subroutine interp_z

    implicit none

    real(kind=8), save :: p0_hpa, missing
    real(kind=8), save, allocatable, dimension(:) :: tmp1d, plevs_hpa

    if (.not. allocated(tmp1d)) then
      allocate(tmp1d(ldm + 1), stat=status)
      if (status /= 0) stop 'cannot allocate enough memory'
      missing = 1.e20
      p0_hpa = p0 * 1e-2
    end if
    if (allocated(plevs_hpa)) deallocate(plevs_hpa)
    allocate(plevs_hpa(pdm), stat=status)
    if (status /= 0) stop 'cannot allocate enough memory'
    plevs_hpa = plevs * 1e-2

    if (.not. plevdummy) then
      if (trim(ivnm) == 'T') then
        call vinth2pecmwf(ifld, ofld, hyam, hybm, p0_hpa, tmp1d, plevs_hpa, &
          1, ps, missing, 1, idm, jdm, ldm, ldm, pdm, 1, tbot, phis)
      else if (trim(ivnm) == 'Z3') then
        call vinth2pecmwf(ifld, ofld, hyam, hybm, p0_hpa, tmp1d, plevs_hpa, &
          1, ps, missing, 1, idm, jdm, ldm, ldm, pdm, -1, tbot, phis)
      else
        if (index(table, 'EmonZ') > 0) then
          call vinth2pecmwf(ifld, ofld, ilev / 1000., hybi * 0., p0_hpa, tmp1d, &
            plevs_hpa, 1, ps, missing, 0, 1, jdm, ldm + 1, ldm + 1, pdm, 0, tbot, phis)
        else
          call vinth2pecmwf(ifld, ofld, hyam, hybm, p0_hpa, tmp1d, plevs_hpa, &
            1, ps, missing, 1, idm, jdm, ldm, ldm, pdm, 0, tbot, phis)
        end if
      end if
    else
      ofld = 0.
    end if

  end subroutine interp_z

end module m_modelsatm
