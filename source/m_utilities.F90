module m_utilities

  use netcdf
  use m_namelists, only: itag, verbose, ibasedir, slenmax, r8, &
    casename, fnm, year1, month1, yearn, monthn, exprefyear, rec, tval, &
    tbnd, mbnd, year, month, forcefilescan, funit, scanallfiles, membertag

  implicit none

contains


  ! -----------------------------------------------------------------

  logical function var_in_file(fnm, vnm)

    implicit none

    character(len=*), intent(in) :: fnm, vnm

    integer :: ncid, rhid, status, n

    var_in_file = .false.
    if (len_trim(fnm) > 0) then
      status = nf90_open(trim(fnm), nf90_nowrite, ncid)
      if (status /= nf90_noerr) then
        if (verbose) write(*, *) &
          'WARNING: error open file ', trim(fnm), &
          '. Will skip respective output group.'
        return
      end if
    else
      if (verbose) write(*, *) &
        'WARNING: no file found for tag=', trim(itag), &
        '. Will skip respective output group.'
      return
    end if

    var_in_file = .true.
    
    status = nf90_inq_varid(ncid, trim(vnm), rhid)
    if (status /= nf90_noerr) then
      var_in_file = .false.
      if (verbose) write(*, *) &
        'skipping variable. ' // trim(vnm) // ' not in ' // trim(fnm)
    end if

    status = nf90_close(ncid)
    call handle_ncerror(status)

  end function var_in_file

  ! -----------------------------------------------------------------

  integer function get_nrec(fnm)

    implicit none

    character(len=*), intent(in) :: fnm

    integer :: ncid, dimid, status

    status = nf90_open(trim(fnm), nf90_nowrite, ncid)
    call handle_ncerror(status)
    status = nf90_inq_dimid(ncid, 'time', dimid)
    call handle_ncerror(status)
    status = nf90_inquire_dimension(ncid, dimid, len=get_nrec)
    call handle_ncerror(status)
    status = nf90_close(ncid)
    call handle_ncerror(status)

  end function get_nrec

  ! -----------------------------------------------------------------

  subroutine handle_ncerror(status)

    implicit none

    integer, intent(in) :: status

    if (status /= nf90_noerr) then
      write(*, *) trim(nf90_strerror(status))
      stop
    end if

  end subroutine handle_ncerror

  ! -----------------------------------------------------------------

  subroutine scan_files(reset)

    implicit none

    logical, intent(in), optional :: reset
    logical :: isloop, isflag

    isloop = .true.

    do while (isloop)
      if (present(reset)) then
        call get_file_info(ibasedir, casename, itag, fnm, year1, month1, &
          yearn, monthn, exprefyear, reset, rec, tval(1), tbnd, &
          mbnd, year, month)
      else
        call get_file_info(ibasedir, casename, itag, fnm, year1, month1, &
          yearn, monthn, exprefyear, .false., rec, tval(1), tbnd, &
          mbnd, year, month)
      end if

      isloop = .false.
    end do

  end subroutine scan_files

  ! -----------------------------------------------------------------

  subroutine get_file_info(idir, cnam, ftag, fnam, y1, m1, y2, m2, yr, &
      lreset, irec, tval, tbnd, mbnd, year, month)

    use netcdf
    implicit none

    character(len=1024), intent(in)     :: idir, cnam, ftag
    character(len=1024), intent(out)    :: fnam
    integer, intent(in)     :: y1, m1, y2, m2, yr
    logical, intent(in)     :: lreset
    integer, intent(out)    :: irec
    real(r8), intent(out)   :: tval, tbnd(2), mbnd(2)
    integer, intent(out)        :: year, month

    character(len=1024) :: fpre, str1, str2
    integer, save       :: y1old = 0, m1old = 0, y2old = 0, m2old = 0
    integer             :: idx, yyyy, yref, mref, dref

    character(len=1024), save   :: fpreold = 'xxx'
    character(len=1024)         :: units, calendar
    integer, save       :: fstat, n, nrec, unlimdimid, ncid, rhid, rhid2, &
                           status, sstartend(2,10), nskip, irecold, toff, ndays
    real(r8), save  :: dnumlo, dnumhi, tbndold(2) = (/-999999., -999999./)
    logical, save       :: ltimebnds, ldone

    ! check whether to start file list read from beginning or to resume
    idx = index(ftag, '.')
    str1 = ftag(1:idx - 1)
    str2 = ftag(idx + 1:)
    if (len_trim(membertag) > 0) then
      fpre = trim(cnam) // '.' // trim(str1) // '_' // trim(membertag) // '.' // trim(str2) // '.'
    else
      fpre = trim(cnam) // '.' // trim(ftag) // '.'
    end if

    if (y1old /= y1 .or. y2old /= y2 .or. m1old /= m1 .or. &
        m2old /= m2 .or. index(fpreold, trim(fpre)) <= 0 .or. &
        tbndold(2) <= -999998. .or. lreset) then
      nskip = 0
      irecold = 0
    end if

    if (lreset) then
      y1old = 0
      y2old = 0
      m1old = 0
      m2old = 0
      fpreold = ' '
    else
      y1old = y1
      y2old = y2
      m1old = m1
      m2old = m2
      fpreold = fpre
    end if

    ! loop over files
    open(funit, file=trim('filelist_' // cnam) // trim(membertag), form='formatted')

    ! skip all files that have already been scanned
    do n = 1, nskip
      read(funit, '(a1024)', iostat=status) fnam
      if (status /= 0) stop 'unexpected end of list file'
    end do

    ! continue with files that have not yet been scanned
    ldone = .false.
    do
      read(funit, '(a1024)', iostat=status) fnam
      if (status /= 0) then
        fnam = ' '
        irec = 0
        exit
      end if

      ! pick out files that match fpre
      if (index(fnam, trim(fpre)) <= 0) then
        nskip = nskip + 1
        irecold = 0
        cycle
      end if

      ! skip files beyond year1 - yearn (with one-year halo)
      if (.not. scanallfiles) then
        idx = scan(fnam, ".", back=.true.)
        str1 = fnam(1:idx - 1)
        idx = scan(str1, ".", back=.true.)
        str2 = str1(idx + 1:idx + 4)
        read(str2, "(i4)") yyyy
        if (yyyy < y1 - 1 .or. yyyy > y2 + 1) then
          nskip = nskip + 1
          irecold = 0
          cycle
        end if
      end if

      !write(*,*) 'open_ofile,fnam:',trim(fnam)
      ! read time information from files
      call handle_ncerror(nf90_open(trim(fnam), nf90_nowrite, ncid))

      call handle_ncerror(nf90_inquire(ncid, unlimiteddimid=unlimdimid))
      call handle_ncerror(nf90_inquire_dimension(ncid, unlimdimid, len=nrec))

      call handle_ncerror(nf90_inq_varid(ncid, 'time', rhid))
      units = ' '
      call handle_ncerror(nf90_get_att(ncid, rhid, 'units', units))
      call ncsevl(units, n, sstartend)
      read(units(sstartend(1, 3):sstartend(1, 3) + 3), *) yref
      read(units(sstartend(1, 4):sstartend(1, 4) + 1), *) mref
      read(units(sstartend(1, 5):sstartend(1, 5) + 1), *) dref
      calendar = ' '
      call handle_ncerror(nf90_get_att(ncid, rhid, 'calendar', calendar))

      ! translate time limits to date numbers
      call nccaln(trim(calendar), y1, m1, 1, yref, mref, dref, ndays)
      dnumlo = ndays
      if (m2 < 12) then
        call nccaln(trim(calendar), y2, m2 + 1, 1, yref, mref, dref, ndays)
      else
        call nccaln(trim(calendar), y2 + 1, 1, 1, yref, mref, dref, ndays)
      end if
      dnumhi = ndays
      call nccaln(trim(calendar), y1, m1, 1, yr, 1, 1, ndays)
      toff = dnumlo - ndays

      ltimebnds = .true.
      status = nf90_inq_varid(ncid, 'time_bnds', rhid2)
      if (status /= nf90_noerr) then
        status = nf90_inq_varid(ncid, 'time_bounds', rhid2)
        if (status /= nf90_noerr) ltimebnds = .false.
      end if

      do n = irecold + 1, nrec
        irecold = n
        call handle_ncerror(nf90_get_var(ncid, rhid, tval, start=(/n/)))
        if (ltimebnds) then
          call handle_ncerror(nf90_get_var(ncid, rhid2, tbnd, &
            start=(/1, n/), count=(/2, 1/)))
        else
          tbnd = tval
        end if

        ! check if time value lies within limits
        call nccaln(trim(calendar), yref, mref, dref, y1, m1, 1, ndays)
        dnumlo = ndays
        if (m2 < 12) then
          call nccaln(trim(calendar), yref, mref, dref, y2, m2 + 1, 1, ndays)
        else
          call nccaln(trim(calendar), yref, mref, dref, y2 + 1, 1, 1, ndays)
        end if
        dnumhi = ndays

        if (tbnd(1) >= dnumlo .and. tbnd(2) <= dnumhi) then
          irec = n
          tval = tval + toff
          tbnd = tbnd + toff
          tbndold = tbnd + toff
          ldone = .true.
          exit
        elseif (tbnd(1) >= dnumhi) then
          irec = 0
          tval = -999.
          tbnd = -999.
          ldone = .true.
          exit
        end if
      end do

      if (irecold == nrec) then
        irecold = 0
        nskip = nskip + 1
      end if

      call handle_ncerror(nf90_close(ncid))

      if (ldone) exit
    end do

    close(funit)

    ! get bounds of month
    call tbndmon(calendar, yr, 0.5 * (tbnd(1) + tbnd(2)), mbnd, year, month)

  end subroutine get_file_info

  ! -----------------------------------------------------------------

  subroutine tbndmon(calendar, yref, tval, tbnd, year, month)

    implicit none

    character(len=*), intent(in) :: calendar
    integer, intent(in) :: yref
    real(r8), intent(in) :: tval
    real(r8), intent(out) :: tbnd(2)
    integer, intent(out) :: year, month

    integer :: y1, m1, y2, m2, ndays, ndayslast

    if (tval < 0.) return

    ndayslast = 0
    ndays = 0
    y2 = yref
    m2 = 1
    year = yref
    month = 1


    do while (ndays < tval)
      ndayslast = ndays

      year = y2
      month = m2
      if (m2 == 12) then
        y2 = y2 + 1
        m2 = 1
      else
        m2 = m2 + 1
      end if
      call nccaln(calendar, yref, 1, 1, y2, m2, 1, ndays)
    end do
    tbnd(1) = ndayslast
    tbnd(2) = ndays

  end subroutine tbndmon

  ! -----------------------------------------------------------------

  subroutine ncsevl(strg, strgn, strgind)

    ! Description:
    !   Finds the number and the locations of sub-strings.
    !   Valid deliminators are ' ', '-' and ':'.
    !
    ! Arguments:
    !   char(*) strg   (in)  - input string
    !   int strgn      (out) - number of sub-strings
    !   int strgind(*) (out) - start/end locations of sub-strings (the
    !                         dimension must at least equal to strgn*2)

    implicit none

    character(len=*), intent(in) :: strg
    integer, intent(out) :: strgn, strgind(*)

    character :: charold, charnew
    integer :: i

    charold = ' '
    strgn = 0
    do i = 1, len(strg)
      charnew = strg(i:i)
      if ((charold == ' ' .or. charold == '-' .or. charold == ':') .and. &
          (charnew /= ' ' .and. charnew /= '-' .and. charnew /= ':')) then
        strgn = strgn + 1
        strgind(strgn) = i
      elseif ((charnew == ' ' .or. charnew == '-' .or. charnew == ':') &
          .and. (charold /= ' ' .and. charold /= '-' .and. charold /= ':')) then
        strgn = strgn + 1
        strgind(strgn) = i - 1
      end if
      charold = charnew
    end do
    if (mod(strgn, 2) == 1) then
      strgn = strgn + 1
      strgind(strgn) = len(strg)
    end if
    strgn = strgn / 2

  end subroutine ncsevl

  ! -----------------------------------------------------------------

  subroutine ncdnum(calendar, year, month, day, ndays)

    ! Description:
    !   Auxilary subroutine for nccaln which calculates the number
    !   of days between a given date and a reference date. Thereby the
    !   reference date varies from calendar to calendar and is not
    !   standerized. The purpose is more of differential nature rather
    !   than to create an absolute time axis.

    implicit none

    character(len=*), intent(in) :: calendar
    integer, intent(in) :: year, month, day
    integer, intent(out) :: ndays

    integer :: y
    integer, parameter :: yoffset = 1000000, yr1 = 365, yr4 = 1461, &
      yr100 = 36524, yr400 = 146097
    integer :: acc_leap(12), acc_noleap(12)
    data acc_leap   /0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335/
    data acc_noleap /0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/

    y = year + yoffset
    if (calendar(1:3) == 'pro') then
      ndays = int((y - 1) / 400) * yr400 &
        + int(mod((y - 1), 400) / 100) * yr100 &
        + int(mod(mod((y - 1), 400), 100) / 4) * yr4 &
        + mod(mod(mod((y - 1), 400), 100), 4) * yr1
      if (month > 2 .and. ((mod(y, 4) == 0 .and. mod(y, 100) /= 0) .or. mod(y, 400) == 0)) then
        ndays = ndays + acc_leap(month) + day
      else
        ndays = ndays + acc_noleap(month) + day
      end if
    elseif (calendar(1:3) == 'jul') then
      ndays = int((y - 1) / 4) * yr4 + mod(y, 4) * yr1
      if (month > 2 .and. mod(y, 4) == 0) then
        ndays = ndays + acc_leap(month) + day
      else
        ndays = ndays + acc_noleap(month) + day
      end if
    elseif (calendar(1:3) == '360') then
      ndays = (y - 1) * 360 + 30 * (month - 1) + day
    elseif (calendar(1:2) == 'no' .or. calendar(1:3) == '365') then
      ndays = (y - 1) * 365 + acc_noleap(month) + day
    elseif (calendar(1:3) == 'all' .or. calendar(1:3) == '366') then
      ndays = (y - 1) * 366 + acc_leap(month) + day
    end if

  end subroutine ncdnum

  ! -----------------------------------------------------------------

  subroutine nccaln(calendar, y1, m1, d1, y2, m2, d2, ndays)

    ! Description:
    !   nccaln calculates the number of days between two
    !   dates. Valid calendars are '360_day', 'noleap' = '365_day',
    !   'all_leap' = '366_day', 'julian', 'proleptic_gregorian' and
    !   'standard' = 'mixed' = 'gregorian'.
    !
    ! Arguments:
    !   int nccaln       (out) -  number of days
    !   char(*) calendar (in)  -  choice of calendar
    !   int y1           (in)  -  start year
    !   int m1           (in)  -  start month
    !   int d1           (in)  -  start day
    !   int y2           (in)  -  last year
    !   int m2           (in)  -  last month
    !   int d2           (in)  -  last day

    implicit none

    character(len=*), intent(in) :: calendar
    integer, intent(in) :: y1, m1, d1, y2, m2, d2
    integer, intent(out) :: ndays

    integer :: n1, n2, nd1, nd2, nd3

    if (calendar(1:3) == 'sta' .or. calendar(1:3) == 'mix' .or. &
        calendar(1:3) == 'gre') then
      if (y1 < 1582 .or. (y1 == 1582 .and. m1 < 10) .or. &
          (y1 == 1582 .and. m1 == 10 .and. d1 <= 15)) then
        call ncdnum('julian', y1, m1, d1, n1)
      else
        call ncdnum('julian', 1582, 10, 15, nd1)
        call ncdnum('proleptic', 1582, 10, 15, nd2)
        call ncdnum('proleptic', y1, m1, d1, nd3)
        n1 = nd1 - nd2 + nd3 - 10
      end if
      if (y2 < 1582 .or. (y2 == 1582 .and. m2 < 10) .or. &
          (y2 == 1582 .and. m2 == 10 .and. d2 <= 15)) then
        call ncdnum('julian', y2, m2, d2, n2)
      else
        call ncdnum('julian', 1582, 10, 15, nd1)
        call ncdnum('proleptic', 1582, 10, 15, nd2)
        call ncdnum('proleptic', y2, m2, d2, nd3)
        n2 = nd1 - nd2 + nd3 - 10
      end if
    else
      call ncdnum(calendar, y1, m1, d1, n1)
      call ncdnum(calendar, y2, m2, d2, n2)
    end if
    ndays = n2 - n1

  end subroutine nccaln

  ! -----------------------------------------------------------------

  subroutine read_secindex(fname, fexist, seclen, iind, jind, iflg, jflg)

    implicit none

    logical, intent(out) :: fexist
    character(len=*), intent(in) :: fname
    integer, intent(out) :: iind(*), jind(*), iflg(*), jflg(*)
    integer, intent(out) :: seclen

    character(len=100) :: c100
    integer :: iostatus

    inquire(file=trim(fname), exist=fexist)
    if (.not. fexist) return

    seclen = 0
    open(funit, file=trim(fname))
    do
      c100 = ' '
      read(funit, '(a100)', iostat=iostatus) c100
      if (iostatus < 0) exit
      if (index(c100, 'fram_strait') > 0) then
        do
          c100 = ' '
          read(funit, '(a100)', iostat=iostatus) c100
          if (iostatus < 0) exit
          if (index(c100, 'Name') > 0) exit
          seclen = seclen + 1
          read(c100, *, iostat=iostatus) iind(seclen), jind(seclen), &
            iflg(seclen), jflg(seclen)
          if (iostatus < 0) exit
        end do
      end if
    end do
    close(funit)

  end subroutine read_secindex

  ! -----------------------------------------------------------------

  real(r8) function transifs(seclen, iind, jind, iflg, jflg, fldx, fldy)

    implicit none

    integer, intent(in) :: iind(*), jind(*), iflg(*), jflg(*), seclen
    real(r8), dimension(:, :), intent(in) :: fldx, fldy

    integer :: n

    transifs = 0
    do n = 1, seclen
      transifs = transifs - &
        fldx(iind(n) - 1, jind(n)) * iflg(n) - &
        fldy(iind(n), jind(n) - 1) * jflg(n)
    end do

  end function transifs

  ! -----------------------------------------------------------------

  subroutine strmf_eval(idm, jdm, kdm, umflx, vmflx, strmf)

    ! Description: computes stream function from h. velocity components
    ! Author:      Mats Bentsen
    ! Date:        feb2005

    implicit none

    integer, intent(in) :: idm, jdm, kdm
    real(r8), dimension(idm, jdm, kdm), intent(inout) :: umflx, vmflx
    real(r8), dimension(idm, jdm, kdm), intent(out) :: strmf

    integer :: i, j, ip1, jp1

    ! ------------------------------------------------------------------
    ! integrate the stream function with boundary condition strmf(1,1)=0
    ! ------------------------------------------------------------------

    where (umflx == 1e20) umflx = 0
    where (vmflx == 1e20) vmflx = 0
    umflx(:, :, 1) = sum(umflx, 3)
    vmflx(:, :, 1) = sum(vmflx, 3)
    strmf(1, 1, 1) = 0.

    do j = 2, jdm
      strmf(1, j, 1) = strmf(1, j - 1, 1) - umflx(1, j - 1, 1)
    end do

    do j = jdm, 1, -1
      do i = 2, idm
        strmf(i, j, 1) = strmf(i - 1, j, 1) + vmflx(i - 1, j, 1)
      end do
    end do

    ! Move reference to Greenland
    do j = 1, jdm
      do i = 1, idm
        strmf(i, j, 1) = strmf(i, j, 1) - strmf(idm, jdm, 1)
      end do
    end do

    ! ------------------------------------------------------------------
    ! interpolate the streamfunction to the p-point (also smooths)
    ! ------------------------------------------------------------------

    do j = 1, jdm
      jp1 = mod(j, jdm) + 1
      do i = 1, idm - 1
        ip1 = mod(i, idm) + 1
        strmf(i, j, 1) = 0.25 * (strmf(i, j, 1) + strmf(ip1, j, 1) &
          + strmf(i, jp1, 1) + strmf(ip1, jp1, 1))
      end do
    end do

  end subroutine strmf_eval

  ! -----------------------------------------------------------------

  subroutine rotate_uv(idm, jdm, angle, u, v)

    implicit none

    integer, intent(in) :: idm, jdm
    real(r8), dimension(idm, jdm), intent(in) :: angle
    real(r8), dimension(idm, jdm), intent(inout) :: u, v

    integer :: i, j
    real(r8) :: urot

    do j = 1, jdm
      do i = 1, idm
        urot = u(i, j) * cos(angle(i, j)) - v(i, j) * sin(angle(i, j))
        v(i, j) = u(i, j) * sin(angle(i, j)) + v(i, j) * cos(angle(i, j))
        u(i, j) = urot
      end do
    end do

  end subroutine rotate_uv

  ! -----------------------------------------------------------------

  real(r8) function rho(p, th, s)

    ! Description: computes in-situ density from potential temperature
    !              and salinity
    ! Comment: units are in cgs

    implicit none

    real(r8), intent(in) :: p
    real(r8), intent(in) :: th, s

    real(r8), parameter :: &
      a11 = 9.9985372432159340e-01_r8, a12 = 1.0380621928183473e-02_r8, &
      a13 = 1.7073577195684715e-03_r8, a14 = -3.6570490496333680e-05_r8, &
      a15 = -7.3677944503527477e-06_r8, a16 = -3.5529175999643348e-06_r8, &
      b11 = 1.7083494994335439e-10_r8, b12 = 7.1567921402953455e-13_r8, &
      b13 = 1.2821026080049485e-13_r8, a21 = 1.0, &
      a22 = 1.0316374535350838e-02_r8, a23 = 8.9521792365142522e-04_r8, &
      a24 = -2.8438341552142710e-05_r8, a25 = -1.1887778959461776e-05_r8, &
      a26 = -4.0163964812921489e-06_r8, b21 = 1.1995545126831476e-10_r8, &
      b22 = 5.5234008384648383e-13_r8, b23 = 8.4310335919950873e-14_r8

    rho = (a11 + (a12 + a14 * th + a15 * s) * th + (a13 + a16 * s) * s &
      + (b11 + b12 * th + b13 * s) * p) &
      / (a21 + (a22 + a24 * th + a25 * s) * th + (a23 + a26 * s) * s &
      + (b21 + b22 * th + b23 * s) * p)

  end function rho

  ! -----------------------------------------------------------------

  real(r8) function p_alpha(p1, p2, th, s)

    ! Description: integrate specific volume with respect to pressure
    ! Comment: units are in cgs

    implicit none

    real(r8), intent(in) :: p1, p2, th, s

    real(r8), parameter :: &
      a11 = 9.9985372432159340e-01_r8, a12 = 1.0380621928183473e-02_r8, &
      a13 = 1.7073577195684715e-03_r8, a14 = -3.6570490496333680e-05_r8, &
      a15 = -7.3677944503527477e-06_r8, a16 = -3.5529175999643348e-06_r8, &
      b11 = 1.7083494994335439e-10_r8, b12 = 7.1567921402953455e-13_r8, &
      b13 = 1.2821026080049485e-13_r8, a21 = 1.0, &
      a22 = 1.0316374535350838e-02_r8, a23 = 8.9521792365142522e-04_r8, &
      a24 = -2.8438341552142710e-05_r8, a25 = -1.1887778959461776e-05_r8, &
      a26 = -4.0163964812921489e-06_r8, b21 = 1.1995545126831476e-10_r8, &
      b22 = 5.5234008384648383e-13_r8, b23 = 8.4310335919950873e-14_r8

    real(r8), parameter :: r1_3 = 1. / 3., r1_5 = 1. / 5., r1_7 = 1. / 7., r1_9 = 1. / 9.

    real(r8) :: a1, a2, b1, b2, pm, r, q, qq

    a1 = a11 + (a12 + a14 * th + a15 * s) * th + (a13 + a16 * s) * s
    a2 = a21 + (a22 + a24 * th + a25 * s) * th + (a23 + a26 * s) * s
    b1 = b11 + b12 * th + b13 * s
    b2 = b21 + b22 * th + b23 * s

    ! the analytic solution of the integral is
    !   p_alpha=(b2*(p2-p1)
    !           +(a2-a1*b2/b1)*log((a1+b1*p2)/(a1+b1*p1)))/b1
    ! a truncated series expansion of the integral is used that provide
    ! better computational efficiency and accuarcy for most relevant
    ! parameters

    pm = 0.5 * (p2 + p1)
    r = 0.5 * (p2 - p1) / (a1 + b1 * pm)
    q = b1 * r
    qq = q * q

    p_alpha = 2. * r * (a2 + b2 * pm &
      + (a2 - a1 * b2 / b1) * qq * (r1_3 + qq * (r1_5 + qq * (r1_7 + qq * r1_9))))

  end function p_alpha

  ! -----------------------------------------------------------------

  real(r8) function getlpi(temp, saln, phiu, phil, pu)

    ! get lower pressure interface of a layer knowing the temperature,
    ! salinity of the layer and the geopotential at upper and lower
    ! interface

    implicit none

    real(r8), intent(in) :: temp, saln, phiu, phil, pu

    real(r8) :: pl, q, dphi, alpu, alpl

    ! first guess on pressure interface
    pl = pu - rho(pu, temp, saln) * (phil - phiu)

    ! improve the accuracy of the pressure interface by an
    ! iterative procedure
    q = 1.
    do while (abs(q) > 1.e-4)
      call delphi(pu, pl, temp, saln, dphi, alpu, alpl)
      q = (phil - phiu - dphi) / alpl
      pl = pl - q
    end do

    getlpi = pl

  end function getlpi

  ! -----------------------------------------------------------------

  subroutine delphi(p1, p2, th, s, dphi, alp1, alp2)

    ! integrate specific volume with respect to pressure to find the
    ! difference in geopotential between two pressure levels

    implicit none

    real(r8), intent(in) :: p1, p2, th, s
    real(r8), intent(out) :: dphi, alp1, alp2

    real(r8), parameter :: &
      a11 = 9.9985372432159340e-01_r8, a12 = 1.0380621928183473e-02_r8, &
      a13 = 1.7073577195684715e-03_r8, a14 = -3.6570490496333680e-05_r8, &
      a15 = -7.3677944503527477e-06_r8, a16 = -3.5529175999643348e-06_r8, &
      b11 = 1.7083494994335439e-10_r8, b12 = 7.1567921402953455e-13_r8, &
      b13 = 1.2821026080049485e-13_r8, a21 = 1.0, &
      a22 = 1.0316374535350838e-02_r8, a23 = 8.9521792365142522e-04_r8, &
      a24 = -2.8438341552142710e-05_r8, a25 = -1.1887778959461776e-05_r8, &
      a26 = -4.0163964812921489e-06_r8, b21 = 1.1995545126831476e-10_r8, &
      b22 = 5.5234008384648383e-13_r8, b23 = 8.4310335919950873e-14_r8

    real(r8), parameter :: r1_3 = 1. / 3., r1_5 = 1. / 5., r1_7 = 1. / 7., r1_9 = 1. / 9.

    real(r8) :: a1, a2, b1, b2, pm, r, q, qq

    a1 = a11 + (a12 + a14 * th + a15 * s) * th + (a13 + a16 * s) * s
    a2 = a21 + (a22 + a24 * th + a25 * s) * th + (a23 + a26 * s) * s
    b1 = b11 + b12 * th + b13 * s
    b2 = b21 + b22 * th + b23 * s

    ! the analytic solution of the integral is
    !   dphi=-(b2*(p2-p1)
    !         +(a2-a1*b2/b1)*log((a1+b1*p2)/(a1+b1*p1)))/b1
    ! a truncated series expansion of the integral is used that provide
    ! better computational efficiency and accuarcy for most relevant
    ! parameters

    pm = 0.5 * (p2 + p1)
    r = 0.5 * (p2 - p1) / (a1 + b1 * pm)
    q = b1 * r
    qq = q * q

    dphi = -2. * r * (a2 + b2 * pm &
      + (a2 - a1 * b2 / b1) * qq * (r1_3 + qq * (r1_5 + qq * (r1_7 + qq * r1_9))))

    alp1 = (a2 + b2 * p1) / (a1 + b1 * p1)
    alp2 = (a2 + b2 * p2) / (a1 + b1 * p2)

  end subroutine delphi

  ! -----------------------------------------------------------------

  subroutine sphmidpnt(lambda_a, theta_a, lambda_b, theta_b, lambda_c, theta_c)

    ! This subroutine computes a point (theta_c, lambda_c) on the unit
    ! sphere so that it lies on the geodesic curve defined by the points
    ! (theta_a, lambda_a) and (theta_b, lambda_b) and midway between
    ! (theta_b, lambda_b) as (theta_a, lambda_a). (Mats Bentsen 2010)

    implicit none

    real(r8), intent(in) :: theta_a, lambda_a, theta_b, lambda_b
    real(r8), intent(out) :: theta_c, lambda_c

    real(r8) :: x_a, y_a, z_a, x_b, y_b, z_b, beta, x_c, y_c, z_c

    real(r8), parameter :: deg2rad = 3.141592654_r8 / 180._r8, rad2deg = 1._r8 / deg2rad

    ! Represent the spherical coordinates as Cartesian coordinates on a
    ! unit sphere.
    x_a = cos(lambda_a * deg2rad) * cos(theta_a * deg2rad)
    y_a = cos(lambda_a * deg2rad) * sin(theta_a * deg2rad)
    z_a = sin(lambda_a * deg2rad)

    x_b = cos(lambda_b * deg2rad) * cos(theta_b * deg2rad)
    y_b = cos(lambda_b * deg2rad) * sin(theta_b * deg2rad)
    z_b = sin(lambda_b * deg2rad)

    x_c = x_a + x_b
    y_c = y_a + y_b
    z_c = z_a + z_b

    ! Convert from Cartesian coordinates to spherical coordinates
    theta_c = atan2(y_c, x_c) * rad2deg
    lambda_c = atan2(z_c, sqrt(x_c * x_c + y_c * y_c)) * rad2deg

  end subroutine sphmidpnt

  ! -----------------------------------------------------------------

  subroutine sphextpnt(lambda_a, theta_a, lambda_b, theta_b, lambda_c, theta_c)

    ! Description: This subroutine computes a point (theta_c, lambda_c)
    ! on the unit sphere so that it lies on the geodesic curve defined
    ! by the points (theta_a, lambda_a) and (theta_b, lambda_b), at
    ! the same distance from (theta_b, lambda_b) as (theta_a, lambda_a),
    ! but on the opposite side. (Mats Bentsen 2010)

    implicit none

    real(r8), intent(in) :: theta_a, lambda_a, theta_b, lambda_b
    real(r8), intent(out) :: theta_c, lambda_c

    real(r8) :: x_a, y_a, z_a, x_b, y_b, z_b, beta, x_c, y_c, z_c

    real(r8), parameter :: deg2rad = 3.141592654_r8 / 180._r8, rad2deg = 1._r8 / deg2rad

    ! Represent the spherical coordinates as Cartesian coordinates on
    ! a unit sphere.
    x_a = cos(lambda_a * deg2rad) * cos(theta_a * deg2rad)
    y_a = cos(lambda_a * deg2rad) * sin(theta_a * deg2rad)
    z_a = sin(lambda_a * deg2rad)

    x_b = cos(lambda_b * deg2rad) * cos(theta_b * deg2rad)
    y_b = cos(lambda_b * deg2rad) * sin(theta_b * deg2rad)
    z_b = sin(lambda_b * deg2rad)

    beta = 2 * (x_a * x_b + y_a * y_b + z_a * z_b)

    x_c = beta * x_b - x_a
    y_c = beta * y_b - y_a
    z_c = beta * z_b - z_a

    ! Convert from Cartesian coordinates to spherical coordinates
    theta_c = atan2(y_c, x_c) * rad2deg
    lambda_c = atan2(z_c, sqrt(x_c * x_c + y_c * y_c)) * rad2deg

  end subroutine sphextpnt

  ! -----------------------------------------------------------------

  subroutine vinth2pecmwf(dati, dato, hbcofa, hbcofb, p0, plevi, plevo, &
      intyp, psfc, spvl, kxtrp, imax, nlat, nlevi, nlevip1, nlevo, varflg, tbot, phis)

    !     THIS ROUTINE INTERPLOATES CCM2/3 HYBRID COORDINATE DATA
    !     TO PRESSURE COORDINATES USING PRESSURE SURFACES AS THE
    !     COORDINATE SURFACE WHERE THE INTERPOLATION IS DONE.  THE
    !     TYPE OF INTERPOLATION IS CURRENTLY A VARIANT OF TRANSFORMED
    !     PRESSURE COORDINATES WITH THE  INTERPOLATION TYPE
    !     SPECIFIED BY INTYP.  ALL HYBRID COORDINATE VALUES ARE
    !     TRANSFORMED TO PRESSURE VALUES. WHERE THE
    !     FORMULA FOR THE PRESSURE OF A HYBRID SURFACE IS;
    !          P(K) = HBCOFA(LEVH,K)*P0 + HBCOFB(LEVH,K)*PSFC
    !     WHERE,
    !          HBCOFA - IS THE "A" OR PRESSURE HYBRID COEF
    !          LEVH   - IS THE LAYER SURFACE (INTERFACE=1 MIDPOINT=2)
    !          P0     - IS THE BASE PRESSURE IN MB
    !          K      - THE LEVEL INDEX (RUNNING FROM TOP TO BOTTOM)
    !          HBCOFB - IS THE "B" OR SIGMA COEFICIENT
    !          P(K)   - IS THE PRESSURE OF A HYBRID SURFACE IN MB.
    !          PSFC   - IS THE SURFACE PRESSURE IN PASCALS
    !                   (MB = .01*PASCALS)
    !
    !     FOR HYBRID DATA AT LEVEL INTERFACES SINCE THERE IS ONE
    !     MORE VERTICAL LEVEL FOR INTERFACES THAN FOR LEVEL MIDPOINTS
    !     IT IS ASSUNMED THAT THE FIRST INTERFACE LEVEL WITH A DATA
    !     VALUE IS THE SECOND LEVEL FROM THE TOP.
    !
    !     ON INPUT-
    !        DATI    - 3 DIMENSIONAL ARRAY (I,J,KI) CONTAINING DATA
    !                  ON HYBRID SURFACES  WHERE I IS LONGTIUDE, J
    !                  IS LATITUDE AND K IS THE VERTICAL HYBRID
    !                  COORDINATE.  THE VERTICAL DATA RUN TOP TO BOTTOM.
    !                  SIGMA DATA WITH THE DATA ORDERED TOP TO BOTTOM.
    !        HBCOFA  - 2 DIMENSIONAL ARRAY CONTAINING "A" OR PRESSURE
    !                  COEFICIENTS FOR COMPUTING PRESSURE AT A LEVEL.
    !                  ARRAY IS 2XNLEVIP1.  THE 1ST INDEX TAKES ON
    !                  THE VALUE OF EITHER
    !                   1 - FOR LEVEL INTERFACES (OR 1/2 LEVELS) OR;
    !                   2 - FOR LEVEL MIDPOINTS  (OR FULL LEVELS WHERE
    !                       VIRTUALLY ALL VARIABLES ARE LOCATED)
    !                  NOTE THAT COEFICIENTS ARE SCALED TO YIELD A
    !                  PRESSURE IN MB.  THEY ARE ORDERED FROM TOP
    !                  OF THE MODEL TO THE BOTTOM.
    !        HBCOFB  - SAME AS HCOFA BUT FOR THE "B" OR SIGMA COEFICIENT
    !        P0      - BASE PRESSURE IN MB FOR COMPUTING PRESSURE
    !                  OF A HYBRID COORDINATE LEVEL
    !        PLEVI -  1 DIMENSIONAL ARRAY TO HOLD PRESSURE VALUES
    !                  OF HYBRID SURFACES FOR A VERTICAL COLUMN
    !                  SLICE
    !        PLEVO   - LIST OF OUTPUT PRESSURE SURFACES IN MB
    !                  LOW TO HIGH PRESSURE
    !        INTYP   - A FLAG INDICATING INTERPOLATION FOR EACH
    !                  FIELD (1 - LINEAR,2 - LOG ,3 - LOG LOG)
    !                  WHERE EACH INTERPOLATION IS DONE IN TRANSFORMED
    !                  PRESSURE COORDINATES.
    !        PSFC    - MODEL SFC PRESSURE IN PASCALS (WILL BE CONVERTED
    !                  TO MB)
    !        VCOLI   - ARRAY TO STORE A LONGITUDINAL VERTICAL SLICE OF
    !                  INPUT DATA (IMAX BY NLEVI).
    !        VCOLO   - SAME BUT FOR OUTPUT DATA (IMAX BY NLEVO)
    !        IMAX    - LONGITUDINAL DIMENSION OF THE DATA.
    !        NLAT    - LATITUDINAL DIMENSION OF THE DATA.
    !        NLEVI   - NO. OF LEVELS FOR THE HYBRID DATA
    !        NLEVIP1 - NLEVI + 1
    !        NLEVO   - NUMBER OF OUTPUT LEVELS FOR PRESSURE DATA
    !        KXTRP   - FLAG WHICH INDICATES WHETHER OR NOT
    !                  EXTRAPOLATION WILL BE USED WHEN THE OUTPUT
    !                  PRESSURE SURFACE IS BELOW THE LOWEST LEVEL
    !                  OF THE MODEL.
    !                     0 - DON'T EXTRAPOLATE USE SPECIAL VALUE SPVL
    !                     1 - EXTRAPOLATE DATA using ECMWF formulation
    !                         below PSFC
    !        SPVL    - SPECIAL VALUE TO USE WHEN DATA IS NOT
    !                  EXTRAPOLATED
    !        varflg  - flag which indicates the name of the variable
    !                  -1 means geopotential (Z)
    !                  +1 means geopotential (T)
    !                   0 any other variable
    !        tbot    - temperature at level closest to ground
    !        phis    - surface geopotential
    !
    !     ON OUTPUT-
    !        DATO  - 3 DIMENSIONAL ARRAY TO HOLD DATA INTERPOLATED
    !                TO PRESSURE SURFACES.

    implicit none

    integer, intent(in) :: intyp, imax, nlat, nlevi, nlevip1, nlevo, kxtrp, varflg
    real(r8), intent(in) :: dati(imax, nlat, nlevi), hbcofa(nlevip1), hbcofb(nlevip1)
    real(r8), intent(in) :: p0, plevo(nlevo), psfc(imax, nlat)
    real(r8), intent(in) :: spvl, tbot(imax, nlat), phis(imax, nlat)
    real(r8), intent(out) :: dato(imax, nlat, nlevo)
    real(r8), intent(inout) :: plevi(nlevip1)

    integer :: i, j, k, kp, kpi
    real(r8) :: a1, a2ln, a2ln1, a2ln2
    real(r8) :: tstar, hgt, alnp, t0, tplat, tprime0, alpha, alnp3, psfcmb
    real(r8) :: alp, alphp
    real(r8), parameter :: rd = 287.04d0
    real(r8), parameter :: ginv = 1.d0 / 9.80616d0
    real(r8), parameter :: alpha0 = 0.0065d0 * rd * ginv

    ! statement function for double log interpolation
    a2ln(a1) = log(log(a1 + 2.72d0))

    do j = 1, nlat
      do i = 1, imax
        if (psfc(i, j) == spvl) then
          do k = 1, nlevo
            dato(i, j, k) = spvl
          end do
          cycle
        end if

        ! Get pressure values for hybrid surfaces for this point
        ! and for the type of model surface the data is on.
        ! Interface data starts at the second interface level so
        ! if the data is on those levels start the
        do k = 1, nlevi
          kpi = k
          plevi(k) = (hbcofa(kpi) * p0) + hbcofb(kpi) * (psfc(i, j) * 0.01d0)
        end do

        ! Perform vertical interpolation
        do k = 1, nlevo
          ! Check for bracketing level KP will be the input level that
          ! is the upper portion of 2 input bracketing levels.

          ! Branch for model top
          if (plevo(k) <= plevi(1)) then
            kp = 1
            go to 30

            ! Branch for level below lowest hybrid level
          else if (plevo(k) > plevi(nlevi)) then
            if (kxtrp == 0) then
              dato(i, j, k) = spvl
              go to 40
            else if (varflg > 0) then
              ! Variable is "T" and ECMWF extrapolation is desired
              psfcmb = psfc(i, j) * 0.01d0
              tstar = dati(i, j, nlevi) * (1.d0 + alpha0 * (psfcmb / plevi(nlevi) - 1.d0))
              hgt = phis(i, j) * ginv
              if (hgt < 2000.d0) then
                alnp = alpha0 * log(plevo(k) / psfcmb)
              else
                t0 = tstar + 0.0065d0 * hgt
                tplat = min(t0, 298.d0)
                if (hgt <= 2500.d0) then
                  tprime0 = 0.002d0 * ((2500.d0 - hgt) * t0 + (hgt - 2000.d0) * tplat)
                else
                  tprime0 = tplat
                end if
                if (tprime0 < tstar) then
                  alnp = 0.d0
                else
                  alnp = rd * (tprime0 - tstar) / phis(i, j) * log(plevo(k) / psfcmb)
                end if
              end if
              alnp3 = alnp * alnp * alnp
              dato(i, j, k) = tstar * (1.d0 + alnp + 0.5d0 * alnp**2 + 1.d0 / 6.d0 * alnp3)
              go to 40

            else if (varflg < 0) then
              ! Variable is "Z" and ECMWF extrapolation is desired
              psfcmb = psfc(i, j) * 0.01d0
              hgt = phis(i, j) * ginv
              tstar = tbot(i, j) * (1.d0 + alpha0 * (psfcmb / plevi(nlevi) - 1.d0))
              t0 = tstar + 0.0065d0 * hgt

              if (tstar <= 290.5d0 .and. t0 > 290.5d0) then
                alp = rd / phis(i, j) * (290.5d0 - tstar)
              else if (tstar > 290.5d0 .and. t0 > 290.5d0) then
                alp = 0
                tstar = 0.5d0 * (290.5d0 + tstar)
              else
                alp = alpha0
              end if

              if (tstar < 255.d0) then
                tstar = 0.5d0 * (tstar + 255.d0)
              end if
              alnp = alp * log(plevo(k) / psfcmb)
              alnp3 = alnp * alnp * alnp
              dato(i, j, k) = hgt - rd * tstar * ginv * log(plevo(k) / psfcmb) &
                * (1.d0 + 0.5d0 * alnp + 1.d0 / 6.d0 * alnp**2)
              go to 40
            else
              ! Use lowest sigma layer
              dato(i, j, k) = dati(i, j, nlevi)
              go to 40
            end if

            ! Branch to check if output level in between
            ! 2 lowest hybrid levels
          else if (plevo(k) >= plevi(nlevi - 1)) then
            kp = nlevi - 1
            go to 30

            ! Branch for model interior
            ! Loop through input levels till you are bracketing
            ! output level
          else
            kp = 0
20          continue
            kp = kp + 1
            if (plevo(k) <= plevi(kp + 1)) go to 30
            if (kp > nlevi) then
              write(6, fmt=25) kp, nlevi
25            format(' KP.GT.NLEVI IN P2HBD.  KP,NLEVI= ', 2i5)
            end if
            go to 20
          end if

30        continue

          ! Level bracketed, pick type of interpolation
          select case (intyp)
          case (1) ! Linear interpolation
            dato(i, j, k) = dati(i, j, kp) &
              + (dati(i, j, kp + 1) - dati(i, j, kp)) &
              * (plevo(k) - plevi(kp)) / (plevi(kp + 1) - plevi(kp))
          case (2) ! Log interpolation
            dato(i, j, k) = dati(i, j, kp) &
              + (dati(i, j, kp + 1) - dati(i, j, kp)) &
              * log(plevo(k) / plevi(kp)) / log(plevi(kp + 1) / plevi(kp))
          case (3) ! Log-log interpolation
            dato(i, j, k) = dati(i, j, kp) &
              + (dati(i, j, kp + 1) - dati(i, j, kp)) &
              * (a2ln(plevo(k)) - a2ln(plevi(kp))) &
              / (a2ln(plevi(kp + 1)) - a2ln(plevi(kp)))
          end select

40        continue
        end do
      end do
    end do

  end subroutine vinth2pecmwf

  ! -----------------------------------------------------------------

  logical function skip_variable(n, nmax)

    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    integer, intent(in) :: n, nmax

    integer :: mpirank, mpisize, mpierror

#ifdef MPI
    call mpi_comm_size(mpi_comm_world, mpisize, mpierror)
    call mpi_comm_rank(mpi_comm_world, mpirank, mpierror)
#else
    mpirank = 0
    mpisize = 1
#endif

    if (mod(n - 1, mpisize) == mpirank) then
      skip_variable = .false.
    else
      skip_variable = .true.
    end if

  end function skip_variable

end module m_utilities
