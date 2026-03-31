module m_modelsocn

  use netcdf
  use cmor_users_functions
  use m_utilities
  use m_namelists
  use m_jsons

  implicit none

  ! Netcdf variables
  integer, save         :: ii, jj, kk, ncid, rhid, dimid, status

  ! Grid dimensions and variables
  real(r8), save                                :: voglb, aoglb
  integer, save                                     :: idm, jdm, kdm = 0, ddm = 0, ldm = 0, rdm = 0, secdm = 0, slenmax2
  integer, parameter                                :: ncrns = 4
  integer, allocatable, save, dimension(:, :)       :: basin
  real(r8), allocatable, save, dimension(:)     :: xvec, yvec, kvec, kvechalf,&
    sigma, sigmahalf, depth, slat
  real(r8), allocatable, save, dimension(:, :)  :: parea, pmask, pdepth, plon,&
    plat, ulon, ulat, vlon, vlat, slat_bnds, sigma_bnds, sigmahalf_bnds, &
    depth_bnds, bpini, bpinit, uscaley, vscalex, udepth, vdepth
  real(r8), allocatable, save, dimension(:, :, :)   :: plon_crns, plat_crns, &
    ulon_crns, ulat_crns, vlon_crns, vlat_crns, plon_crnsp, plat_crnsp, &
    ulon_crnsp, ulat_crnsp, vlon_crnsp, vlat_crnsp, dzini, sini, tini
  character(len=slenmax), allocatable, save, dimension(:)   :: region1, section1
  character, allocatable, save, dimension(:, :)             :: region, section
  character(len=slenmax), save                              :: tcoord, zcoord, s1
  character(len=slenmax), save                              :: grid, grid_label

  ! Gravity
  real(r8), parameter :: g = 9.80665, ginv = 1.0 / g

  ! Dataset related variables
  character(len=slenmax), save          :: ivnm, ovnm, vunits, vpositive, vtype
  character(len=slenmax), save          :: bvnm, cvnm
  character(len=slenmax * 10), save     :: vcomment
  !character(len=slenmax)               :: key, value
  logical, save :: lsumz
  logical       :: found

  ! Table related variables
  character(len=slenmax), save          :: table, tablepath
  character(len=slenmax), save          :: tabledir_mapping, table_mapping

  ! Cmor parameters
  character(len=1024)   :: fnmo
  integer, save         :: iaxid, jaxid, kaxid, laxid, raxid, saxid, taxid, &
    grdid, varid, table_id, table_id_grid, error_flag

  ! String for module special
  character(len=slenmax), save          :: special

  ! Data fields
  real(r8), allocatable, save, dimension(:, :, :)   :: fld, fld2, fldtmp, &
    fldacc, fldhalf, dp
  real(r8), allocatable, save, dimension(:, :)      :: sealv, pbot
  real(r8)                                          :: sfac, offs, fill

  ! Auxillary variables for special operations
  character(len=slenmax), save                          :: str1, str2

    character(len=slenmax), dimension(:), allocatable  :: vars,preproc_keys
    integer, dimension(:), allocatable                 :: idx
    real(r8), dimension(:), allocatable            :: facs

contains

  ! -----------------------------------------------------------------

  subroutine ocn2cmor

    implicit none

    logical :: badrec, last, first
    !logical :: badrec
    !integer :: k, m, n, nrec
    integer :: k, m, n
    integer :: romon = 365*10*2
    !character(len=slenmax) :: realm, frequency


    badrec = .false.

    ! Print start information
    if (verbose) then
      write(*, *)
      write(*, *) '----------------------------'
      write(*, *) '--- Process ocean output ---'
      write(*, *) '----------------------------'
      write(*, *)
    end if

!   ! Process table Omon
!   write(*, *) 'Process table Omon'
    pomon = ''
    fnm = pomon

    tabledir_mapping='/diagnostics/CMOR/esm2cmor/recipes/template/'
    table_mapping='variable_mapping_NorESM3_to_CMIP7.json'

    ! filter only ocean variables, facilitate parallisation
    n = count(realms == 'ocean' .or. realms == 'ocnBgchem')
    allocate(idx(n))
    k=1
    do n = 1, n_variables
      if (trim(realms(n)) == 'ocean' .or. trim(realms(n)) == 'ocnBgchem' ) then
        idx(k) = n
        k = k + 1
      end if
    end do
    write(*,*) 'idx:',idx

    write(*,*) 'realms:'
    !do n = 1, size(idx)
      !write(*,*) trim(realms(idx(n)))
    !end do

    !do n = 1, n_variables
    main_loop: do n = 1, size(idx)
      realm = trim(realms(idx(n)))
      !if (realm /= 'ocean') cycle
!     !if (skip_variable(n, nomon, domon)) cycle
      if (skip_variable(n, n_variables)) cycle

!     ! Map namelist variables
      bvnm = trim(branded_names(idx(n)))
      cvnm = trim(compound_names(idx(n)))
      frequency = trim(frequencies(idx(n)))
      ovnm = bvnm
      table = 'CMIP7_'//trim(realm)//'.json'

      select case (frequency)
      case('mon')
        if (realm == 'ocean') then
          itag = tagomon
        else
          itag = tagomonbgc
        end if
      case('day')
        if (realm == 'ocean') then
          itag = tagoday
        else
          itag = tagodaybgc
        endif
      case('yr')
        if (realm == 'ocean') then
          itag = tagoyr
        else
          itag = tagoyrbgc
        end if
      end select
      write(*,*) 'realm:',trim(realm)
      write(*,*) 'cvnm:',trim(cvnm)
      write(*,*) 'frequency:',trim(frequency)
      write(*,*) 'itag:',trim(itag)


      !! STORE file list for each tag
    ! Read grid information from input files
    !write(*, *) 'Read grid information from input files'
    call scan_files(reset=.true.)

    if (len_trim(fnm) == 0) then
        if (verbose) write(*, *) &
          'WARNING: no file found for case dir|tag|year1|month1|yearn|monthn: ', &
          trim(ibasedir) // '/' // trim(casename), '|', trim(itag), '|', &
          year1, '|', month1, '|', yearn, '|', monthn
      !else

        !isloop = .false.
        !end if
        cycle
      end if
    !if (len_trim(fnm) == 0) return
    call read_gridinfo_ifile
    !stop 'l168'


      !ivnm = 
!     special = vomon(3, n)
      !vunits = ' '
      !call json_get_vunits(trim(tabledir)//trim(table),trim(ovnm),vunits)
      !call json_get_vunits('/diagnostics/CMOR/esm2cmor/tables/CMIP7_ocean.json','tos_tavg-u-hxy-sea',vunits)
      !call json_get_value(trim(tabledir)//trim(table), &
           !'variable_entry.' // trim(ovnm) // '.units',vunits,found)
      call json_get_units(trim(tabledir)//trim(table), trim(ovnm),vunits)
      !write(*,*) 'get_vunits:',trim(vunits)
!     vpositive = ' '
!     vcomment = ' '

!     ! Check if vertical coordinate required
      !write(*, *) 'l381,tabledir/table:', trim(tabledir)//trim(table)
      !write(*, *) 'ovnm:', trim(ovnm)
      call json_get_vertcoord(trim(tabledir)//trim(table), bvnm, zcoord,lfound=found)
      !write(*, *) 'l382,zcoord:', trim(zcoord)

!     ! Choose history file
!     !if (index(special, 'day2mon') > 0) then
        !itag = tagoday
!     !else
!       itag = tagomon
!     !end if

!     ! Check if input variable is present
      !if (len_trim(pomon) == 0) call scan_files(reset=.true.)
      ! call scan_files(reset=.true.)
      !if (.not. var_in_file(fnm, ivnm)) cycle

      write(*,*) 'cvnm:',trim(cvnm)
      !call json_get_array_string(trim(tabledir_mapping)//trim(table_mapping),&
        !'variable_entry:'//trim(cvnm)//':sources:vars',&
       !vars, separator=':',lfound=found) 
      call json_get_vars(trim(tabledir_mapping)//trim(table_mapping),&
                trim(cvnm), vars, lfound=found) 
      if (found) then
        !call json_get_array_real(trim(tabledir_mapping)//trim(table_mapping),&
          !'variable_entry:'//trim(cvnm)//':sources:facs',&
         !facs,separator=':', lfound=found) 
      call json_get_facs(trim(tabledir_mapping)//trim(table_mapping),&
                trim(cvnm), facs, lfound=found) 
        !write(*,*) 'size(vars):',size(vars)
        !write(*,*) 'vars:',vars
        do k =1, size(vars)
          !write(*,*) 'vars(k):',trim(vars(k))
          if (.not. var_in_file(fnm, vars(k))) cycle main_loop
        end do
        ivnm = vars(1)
      else
        !call json_get_val_str(trim(tabledir_mapping)//trim(table_mapping),&
          !'variable_entry:'//trim(cvnm)//':original_name',&
          !ivnm,separator=':',lfound=found) 
        call json_get_original_name(trim(tabledir_mapping)//trim(table_mapping), &
                trim(cvnm), ivnm, lfound=found) 
        if (.not. found) cycle
        !write(*,*) 'original_name:',trim(ivnm)
        !ivnm = value_json
        if (.not. var_in_file(fnm, ivnm)) cycle
      end if

      !!write(*,*) 'facs:',facs
      !if (found) then
        !!write(*,*) 'size(facs):',size(facs)
        !do k =1, size(facs)
          !!write(*,*) 'facs(k):',trim(facs(k))
          !if (.not. var_in_file(fnm, vars(k))) cycle main_loop
        !end do
        !!ivnm = vars(1)
      !else
        !call json_get_val_str(trim(tabledir_mapping)//trim(table_mapping),&
          !'variable_entry:'//trim(cvnm)//':original_name',&
          !ivnm,separator=':',lfound=found) 
        !if (.not. found) cycle
        !!write(*,*) 'original_name:',trim(ivnm)
        !!ivnm = value_json
        !if (.not. var_in_file(fnm, ivnm)) cycle
      !end if

!     ! Prepare output file
      call special_pre

!     ! Loop over input files
      m = 0
      do
        m = m + 1

!       ! Open output file
        !write(*,*) 'm:',m
        if (mod(m - 1, romon) == 0) then
            call open_ofile(ivnm,ovnm)
        end if

!       ! Read variable into buffer (average if necessary)
        rec = 0
        !nrec = 0
        !fldacc = 0.
        !last = .false.
        !do
          !if (len_trim(pomon) == 0) call scan_files(reset=.false.)
          call scan_files(reset=.false.)
          !write(*,*) 'fnm:',trim(fnm)
          !write(*,*) 'rec:',rec
          if (rec == 0) exit
          !write(*,*) 'rec:',rec
          !if (rec == 0) then
            !last = .true.
            !exit
          !end if
          !nrec = nrec + 1
          !write(*,*) 'nrec:',nrec
          !write(*,*) 'fnm:',fnm
          call read_tslice(rec, badrec, fnm)
          !fldacc = fldacc + fld
          !if (index(special, 'day2mon') > 0) then
            !if (tbnd(2) + 0.5 >= mbnd(2)) exit
          !else
            !exit
          !end if
          !write(*,*) 'last:',last
          !if (last) exit
        !end do
        !fld = fldacc / real(nrec)

        !! for monthly
      select case (frequency)
      case('mon')
        tbnds(:, 1) = mbnd
        tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))
      case('day')
        tbnds(1, 1) = tval(1) - 0.5
        tbnds(2, 1) = tval(1) + 0.5
      case('yr')
         tbnds(1, 1) = tval(1) - 365. / 2.
         tbnds(2, 1) = tval(1) + 365. / 2.
      end select

        !write(*,*) 'tval:',tval
        !write(*,*) 'tbnds:',tbnds

        ! Post processing
        call special_post

!       ! Write time slice to output file
        call write_tslice

!       ! Close output file if max rec has been reached
        if (mod(m, romon) == 0) call close_ofile

      end do

!     ! Close output file if still open
      if (mod(m, romon) > 0) call close_ofile


    if (allocated(sigma))          deallocate(sigma)
    if (allocated(sigmahalf))      deallocate(sigmahalf)
    if (allocated(sigma_bnds))     deallocate(sigma_bnds)
    if (allocated(sigmahalf_bnds)) deallocate(sigmahalf_bnds)

    if (allocated(depth))          deallocate(depth)
    if (allocated(depth_bnds))      deallocate(depth_bnds)

    if (allocated(slat))          deallocate(slat)
    if (allocated(slat_bnds))      deallocate(slat_bnds)

    if (allocated(section))          deallocate(section)
    if (allocated(section1))          deallocate(section1)

    if (allocated(region))          deallocate(region)
    if (allocated(region1))          deallocate(region1)

    deallocate(parea, pmask, pdepth, &
      plon, plat, bpini, bpinit, &
      ulon, ulat, vlon, vlat, &
      plon_crns, plat_crns, &
      ulon_crns, ulat_crns, &
      vlon_crns, vlat_crns, &
      plon_crnsp, plat_crnsp, &
      ulon_crnsp, ulat_crnsp, &
      vlon_crnsp, vlat_crnsp, &
      sealv, xvec, yvec, kvec, pbot, &
      dzini, sini, tini, &
      kvechalf, uscaley, vscalex, &
      udepth, vdepth, basin, stat=status)

    if (allocated(vars)) deallocate(vars)
    if (allocated(facs)) deallocate(facs)
    if (allocated (preproc_keys)) deallocate(preproc_keys)

    end do main_loop

    if (allocated(idx)) deallocate(idx)

!   ! Process table fx
!   write(*, *) 'Process table fx'
!   fnm = trim(griddata)//trim(ocngridfile)
!   table = tfx
!   do n = 1, nfx
!     if (skip_variable(n, nfx, dfx)) cycle

!     ! Map namelist variables
!     ovnm = vfx(ovnmpos, n)
!     ivnm = vfx(ivnmpos, n)
!     special = vfx(3, n)
!     vunits = ' '
!     vpositive = ' '
!     vcomment = ' '

!     ! Check if input variable is present
!     if (.not. var_in_file(fnm, ivnm)) cycle

!     ! Check if vertical coordinate required
!     call json_get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)

!     ! Prepare output file
!     call special_pre
!     call open_ofile(fx=.true.)

!     ! Read field
!     call read_field

!     ! Post Processing
!     call special_post

!     ! Write field
!     call write_field

!     ! Close output file
!     call close_ofile

!   end do

!   ! Process table Ofx
!   write(*, *) 'Process table Ofx'
!   fnm = trim(griddata)//trim(ocngridfile)
!   table = tofx
!   do n = 1, nofx
!     if (skip_variable(n, nofx, dofx)) cycle

!     ! Map namelist variables
!     ovnm = vofx(ovnmpos, n)
!     ivnm = vofx(ivnmpos, n)
!     special = vofx(3, n)
!     vunits = ' '
!     vpositive = ' '
!     vcomment = ' '

!     ! Use oceanregnfile for region
!     if (ovnm == 'basin') then
!       fnm = trim(griddata)//trim(ocnregnfile)
!     else
!       fnm = trim(griddata)//trim(ocngridfile)
!     end if

!     ! Check if input variable is present
!     if (.not. var_in_file(fnm, ivnm)) cycle

!     ! Check if vertical coordinate required
!     call json_get_vertcoord(trim(tabledir)//trim(table), ovnm, zcoord)

!     ! Prepare output file
!     call special_pre
!     call open_ofile(fx=.true.)

!     ! Read field
!     call read_field

!     ! Post Processing
!     call special_post

!     ! Write field
!     call write_field

!     ! Close output file
!     call close_ofile

!   end do

  end subroutine ocn2cmor

  ! -----------------------------------------------------------------

  subroutine special_pre

    implicit none

    integer :: i, j, k, n

    character(len=slenmax), dimension(:), allocatable  :: preproc_keys
    character(len=slenmax)        :: preproc_key, preproc_val

    lsumz = .false.
    !str2 = special
      !tabledir='/diagnostics/CMOR/esm2cmor/recipes/template/'
    !table='variable_mapping_NorESM3_to_CMIP7.json'
    !call json_get_keys(trim(tabledir_mapping)//trim(table_mapping),&
        !'variable_entry:'//trim(cvnm)//':preproc',&
        !preproc_keys,separator=':',lfound=found)
    call json_get_preproc_keys(trim(tabledir_mapping)//trim(table_mapping),&
        trim(cvnm), preproc_keys, lfound=found)

    !write(*,*) 'preproc:'
    !write(*,*) 'len_trim(preproc):',len_trim(preproc)
    !write(*,*) 'size(preproc):',size(preproc)

    if (.not. found) return
    write(*,*) 'preproc_keys:', preproc_keys

    do n=1,size(preproc_keys)
      write(*,*) 'n:',n
      preproc_key = preproc_keys(n)
      write(*,*) 'preproc_key:', trim(preproc_key)
      !call json_get_val_str(trim(tabledir_mapping)//trim(table_mapping),&
          !'variable_entry:'//trim(cvnm)//':preproc:'//trim(preproc_key),&
          !preproc_val,lfound=found)
      call json_get_preproc_val(trim(tabledir_mapping)//trim(table_mapping),&
          trim(cvnm),trim(preproc_key), preproc_val, lfound=found)
      if (found) then
        write(*,*) trim(preproc_key),":",trim(preproc_val)
      else
        cycle
      end if

    !do
      !if (index(str2, ';') > 0) then
        !str1 = str2(1:index(str2, ';') - 1)
        !str2 = str2(index(str2, ';') + 1:)
      !else
        !str1 = str2
      !end if
      !select case (str1)
      select case (preproc_val)

        ! atm to Pa
      case ('atm2Pa')
        vunits = 'Pa'

        ! uatm to Pa
      case ('uatm2Pa')
        vunits = 'Pa'

        ! Unit transformation: mol cfcXX m-3 -> mol cfcXX kg-1
      case ('cfcunits')
        vunits = 'mol kg-1'

        ! CFC11 comment
      case ('cfc11comment')
        vcomment = 'In this simulation, annual means of reconstructed' &
          // ' Northern Hemisphere CFC-11 are applied globally to the ocean.' &
          // ' Reference: Walker S.J., Weiss R.F., Salameh P.K. (2000)' &
          // ' Reconstructed histories of the annual mean atmospheric mole' &
          // ' fractions for the halocarbons CFC-11, CFC-12, CFC-113 and' &
          // ' carbon tetrachloride. J. Geophys. Res. 105(C6): 14285-14296.' &
          // ' CFC-11 data in ppt (1910.5-2008.5):' &
          // ' 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,' &
          // ' 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,' &
          // ' 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,' &
          // ' 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.2, 0.4, 0.7, 1.0,' &
          // ' 1.5, 2.2, 3.0, 4.1, 5.4, 6.8, 8.1, 9.4, 11.1, 13.3, 16.1,' &
          // ' 19.6, 23.8, 28.4, 33.6, 39.5, 46.1, 53.7, 62.5, 72.0, 82.7,' &
          // ' 94.9, 108.4, 121.4, 133.9, 145.9, 156.6, 168.3, 176.7, 184.3,' &
          // ' 191.4, 199.4, 208.1, 218.1, 229.5, 241.7, 253.0, 259.5,' &
          // ' 266.0, 268.4, 268.3, 269.7, 269.8, 268.5, 267.3, 265.9,' &
          // ' 264.4, 262.9, 262.0, 260.3, 258.1, 256.0, 254.1, 252.0,' &
          // ' 248.9, 246.9, 245.3'

        ! CFC12 comment
      case ('cfc12comment')
        vcomment = 'In this simulation, annual means of reconstructed' &
          // ' Northern Hemisphere CFC-12 are applied globally to the ocean.' &
          // ' Reference: Walker S.J., Weiss R.F., Salameh P.K. (2000)' &
          // ' Reconstructed histories of the annual mean atmospheric mole' &
          // ' fractions for the halocarbons CFC-11, CFC-12, CFC-113 and' &
          // ' carbon tetrachloride. J. Geophys. Res. 105(C6): 14285-14296.' &
          // ' CFC-12 data in ppt (1910.5-2008.5):' &
          // ' 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,' &
          // ' 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,' &
          // ' 0.0, 0.0, 0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1.2, 1.7,' &
          // ' 2.3, 3.4, 4.8, 6.1, 7.6, 9.2, 11.0, 12.8, 15.0, 17.4, 20.2,' &
          // ' 23.4, 26.8, 30.5, 35.0, 40.0, 45.8, 52.5, 60.4, 69.3, 79.2,' &
          // ' 90.3, 102.8, 116.7, 132.0, 148.3, 166.1, 185.7, 207.1, 228.1,' &
          // ' 248.0, 266.9, 284.2, 305.9, 323.0, 339.6, 353.3, 369.0,' &
          // ' 385.7, 403.4, 424.3, 444.0, 465.4, 483.6, 497.7, 506.0,' &
          // ' 516.4, 523.2, 528.9, 533.9, 537.6, 540.5, 542.5, 544.1,' &
          // ' 546.4, 546.9, 546.8, 546.5, 547.2, 546.5, 544.5, 541.8, 540.6'

        ! mol m-3
      case ('mol m-3')
        vunits = 'mol m-3'

        ! mol m-2 s-1
      case ('mol m-2 s-1')
        vunits = 'mol m-2 s-1'

        ! mol m-3 s-1
      case ('mol m-3 s-1')
        vunits = 'mol m-3 s-1'

        ! Salinity has to be in practical salinity units
      case ('psu')
        vunits = 'psu'

        ! Fix unitless units
      case ('unitless')
        vunits = '1'

        ! Set correct units for percentage
      case ('percent')
        vunits = '%'

        ! Set correct units for fraction
      case ('fraction')
        vunits = '1'

        ! Unit transformation: kg m-2
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

        ! Fix wo units
      case ('wflx2wo')
        vunits = 'm s-1'

        ! Set units to kg
      case ('kg')
        vunits = 'kg'

        ! Set units to m^3
      case ('m3')
        vunits = 'm3'

        ! Convert units from radians2 to m2
      case ('rad2m')
        vunits = 'm2'

        ! Set units for streamfunction
      case ('strmf')
        vunits = 'kg s-1'
        lsumz = .true.

        ! Set positive attribute
      case ('positiveup')
        vpositive = 'up'
      case ('positivedo')
        vpositive = 'down'

        ! Compute vertical sum
      case ('sumz')
        lsumz = .true.

        ! Compute density
      case ('ts2rho0')
        vunits = 'kg m-3'

        ! Compute steric sea level change from density
      case ('ts2zossga', 't2zostoga')
        vunits = 'm'

        ! Compute fixed cell volume of interpolated grid
      case ('volcello')
        vunits = 'm3'

        ! Compute fixed cell thickness
      case ('thkcello')
        vunits = 'm'

        ! Compute depth of local mimina
      case ('locminz')
        vunits = 'm'

        ! Compute depth of local mimina
      case ('omega2z')
        vunits = 'm'

        ! Convert units from [mol P/m3] to [kg Chl/m3] using [60 gC/gChl]
      case ('kg m-3')
        vunits = 'kg m-3'

        ! Convert units from [m3] to [1e3 km3]
      case ('1e3 km3')
        vunits = '1e3 km3'

        ! Convert units from [s-1] to [s-2], fix bug for the units of bfsq in micom
      case ('s-2')
        vunits = 's-2'

        ! Set unit g m-2
      case ('g m-2')
        vunits = 'g m-2'

        ! Set unit degC kg m-2
      case ('degC kg m-2')
        vunits = 'degC kg m-2'

      end select
      !if (str1 == str2) exit
    end do

    if (allocated(preproc_keys)) deallocate(preproc_keys)

  end subroutine special_pre

  ! -----------------------------------------------------------------

  subroutine special_post

    implicit none

    integer         :: i, j, k, n
    real(r8)    :: r, rd, p, ptoptmp, pbottmp, sref = 35.0

    !character(len=slenmax), dimension(:), allocatable  :: preproc
    !character(len=:), allocatable  :: preproc_key, preproc_val
    character(len=slenmax), dimension(:), allocatable  :: postproc_keys
    !character(len=:), allocatable  :: preproc_key, preproc_val
    character(len=slenmax) :: postproc_key, postproc_val

    !call json_get_keys(trim(tabledir_mapping)//trim(table_mapping),&
        !'variable_entry:'//cvnm//':postproc',&
        !postproc,separator=':',lfound=found)
    call json_get_postproc_keys(trim(tabledir_mapping)//trim(table_mapping),&
        trim(cvnm), postproc_keys, lfound=found)
    !write(*,*) 'postproc:'
    !write(*,*) 'len_trim(postproc):',len_trim(postproc)
    !write(*,*) 'size(postproc):',size(postproc)

    if (.not. found) return

    do n=1,size(postproc_keys)
      postproc_key = postproc_keys(n)
      !call json_get_val_str(trim(tabledir_mapping)//trim(table_mapping),&
          !'variable_entry:'//cvnm//':postproc:'//postproc_key,&
          !postproc_val,separator=':',lfound=found)
      call json_get_postproc_val(trim(tabledir_mapping)//trim(table_mapping),&
          trim(cvnm),trim(postproc_key), postproc_val, lfound=found)
      if (found) then
      !if (len_trim(postproc_val) >0 ) then
        write(*,*) trim(postproc_key),":",trim(postproc_val)
      else
        cycle
      end if

    !str2 = special
    !do
      !if (index(str2, ';') > 0) then
        !str1 = str2(1:index(str2, ';') - 1)
        !str2 = str2(index(str2, ';') + 1:)
      !else
        !str1 = str2
      !end if
      select case (postproc_val)

        ! Compute depth below geoid from dz or pddpo
      case ('dz2zfull')
        fldacc(:, :, 1) = fld(:, :, 1) * 0.5 - sealv
        do k = 2, kdm
          fldacc(:, :, k) = fldacc(:, :, k - 1) + (fld(:, :, k - 1) + fld(:, :, k)) * 0.5
        end do
        fld = fldacc

        ! Compute depth below geoid at interfaces from dz
      case ('dz2zhalf')
        fldacc(:, :, 1) = -sealv
        do k = 2, kdm
          fldacc(:, :, k) = fldacc(:, :, k - 1) + fld(:, :, k - 1)
        end do
        fld = fldacc

        ! Convert units from radians2 to m2
      case ('rad2m')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) < 1e20) &
                fld(i, j, k) = fld(i, j, k) * 6.37122e6**2
            end do
          end do
        end do

        ! Set ice free points to missing value
      case ('zero2missing')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (abs(fld(i, j, k)) < 1e-6) fld(i, j, k) = 1e20
            end do
          end do
        end do

        ! Set ice free points to missing value
      case ('pmask')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (pmask(i, j) < 0.5) fld(i, j, k) = 1e20
            end do
          end do
        end do

        ! Compute vertical sum
      case ('sumz')
        do k = 2, kk
          do j = 1, jj
            do i = 1, ii
              if (abs(fld(i, j, k)) < 1e20) &
                fld(i, j, 1) = fld(i, j, 1) + fld(i, j, k)
            end do
          end do
        end do

        ! Compute local minima
      case ('locmin')
        do j = 1, jj
          do i = 1, ii
            do k = 2, kk
              if (fld(i, j, 1) > fld(i, j, k)) then
                fld(i, j, 1) = fld(i, j, k)
              else
                exit
              end if
            end do
          end do
        end do

        ! Compute local minima depth
      case ('locminz')
        fld2 = fld
        do j = 1, jj
          do i = 1, ii
            fld2(i, j, 1) = depth(1)
            do k = 2, kk
              if (fld(i, j, 1) > fld(i, j, k)) then
                fld(i, j, 1) = fld(i, j, k)
                fld2(i, j, 1) = depth(k)
              else
                exit
              end if
            end do
          end do
        end do
        fld = fld2

        ! Compute Aragonite Saturation Depth
      case ('omega2z')
        fldtmp = 1e20
        do j = 1, jj
          do i = 1, ii
            do k = 1, kk
              if (fld(i, j, k) < 1.0) then
                if (k == 1) then
                  fldtmp(i, j, 1) = 0.0
                else
                  fldtmp(i, j, 1) = depth(k)
                end if
              end if
            end do
          end do
        end do
        fld = fldtmp

        ! Multiply with global ocean area
      case ('xglbarea')
        fld = fld * aoglb

        ! Multiply gravity constant
      case ('xg')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) < 1e20) fld(i, j, k) = fld(i, j, k) * 9.806
            end do
          end do
        end do

        ! Devide by gravity constant
      case ('xginv')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) < 1e20) fld(i, j, k) = fld(i, j, k) / 9.806
            end do
          end do
        end do

        ! Flip sign
      case ('flipsign')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) < 1e20) fld(i, j, k) = -fld(i, j, k)
            end do
          end do
        end do

        ! Compute global 2d average
      case ('glbave2d')
        fld(1, 1, 1) = sum(fld(:, :, 1) * parea) / sum(parea)

        ! Compute global 3d average
      case ('glbave3d')
        r = 0.
        rd = 0.
        do k = 1, kk
          r = r + sum(fld(:, :, k) * dp(:, :, k) * parea)
          rd = rd + sum(dp(:, :, k) * parea)
        end do
        fld(1, 1, 1) = r / max(1e-10, rd)

        ! Compute potential density
      case ('ts2rho0')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) /= 1e20) &
                fld(i, j, k) = 1e3 * rho(0., fld(i, j, k), fld2(i, j, k))
            end do
          end do
        end do

        ! Compute steric sea level
      case ('ts2zossga')
        if (vtype == 'layer') then
          r = 0.
          do j = 1, jj
            do i = 1, ii
              pbottmp = 0.
              do k = 1, kk
                ptoptmp = pbottmp
                pbottmp = ptoptmp + dp(i, j, k)
                if (fld(i, j, k) /= 1e20) then
                  r = r + 1e-2 * (1e-2 * ginv) * parea(i, j) * &
                    p_alpha(ptoptmp, pbottmp, fld(i, j, k), fld2(i, j, k))
                end if
              end do
            end do
          end do
          fld(1, 1, 1) = (r - voglb) / aoglb
          write(*, *) 'zossga=', fld(1, 1, 1)
        else
          stop 'input variables for zossga must be of type layer'
        end if

        ! Compute thermo-steric sea level
      case ('t2zostoga')
        if (vtype == 'layer') then
          r = 0.
          do j = 1, jj
            do i = 1, ii
              pbottmp = 0.
              do k = 1, kk
                ptoptmp = pbottmp
                pbottmp = ptoptmp + dp(i, j, k)
                if (fld(i, j, k) /= 1e20) &
                  r = r + 1e-2 * (1e-2 * ginv) * parea(i, j) * &
                    p_alpha(ptoptmp, pbottmp, fld(i, j, k), sref)
              end do
            end do
          end do
          fld(1, 1, 1) = (r - voglb) / aoglb
          write(*, *) 'zostoga', fld(1, 1, 1)
        else
          stop 'input variable for zostoga must be of type layer'
        end if

        ! Set land mask of streamfunction
      case ('strmf')
        do j = 1, jj
          do i = 1, ii
            if (pmask(i, j) == 0) fld(i, j, 1) = 1e20
          end do
        end do

        ! Compute fixed cell volume of interpolated grid
      case ('volcello')
        do j = 1, jj
          do i = 1, ii
            do k = ddm, 1, -1
              if (fld(i, j, 1) == 0.) then
                fld(i, j, k) = 1e20
              else
                ptoptmp = min(depth_bnds(1, k), fld(i, j, 1))
                pbottmp = min(depth_bnds(2, k), fld(i, j, 1))
                fld(i, j, k) = (pbottmp - ptoptmp) * fld2(i, j, 1)
              end if
            end do
          end do
        end do

        ! Compute vertical velocity form vertical mass flux
      case ('wflx2wo')
        do j = 1, jj
          do i = 1, ii
            do k = 1, kk
              if (fld(i, j, k) /= 1e20) then
                fld(i, j, k) = fld(i, j, k) / (1035. * parea(i, j))
              end if
            end do
          end do
        end do

        ! Compute fixed cell thickness
      case ('thkcello')
        do j = 1, jj
          do i = 1, ii
            do k = ddm, 1, -1
              if (fld(i, j, 1) == 0.) then
                fld(i, j, k) = 1e20
              else
                fld(i, j, k) = min(fld(i, j, 1), depth_bnds(2, k)) - &
                  min(fld(i, j, 1), depth_bnds(1, k))
              end if
            end do
          end do
        end do

        ! Compute basin index
      case ('basin')
        open(10, file=trim(griddata)//trim(ocnmertfile))
        read(10, '(2i6)') i, j
        if (i /= idm .or. j /= jdm) &
          stop 'mertraocean: incorrect indexes in mertraoceans.dat!'
        str1 = ' '
        write(str1, *) '(', jdm, 'i1)'
        read(10, str1) ((basin(i, j), j = 1, jdm), i = 1, idm)
        close(10)
        fld = 0
        do j = 1, jdm
          do i = 1, idm
            ! Southern Ocean
            if (plat(i, j) < 0. .and. basin(i, j) == 1) fld(i, j, 1) = 1
            ! Pacific Ocean
            if (basin(i, j) == 3) fld(i, j, 1) = 3
            ! Arctic Ocean
            if (plat(i, j) > 60. .and. basin(i, j) == 1) fld(i, j, 1) = 4
            ! Indian Ocean
            if (basin(i, j) == 4) fld(i, j, 1) = 5
            ! Mediterranean Sea
            if (basin(i, j) == 2 .and. plat(i, j) > 30.5 .and. &
              plat(i, j) < 40.5 .and. (plon(i, j) > 354.5 .or. &
              plon(i, j) < 37)) fld(i, j, 1) = 6
            if (basin(i, j) == 2 .and. plat(i, j) > 40.5 .and. &
              plat(i, j) < 46. .and. (plon(i, j) > 359. .or. &
              plon(i, j) < 27.5)) fld(i, j, 1) = 6
            ! Black Sea
            if (basin(i, j) == 1 .and. plat(i, j) > 40.5 .and. &
              plat(i, j) < 48. .and. plon(i, j) > 27.5 .and. &
              plon(i, j) < 45) fld(i, j, 1) = 7
            ! Hudson Bay
            if (basin(i, j) == 2 .and. plat(i, j) > 50. .and. &
              plat(i, j) < 70. .and. plon(i, j) > 265. .and. &
              plon(i, j) < 295) fld(i, j, 1) = 8
            ! Baltic Sea
            if (basin(i, j) == 2 .and. plat(i, j) > 53. .and. &
              plat(i, j) < 62. .and. plon(i, j) > 10. .and. &
              plon(i, j) < 30) fld(i, j, 1) = 9
            if (basin(i, j) == 2 .and. plat(i, j) > 62. .and. &
              plat(i, j) < 66.5 .and. plon(i, j) > 17. .and. &
              plon(i, j) < 30) fld(i, j, 1) = 9
            ! Red Sea
            if (basin(i, j) == 4 .and. plat(i, j) > 13. .and. &
              plat(i, j) < 30. .and. plon(i, j) > 31. .and. &
              plon(i, j) < 44) fld(i, j, 1) = 10
            ! Atlantic ocean
            if (basin(i, j) == 2 .and. fld(i, j, 1) == 0) fld(i, j, 1) = 2
          end do
        end do

        ! Multiple a second field
      case ('Xfield2')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) /= 1e20) &
                fld(i, j, k) = fld(i, j, k) * fld2(i, j, k)
            end do
          end do
        end do

        ! Divide a second field
      case ('Dfield2')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) /= 1e20) &
                fld(i, j, k) = fld(i, j, k) / fld2(i, j, k)
            end do
          end do
        end do

        ! Integratal with respect to depth
      case ('dpint')
        do j = 1, jj
          do i = 1, ii
            if (fld(i, j, 1) /= 1e20) &
              fld(i, j, 1) = fld(i, j, 1) * fld2(i, j, 1) / 9.806
            do k = 2, kk
              if (fld(i, j, k) /= 1e20) &
                fld(i, j, k) = fld(i, j, k - 1) + fld(i, j, k) * fld2(i, j, k) / 9.806
            end do
          end do
        end do

        ! Average with respect to pressure
      case ('dpavg')
        do j = 1, jj
          do i = 1, ii
            if (fld(i, j, 1) /= 1e20) &
              fld(i, j, 1) = fld(i, j, 1) * fld2(i, j, 1)
            do k = 2, kk
              if (fld(i, j, k) /= 1e20) then
                fld(i, j, k) = fld(i, j, k - 1) + fld(i, j, k) * fld2(i, j, k)
                fld2(i, j, 1) = fld2(i, j, 1) + fld2(i, j, k)
              end if
            end do
            if (fld(i, j, 1) /= 1e20) &
              fld(i, j, 1) = fld(i, j, 1) / fld2(i, j, 1)
          end do
        end do

        ! Average over upper 300 m
      case ('dzavg300')
        fldtmp = 1e20
        do j = 1, jj
          do i = 1, ii
            if (fld(i, j, 1) /= 1e20) then
              fldtmp(i, j, 1) = (min(300., pdepth(i, j), depth_bnds(2, 1)) &
                - min(300., pdepth(i, j), depth_bnds(1, 1))) &
                * pbot(i, j) / pdepth(i, j)
              fld(i, j, 1) = fld(i, j, 1) * fldtmp(i, j, 1)
            end if
            do k = 2, kk
              if (fld(i, j, k) /= 1e20) then
                fldtmp(i, j, k) = (min(300., pdepth(i, j), depth_bnds(2, k)) &
                  - min(300., pdepth(i, j), depth_bnds(1, k))) &
                  * pbot(i, j) / pdepth(i, j)
                fld(i, j, 1) = fld(i, j, 1) + fld(i, j, k) * fldtmp(i, j, k)
                fldtmp(i, j, 1) = fldtmp(i, j, 1) + fldtmp(i, j, k)
              end if
            end do
            if (fld(i, j, 1) /= 1e20) &
              fld(i, j, 1) = fld(i, j, 1) / fldtmp(i, j, 1)
          end do
        end do

        ! Average over upper 700 m
      case ('dzavg700')
        fldtmp = 1e20
        do j = 1, jj
          do i = 1, ii
            if (fld(i, j, 1) /= 1e20) then
              fldtmp(i, j, 1) = (min(700., pdepth(i, j), depth_bnds(2, 1)) &
                - min(700., pdepth(i, j), depth_bnds(1, 1))) &
                * pbot(i, j) / pdepth(i, j)
              fld(i, j, 1) = fld(i, j, 1) * fldtmp(i, j, 1)
            end if
            do k = 2, kk
              if (fld(i, j, k) /= 1e20) then
                fldtmp(i, j, k) = (min(700., pdepth(i, j), depth_bnds(2, k)) &
                  - min(700., pdepth(i, j), depth_bnds(1, k))) &
                  * pbot(i, j) / pdepth(i, j)
                fld(i, j, 1) = fld(i, j, 1) + fld(i, j, k) * fldtmp(i, j, k)
                fldtmp(i, j, 1) = fldtmp(i, j, 1) + fldtmp(i, j, k)
              end if
            end do
            if (fld(i, j, 1) /= 1e20) &
              fld(i, j, 1) = fld(i, j, 1) / fldtmp(i, j, 1)
          end do
        end do

        ! Average over upper 2000 m
      case ('dzavg2000')
        fldtmp = 1e20
        do j = 1, jj
          do i = 1, ii
            if (fld(i, j, 1) /= 1e20) then
              fldtmp(i, j, 1) = (min(2000., pdepth(i, j), depth_bnds(2, 1)) &
                - min(2000., pdepth(i, j), depth_bnds(1, 1))) &
                * pbot(i, j) / pdepth(i, j)
              fld(i, j, 1) = fld(i, j, 1) * fldtmp(i, j, 1)
            end if
            do k = 2, kk
              if (fld(i, j, k) /= 1e20) then
                fldtmp(i, j, k) = (min(2000., pdepth(i, j), depth_bnds(2, k)) &
                  - min(2000., pdepth(i, j), depth_bnds(1, k))) &
                  * pbot(i, j) / pdepth(i, j)
                fld(i, j, 1) = fld(i, j, 1) + fld(i, j, k) * fldtmp(i, j, k)
                fldtmp(i, j, 1) = fldtmp(i, j, 1) + fldtmp(i, j, k)
              end if
            end do
            if (fld(i, j, 1) /= 1e20) &
              fld(i, j, 1) = fld(i, j, 1) / fldtmp(i, j, 1)
          end do
        end do

        ! Unit transformation: mol cfcXX m-3 -> mol cfcXX kg-1
      case ('cfcunits')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) /= 1e20) fld(i, j, k) = fld(i, j, k) / 1027.
            end do
          end do
        end do

        ! atm to Pa
      case ('atm2Pa')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) /= 1e20) fld(i, j, k) = fld(i, j, k) * 101325
            end do
          end do
        end do

        ! uatm to Pa
      case ('uatm2Pa')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) /= 1e20) fld(i, j, k) = fld(i, j, k) * 0.101325
            end do
          end do
        end do

        ! Convert units from [mol P/m3] to [kg Chl/m3] using [60 gC/gChl]
      case ('mol P m-3 -> kg Chl m-3')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) /= 1e20) &
                fld(i, j, k) = fld(i, j, k) * 122 * 12 / 60 / 1000
            end do
          end do
        end do

        ! Convert units from [m3] to [1e3 km3]
      case ('1e3 km3')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) /= 1e20) &
                fld(i, j, k) = fld(i, j, k) / 1e3 / 1e9
            end do
          end do
        end do

        ! Mask grid points in the southern hemisphere
      case ('masks')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (plat(i, j) < 0.0) fld(i, j, k) = 1e20
            end do
          end do
        end do

        ! Mask grid points in the southern hemisphere
      case ('maskn')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (plat(i, j) > 0.0) fld(i, j, k) = 1e20
            end do
          end do
        end do

        ! Multiple parea
      case ('xparea')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) /= 1e20) fld(i, j, k) = fld(i, j, k) * parea(i, j)
            end do
          end do
        end do

        ! Iron to phosphorous ratio in organic matter
      case ('fe2ph')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) /= 1e20) &
                fld(i, j, k) = fld(i, j, k) * 5. * 122. * 1.e-6
            end do
          end do
        end do

        ! Carbon to iron
      case ('c2fe')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) /= 1e20) fld(i, j, k) = fld(i, j, k) * 5. * 1.e-6
            end do
          end do
        end do

        ! epc100 to epn100
      case ('epc100toepn100')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) /= 1e20) fld(i, j, k) = fld(i, j, k) / 122. * 16.
            end do
          end do
        end do

        ! epc100 to epp100
      case ('epc100toepp100')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) /= 1e20) fld(i, j, k) = fld(i, j, k) / 122.
            end do
          end do
        end do

        ! percent
      case ('percent')
        do k = 1, kk
          do j = 1, jj
            do i = 1, ii
              if (fld(i, j, k) /= 1e20) fld(i, j, k) = fld(i, j, k) * 100.
            end do
          end do
        end do

      end select

      !if (str1 == str2) exit
    end do

    if (allocated(preproc_keys)) deallocate(preproc_keys)

  end subroutine special_post

  ! -----------------------------------------------------------------

  subroutine read_gridinfo_ifile

    implicit none

    logical         :: check
    integer         :: i, j, k, n, fid
    real(r8)    :: missing, phiu, phil

    ! Open first input file
    call scan_files(reset=.true.)
    !write(*,*) 'fnm:',trim(fnm)

    status = nf90_open(fnm, nf90_nowrite, ncid)
    call handle_ncerror(status)

    ! Read dimensions
    status = nf90_inq_dimid(ncid, 'x', dimid)
    call handle_ncerror(status)
    status = nf90_inquire_dimension(ncid, dimid, len=idm)
    call handle_ncerror(status)

    status = nf90_inq_dimid(ncid, 'y', dimid)
    call handle_ncerror(status)
    status = nf90_inquire_dimension(ncid, dimid, len=jdm)
    call handle_ncerror(status)

    status = nf90_inq_dimid(ncid, 'layer', dimid)
    if (status == nf90_noerr) then
      status = nf90_inquire_dimension(ncid, dimid, len=kdm)
      call handle_ncerror(status)
      allocate(sigma(kdm), sigmahalf(kdm + 1), sigma_bnds(2, kdm), &
        sigmahalf_bnds(2, kdm + 1), stat=status)
      if (status /= 0) stop 'cannot ALLOCATE enough memory (1b)'
      status = nf90_inq_varid(ncid, 'sigma', rhid)
      call handle_ncerror(status)
      status = nf90_get_var(ncid, rhid, sigma)
      call handle_ncerror(status)
      sigma_bnds(1, 1) = sigma(1) - 0.5 * (sigma(2) - sigma(1))
      sigma_bnds(2, 1) = 0.5 * (sigma(2) + sigma(1))
      do k = 2, kdm - 1
        sigma_bnds(1, k) = 0.5 * (sigma(k) + sigma(k - 1))
        sigma_bnds(2, k) = 0.5 * (sigma(k) + sigma(k + 1))
      end do
      sigma_bnds(1, kdm) = 0.5 * (sigma(kdm) + sigma(kdm - 1))
      sigma_bnds(2, kdm) = sigma(kdm) + 0.5 * (sigma(kdm) - sigma(kdm - 1))
      sigmahalf(1:kdm) = sigma_bnds(1, 1:kdm)
      sigmahalf(kdm + 1) = sigma_bnds(2, kdm)
      sigmahalf_bnds(1, 2:kdm + 1) = sigma
      sigmahalf_bnds(2, 1:kdm) = sigma
      sigmahalf_bnds(1, 1) = sigmahalf(1)
      sigmahalf_bnds(2, kdm + 1) = sigmahalf(kdm + 1)
    end if

    status = nf90_inq_dimid(ncid, 'depth', dimid)
    if (status == nf90_noerr) then
      status = nf90_inquire_dimension(ncid, dimid, len=ddm)
      call handle_ncerror(status)
      allocate(depth(ddm), depth_bnds(2, ddm), stat=status)
      if (status /= 0) stop 'cannot ALLOCATE enough memory (1c)'
      status = nf90_inq_varid(ncid, 'depth', rhid)
      call handle_ncerror(status)
      status = nf90_get_var(ncid, rhid, depth)
      call handle_ncerror(status)
      status = nf90_inq_varid(ncid, 'depth_bnds', rhid)
      call handle_ncerror(status)
      status = nf90_get_var(ncid, rhid, depth_bnds)
      call handle_ncerror(status)
    end if

    write(*, *) 'read lat'
    status = nf90_inq_dimid(ncid, 'lat', dimid)
    if (status == nf90_noerr) then
      status = nf90_inquire_dimension(ncid, dimid, len=ldm)
      call handle_ncerror(status)
      allocate(slat(ldm), slat_bnds(2, ldm), stat=status)
      if (status /= 0) stop 'cannot ALLOCATE enough memory (1c)'
      status = nf90_inq_varid(ncid, 'lat', rhid)
      call handle_ncerror(status)
      status = nf90_get_var(ncid, rhid, slat)
      call handle_ncerror(status)
      slat_bnds(1, 1) = max(-90., slat(1) - 0.5 * (slat(2) - slat(1)))
      slat_bnds(2, 1) = 0.5 * (slat(2) + slat(1))
      do j = 2, ldm - 1
        slat_bnds(1, j) = 0.5 * (slat(j) + slat(j - 1))
        slat_bnds(2, j) = 0.5 * (slat(j) + slat(j + 1))
      end do
      slat_bnds(1, ldm) = 0.5 * (slat(ldm) + slat(ldm - 1))
      slat_bnds(2, ldm) = min(90., slat(ldm) + 0.5 * (slat(ldm) - slat(ldm - 1)))
    end if

    write(*, *) 'read region'
    status = nf90_inq_varid(ncid, 'region', rhid)
    if (status == nf90_noerr) then
      status = nf90_inq_dimid(ncid, 'region', dimid)
      call handle_ncerror(status)
      status = nf90_inquire_dimension(ncid, dimid, len=rdm)
      call handle_ncerror(status)
      status = nf90_inq_dimid(ncid, 'slenmax', dimid)
      status = nf90_inquire_dimension(ncid, dimid, len=slenmax2)
      call handle_ncerror(status)
      allocate(region(slenmax2, rdm), region1(rdm), stat=status)
      if (status /= 0) stop 'cannot ALLOCATE enough memory (1d)'
      status = nf90_inq_varid(ncid, 'region', rhid)
      call handle_ncerror(status)
      write(*, *) 'read region', slenmax2, rdm, shape(region)
      status = nf90_get_var(ncid, rhid, region, (/1, 1/), (/slenmax2, rdm/))
      call handle_ncerror(status)
      region1 = ' '
      do i = 1, rdm
        do j = 1, slenmax2
          region1(i)(j:j) = region(j, i)
        end do
      end do
    end if

    write(*, *) 'read section'
    status = nf90_inq_varid(ncid, 'section', rhid)
    if (status == nf90_noerr) then
      status = nf90_inq_dimid(ncid, 'section', dimid)
      call handle_ncerror(status)
      status = nf90_inquire_dimension(ncid, dimid, len=secdm)
      call handle_ncerror(status)
      status = nf90_inq_dimid(ncid, 'slenmax', dimid)
      status = nf90_inquire_dimension(ncid, dimid, len=slenmax2)
      call handle_ncerror(status)
      allocate(section(slenmax2, secdm), section1(secdm), stat=status)
      if (status /= 0) stop 'cannot ALLOCATE enough memory (1d)'
      status = nf90_inq_varid(ncid, 'section', rhid)
      call handle_ncerror(status)
      !write(*, *) 'secdm,slenmax2:', secdm, slenmax2
      status = nf90_get_var(ncid, rhid, section, (/1, 1/), (/slenmax2, secdm/))
      call handle_ncerror(status)
      !write(*, *) 'section:', section
      section1 = ' '
      do i = 1, secdm
        !write(*, *) 'i:', i
        s1 = ' '
        do j = 1, slenmax2
          s1(j:j) = section(j, i)
        end do
        !write(*, *) 's1:', s1
        if (trim(s1) == 'taiwan_and_luzon_straits') then
          section1(i) = 'taiwan_luzon_straits'
        else
          section1(i) = trim(s1)
        end if
        !k = k + 1
      end do
    end if

    !write(*, *) 'l1898'
    ! Read calendar information (change reference year)
    status = nf90_inq_varid(ncid, 'time', rhid)
    call handle_ncerror(status)
    status = nf90_get_att(ncid, rhid, 'calendar', calendar)
    call handle_ncerror(status)
    status = nf90_get_att(ncid, rhid, 'units', calunits)
    call handle_ncerror(status)
    read(calunits(12:15), '(i4.4)') exprefyear

    ! Close first file
    status = nf90_close(ncid)
    call handle_ncerror(status)

    ! Read longitudes, latitudes
    allocate(parea(idm, jdm), pmask(idm, jdm), pdepth(idm, jdm), &
      plon(idm, jdm), plat(idm, jdm), bpini(idm, jdm), bpinit(idm, jdm), &
      ulon(idm, jdm), ulat(idm, jdm), vlon(idm, jdm), vlat(idm, jdm), &
      plon_crns(idm, jdm, ncrns), plat_crns(idm, jdm, ncrns), &
      ulon_crns(idm, jdm, ncrns), ulat_crns(idm, jdm, ncrns), &
      vlon_crns(idm, jdm, ncrns), vlat_crns(idm, jdm, ncrns), &
      plon_crnsp(ncrns, idm, jdm), plat_crnsp(ncrns, idm, jdm), &
      ulon_crnsp(ncrns, idm, jdm), ulat_crnsp(ncrns, idm, jdm), &
      vlon_crnsp(ncrns, idm, jdm), vlat_crnsp(ncrns, idm, jdm), &
      sealv(idm, jdm), xvec(idm), yvec(jdm), kvec(kdm), pbot(idm, jdm), &
      dzini(idm, jdm, kdm), sini(idm, jdm, kdm), tini(idm, jdm, kdm), &
      kvechalf(kdm + 1), uscaley(idm, jdm), vscalex(idm, jdm), &
      udepth(idm, jdm), vdepth(idm, jdm), basin(idm, jdm), stat=status)
    if (status /= 0) stop 'cannot ALLOCATE enough memory (1)'

    forall (i = 1:idm) xvec(i) = i
    forall (j = 1:jdm) yvec(j) = j
    forall (k = 1:kdm) kvec(k) = k - 0.5
    forall (k = 1:kdm + 1) kvechalf(k) = k - 1

    ! Open grid file
    status = nf90_open(trim(griddata)//trim(ocngridfile), nf90_nowrite, ncid)
    call handle_ncerror(status)

    ! Read grid cell mask, area and bathymetry
    status = nf90_inq_varid(ncid, 'pdepth', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, pdepth)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'udepth', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, udepth)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'vdepth', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, vdepth)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'pmask', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, pmask)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'parea', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, parea)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'udy', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, uscaley)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'vdx', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, vscalex)
    call handle_ncerror(status)
    parea = parea * pmask

    ! Compute global ocean volume and area
    voglb = sum(parea * pdepth)
    aoglb = sum(parea)

    ! Read coordinates
    write(*, *) 'line1974'
    status = nf90_inq_varid(ncid, 'plon', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, plon)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'plat', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, plat)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'ulon', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, ulon)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'ulat', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, ulat)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'vlon', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, vlon)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'vlat', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, vlat)
    call handle_ncerror(status)

    ! Read grid cell vertices
    status = nf90_inq_varid(ncid, 'pclon', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, plon_crns)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'pclat', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, plat_crns)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'uclon', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, ulon_crns)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'uclat', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, ulat_crns)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'vclon', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, vlon_crns)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'vclat', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, vlat_crns)
    call handle_ncerror(status)
    write(*, *) 'line2025'

    ! Permute to compensate for dimension bug in CMOR
    do j = 1, jdm
      do i = 1, idm
        do n = 1, ncrns
          plon_crnsp(n, i, j) = plon_crns(i, j, n)
          plat_crnsp(n, i, j) = plat_crns(i, j, n)
          ulon_crnsp(n, i, j) = ulon_crns(i, j, n)
          ulat_crnsp(n, i, j) = ulat_crns(i, j, n)
          vlon_crnsp(n, i, j) = vlon_crns(i, j, n)
          vlat_crnsp(n, i, j) = vlat_crns(i, j, n)
          if (plon_crnsp(n, i, j) < 0.) &
            plon_crnsp(n, i, j) = plon_crnsp(n, i, j) + 360
          if (ulon_crnsp(n, i, j) < 0.) &
            ulon_crnsp(n, i, j) = ulon_crnsp(n, i, j) + 360
          if (vlon_crnsp(n, i, j) < 0.) &
            vlon_crnsp(n, i, j) = vlon_crnsp(n, i, j) + 360
        end do
        if (plon(i, j) < 0.) plon(i, j) = plon(i, j) + 360
        if (ulon(i, j) < 0.) ulon(i, j) = ulon(i, j) + 360
        if (vlon(i, j) < 0.) vlon(i, j) = vlon(i, j) + 360
      end do
    end do

    ! Close grid file
    status = nf90_close(ncid)
    call handle_ncerror(status)

    return
    ! Read initial layer profile from inicon.nc
    status = nf90_open(trim(griddata)//trim(ocninitfile), nf90_nowrite, ncid)
    call handle_ncerror(status)

    status = nf90_inq_varid(ncid, 'dz', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, dzini)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'saln', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, sini)
    call handle_ncerror(status)
    status = nf90_inq_varid(ncid, 'temp', rhid)
    call handle_ncerror(status)
    status = nf90_get_var(ncid, rhid, tini)
    call handle_ncerror(status)

    status = nf90_close(ncid)
    call handle_ncerror(status)

    ! Compute initial bottom pressure
    bpini = 0.
    bpinit = 0.
    do j = 1, jdm
      do i = 1, idm
        if (pmask(i, j) > 0.5) then
          phil = 0.
          do k = 1, kdm
            phiu = phil
            phil = phiu - 1e+4 * g * dzini(i, j, k)
            bpini(i, j) = getlpi(tini(i, j, k), sini(i, j, k), phiu, phil, bpini(i, j))
            bpinit(i, j) = getlpi(tini(i, j, k), 35., phiu, phil, bpinit(i, j))
          end do
        end if
      end do
    end do

  end subroutine read_gridinfo_ifile

  ! -----------------------------------------------------------------

  subroutine open_ofile(ivnm,ovnm,fx)

    implicit none

    logical, optional, intent(in)   :: fx
    logical                         :: fxflag

    character(len=*), intent(in)   :: ivnm,ovnm

    !real                            :: fac1, fac2, fac3, fac4, fac5, fac6
    integer, parameter              :: ndimmax = 10
    integer                 :: i, j, k, n, ndims, dimids(ndimmax), dimlens(ndimmax)
    !character(len=slenmax)  :: coord, ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6
    character(len=slenmax)  :: coord

    real(r8), allocatable       :: tmp1d(:), tmp2d(:, :)

    ! Check if output variable should have time coordinate
    fxflag = .false.
    if (present(fx)) then
      if (fx) fxflag = .true.
    end if

    ! Inquire variable units and dimensions in input file
    status = nf90_open(fnm, nf90_nowrite, ncid)
    call handle_ncerror(status)

!   call resolve_vnm(slenmax, ivnm, ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6, &
!     fac1, fac2, fac3, fac4, fac5, fac6)
!   if (verbose) then
!     write(*, *) 'Resolve input variable term ', trim(ivnm), ' ='
!     if (len_trim(ivnm1) > 0) write(*, *) ' ', trim(ivnm1), '*', fac1
!     if (len_trim(ivnm2) > 0) write(*, *) ' + ', trim(ivnm2), '*', fac2
!     if (len_trim(ivnm3) > 0) write(*, *) ' + ', trim(ivnm3), '*', fac3
!     if (len_trim(ivnm4) > 0) write(*, *) ' + ', trim(ivnm4), '*', fac4
!     if (len_trim(ivnm5) > 0) write(*, *) ' + ', trim(ivnm5), '*', fac5
!     if (len_trim(ivnm6) > 0) write(*, *) ' + ', trim(ivnm6), '*', fac6
!   end if
    status = nf90_inq_varid(ncid, trim(ivnm), rhid)
    if (status /= nf90_noerr) then
      write(*, *) 'cannot find input variable ', trim(ivnm)
      stop
    end if
    status = nf90_inquire_variable(ncid, rhid, ndims=ndims)
    call handle_ncerror(status)
    status = nf90_inquire_variable(ncid, rhid, dimids=dimids(1:ndims))
    call handle_ncerror(status)
    dimlens = 1
    do n = 1, ndims
      status = nf90_inquire_dimension(ncid, dimids(n), len=dimlens(n))
      call handle_ncerror(status)
    end do
    if (.not. allocated(dp)) allocate(dp(idm, jdm, kdm))
    if (allocated(fld)) deallocate(fld, fld2, fldacc, fldtmp)
    ii = idm
    jj = jdm
    kk = kdm
    write(*, *) 'dimlens(3):', dimlens(3)
    write(*, *) 'kdm:', kdm
    if (dimlens(3) == kdm .and. kdm > 0) then
      vtype = 'layer'
    else if (dimlens(3) == ddm .and. ddm > 0 .or. &
      index(special, 'volcello') > 0 .or. &
      index(special, 'thkcello') > 0 .or. &
      index(special, 'masscello') > 0) then
      vtype = 'level'
      kk = ddm
    else if (dimlens(1) == idm .and. dimlens(2) == jdm .and. ndims <= 3) then
      vtype = '2d'
      kk = 1
    else if (dimlens(1) == ldm .and. dimlens(2) == kdm) then
      vtype = 'merk'
      ii = ldm
      jj = kdm
      kk = rdm
    else if (dimlens(1) == ldm .and. dimlens(2) == ddm) then
      vtype = 'merd'
      ii = ldm
      jj = ddm
      kk = rdm
    else if (dimlens(1) == ldm .and. dimlens(2) == rdm) then
      vtype = 'mert'
      ii = ldm
      jj = rdm
      kk = 1
    else if (dimlens(1) == secdm .and. ndims == 2) then
      vtype = 'sect'
      ii = secdm
      jj = 1
      kk = 1
    else if (dimlens(1) == 1 .and. ndims == 1) then
      vtype = '1d'
      ii = 1
      jj = 1
      kk = 1
    end if
    allocate(fld(ii, jj, kk), fld2(ii, jj, kk), fldacc(ii, jj, kk), &
      fldtmp(ii, jj, kk), stat=status)
    if (status /= 0) stop 'cannot ALLOCATE enough memory (4)'
    if (index(special, 'half') > 0) then
      if (allocated(fldhalf)) deallocate(fldhalf)
      allocate(fldhalf(idm, jdm, kdm + 1), stat=status)
      if (status /= 0) stop 'cannot ALLOCATE enough memory (5)'
    end if

    if (len_trim(vunits) == 0) then
      status = nf90_get_att(ncid, rhid, 'units', vunits)
      call handle_ncerror(status)
      if (trim(vunits) == 'mm/s') vunits = 'kg m-2 s-1'
    end if

    coord = ' '
    status = nf90_get_att(ncid, rhid, 'coordinates', coord)
    if (status /= nf90_noerr) coord(1:1) = ivnm(1:1)

    status = nf90_close(ncid)
    call handle_ncerror(status)

    ! Derive path of CMOR table
    tablepath = trim(tabledir)//trim(table)

    ! Inquire time dimension of output variable
    !write(*, *) 'tablepath:', trim(tablepath)
    !write(*, *) 'ovm:', trim(ovnm)
    !write(*, *) 'vtype:', trim(vtype)
    if (.not. fxflag) call json_get_timecoord(trim(tablepath), ovnm, tcoord)

    ! Call CMOR setup
    if (verbose) then
      if (createsubdirs) then
        error_flag = cmor_setup(inpath=trim(ibasedir), &
          netcdf_file_action=CMOR_REPLACE_4, set_verbosity=CMOR_NORMAL, &
          exit_control=CMOR_EXIT_ON_WARNING, &
          create_subdirectories=1)
      else
        error_flag = cmor_setup(inpath=trim(ibasedir), &
          netcdf_file_action=CMOR_REPLACE_4, set_verbosity=CMOR_NORMAL, &
          exit_control=CMOR_EXIT_ON_WARNING, &
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

    ! Derive path to CMOR table
    tablepath = trim(tabledir)//trim(table)

    ! Define output dataset
    grid_label = ocngrid_label
    grid = trim(ocngrid)
    if (trim(vtype) == 'layer' .and. .not. &
      (lsumz .or. index(special, 'glbave') > 0 &
      .or. index(special, '2zos') > 0 &
      .or. index(special, 'level1') > 0)) then
      grid = trim(ocngrid)//', vertical density coordinate'
    else
      if (index(special, 'glbave') > 0 &
        .or. index(special, '2zos') > 0) then
        !grid_label = 'gm'
        grid = 'global mean or integral'
      else
        if (vtype(1:3) == 'mer' .or. ovnm(1:7) == 'hfbasin') then
          !grid_label = 'grz'
          grid = 'zonal mean or integral'
        else if (trim(vtype) == 'level') then
          !grid_label = 'gr'
          !grid_label = 'g999'
          grid = trim(ocngrid)//', interpolated to z-levels'
        end if
      end if
    end if
    call write_namelist_json(grid, grid_label, ocngrid_resolution, ovnm)
    error_flag = cmor_dataset_json(namelist_file_json)
    !call system('rm '//trim(namelist_file_json))

    ! Define horizontal axes
    write(*, *) 'Define horizontal axes'
    if (vtype(1:3) /= 'mer' .and. vtype(1:3) /= 'sec') then
      iaxid = cmor_axis( &
        !table=trim(tabledir)//trim(tgrids), &
        table=trim(tabledir)//'CMIP7_grids.json', &
        table_entry='i_index', &
        units='1', &
        length=idm, &
        coord_vals=xvec)
      jaxid = cmor_axis( &
        !table=trim(tabledir)//trim(tgrids), &
        table=trim(tabledir)//'CMIP7_grids.json', &
        table_entry='j_index', &
        units='1', &
        length=jdm, &
        coord_vals=yvec)

      write(*, *) 'Define horizontal grid '//coord(1:1)
      if (coord(1:1) == 'p') then
        grdid = cmor_grid( &
          axis_ids=(/iaxid, jaxid/), &
          latitude=plat, &
          longitude=plon, &
          latitude_vertices=plat_crnsp, &
          longitude_vertices=plon_crnsp)
      else if (coord(1:1) == 'u') then
        grdid = cmor_grid( &
          axis_ids=(/iaxid, jaxid/), &
          latitude=ulat, &
          longitude=ulon, &
          latitude_vertices=ulat_crnsp, &
          longitude_vertices=ulon_crnsp)
      else if (coord(1:1) == 'v') then
        grdid = cmor_grid( &
          axis_ids=(/iaxid, jaxid/), &
          latitude=vlat, &
          longitude=vlon, &
          latitude_vertices=vlat_crnsp, &
          longitude_vertices=vlon_crnsp)
      end if
    end if

    if (vtype(1:3) == 'mer') then
      laxid = cmor_axis( &
        table=trim(tablepath), &
        table_entry='latitude', &
        units='degrees_north', &
        length=ldm, &
        coord_vals=slat, &
        cell_bounds=slat_bnds)
      raxid = cmor_axis( &
        table=trim(tablepath), &
        table_entry='basin', &
        units='1', &
        coord_vals=region1)
    end if

    ! Define vertical axis
    write(*, *) "line 2428"
    write(*, *) 'vtype:', trim(vtype)
    if (trim(vtype) == 'layer' .and. .not. &
      (lsumz .or. index(special, 'glbave') > 0 &
      .or. index(special, '2zos') > 0 &
      .or. index(special, 'level1') > 0)) then
      if (index(special, 'half') > 0) then
        kaxid = cmor_axis( &
          table=trim(tablepath), &
          table_entry='rho', &
          units='kg m-3', &
          length=kdm + 1, &
          coord_vals=1000. + sigmahalf, &
          cell_bounds=1000. + sigmahalf_bnds)
      else
        write(*, *) "line 2451"
        kaxid = cmor_axis( &
          table=trim(tablepath), &
          table_entry='rho', &
          units='kg m-3', &
          length=kdm, &
          coord_vals=1000. + sigma, &
          cell_bounds=1000. + sigma_bnds)
        write(*, *) "line 2459"
      end if
    else if (index(special, 'level1') > 0) then
      kaxid = cmor_axis( &
        table=trim(tablepath), &
        table_entry='depth0m', &
        units='m', &
        length=1, &
        coord_vals=(/0/))
    else if (index(special, 'dzavg300') > 0) then
      kaxid = cmor_axis( &
        table=trim(tablepath), &
        table_entry='depth300m', &
        units='m', &
        length=1, &
        coord_vals=(/150./), &
        cell_bounds=(/0., 300./))
    else if (index(special, 'dzavg700') > 0) then
      kaxid = cmor_axis( &
        table=trim(tablepath), &
        table_entry='depth700m', &
        units='m', &
        length=1, &
        coord_vals=(/350./), &
        cell_bounds=(/0., 700./))
    else if (index(special, 'dzavg2000') > 0) then
      kaxid = cmor_axis( &
        table=trim(tablepath), &
        table_entry='depth2000m', &
        units='m', &
        length=1, &
        coord_vals=(/1000./), &
        cell_bounds=(/0., 2000./))
    else if (trim(vtype) == 'level' .or. vtype(1:4) == 'merd') then
      kaxid = cmor_axis( &
        table=trim(tablepath), &
        table_entry='depth_coord', &
        units='m', &
        length=ddm, &
        coord_vals=depth, &
        cell_bounds=depth_bnds)
    else if (vtype(1:4) == 'merk') then
      kaxid = cmor_axis( &
        table=trim(tablepath), &
        table_entry='rho', &
        units='kg m-3', &
        length=kdm, &
        coord_vals=sigma + 1000., &
        cell_bounds=sigma_bnds + 1000.)
    else if (trim(zcoord) == 'olevel') then
      allocate(tmp1d(1), tmp2d(2, 1))
      tmp1d(:) = (/5.d0/)
      tmp2d(:, 1) = (/0.d0, 10.d0/)
      kaxid = cmor_axis( &
        table=trim(tablepath), &
        table_entry='depth_coord', &
        units='m', &
        length=1, &
        coord_vals=tmp1d, &
        cell_bounds=tmp2d)
      deallocate(tmp1d, tmp2d)
    else if (vtype(1:3) == 'sec') then
      saxid = cmor_axis( &
        table=trim(tablepath), &
        table_entry='oline', &
        units='1', &
        coord_vals=section1)
    end if

    ! Define time axis
    if (.not. fxflag) then
      write(*, *) 'Define time axis '
      write(*, *) 'tablepath:table_entry:', trim(tablepath),':',trim(tcoord)
      write(*, *) 'tcoord:', trim(tcoord)
      write(*, *) 'calunits:', trim(calunits)
      taxid = cmor_axis( &
        table=trim(tablepath), &
        table_entry=trim(tcoord), &
        units=trim(calunits), &
        length=1)
    end if

    ! Define output variable
    write(*, *) 'Define output variable'
    write(*, *) 'zcoord:', trim(zcoord)
    write(*, *) 'vunits:', trim(vunits)
    if (fxflag) then
      if ((trim(vtype) == '2d' .and. .not. (trim(zcoord) == 'olevel' .or. &
        index(special, 'glbave') > 0 .or. index(special, '2zos') > 0.)) &
        .or. lsumz) then
        varid = cmor_variable( &
          table=trim(tablepath), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/grdid/), &
          missing_value=1e20, &
          comment=trim(vcomment), &
          original_name=trim(ivnm))
      else if (index(special, 'glbave') > 0 &
        .or. index(special, '2zos') > 0.) then
        varid = cmor_variable( &
          table=trim(tablepath), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/taxid/), &
          original_name=trim(ivnm), &
          missing_value=1e20, &
          comment=trim(vcomment), &
          positive=trim(vpositive))
      else
        varid = cmor_variable( &
          table=trim(tablepath), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/grdid, kaxid/), &
          missing_value=1e20, &
          comment=trim(vcomment), &
          original_name=trim(ivnm))
      end if
    else
      write(*, *) 'vtype:', trim(vtype)
      write(*, *) 'zcoord:', trim(zcoord)
      if ((trim(vtype) == '2d' .and. .not. (trim(zcoord) == 'ol' .or. &
        index(special, 'glbave') > 0 .or. index(special, '2zos') > 0. &
        )) .or. lsumz .and. .not. index(special, 'glbave') > 0 &
        .or. index(special, 'lvl2srf') > 0 &
        .or. index(special, 'locmin') > 0 &
        .or. index(special, 'dpint') > 0 &
        .or. index(special, 'dpavg') > 0 &
        .or. index(special, 'omega2z') > 0) then
        !write(*, *) 'case l2594'
        !write(*, *) 'vunits:',trim(vunits)
        varid = cmor_variable( &
          table=trim(tablepath), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/grdid, taxid/), &
          original_name=trim(ivnm), &
          missing_value=1e20, &
          comment=trim(vcomment), &
          positive=trim(vpositive))
      else if (index(special, 'dzavg') > 0.) then
        varid = cmor_variable( &
          table=trim(tablepath), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/grdid, taxid, kaxid/), &
          original_name=trim(ivnm), &
          missing_value=1e20, &
          comment=trim(vcomment), &
          positive=trim(vpositive))
      else if (trim(vtype) == 'layer' .and. .not. (trim(ovnm) == 'zfull' &
        .or. trim(ovnm) == 'half' .or. index(special, 'glbave') > 0 &
        .or. index(special, '2zos') > 0.)) then
        write(*, *) 'case l2617'
        varid = cmor_variable( &
          table=trim(tablepath), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/grdid, kaxid, taxid/), &
          original_name=trim(ivnm), &
          missing_value=1e20, &
          positive=trim(vpositive), &
          comment='Please note that the layer depth ' &
          // 'information is stored ' &
          // 'separately in "zfull" ' &
          // 'and "zhalf" while approximate layer ' &
          // 'density values are stored together ' &
          // 'with "msftmrho". '//trim(vcomment))
      else if (vtype(1:4) == 'merd' .or. vtype(1:4) == 'merk') then
        varid = cmor_variable( &
          table=trim(tablepath), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/laxid, kaxid, raxid, taxid/), &
          original_name=trim(ivnm), &
          missing_value=1e20, &
          comment=trim(vcomment), &
          positive=trim(vpositive))
      else if (vtype(1:4) == 'mert') then
        varid = cmor_variable( &
          table=trim(tablepath), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/laxid, raxid, taxid/), &
          original_name=trim(ivnm), &
          missing_value=1e20, &
          comment=trim(vcomment), &
          positive=trim(vpositive))
      else if (vtype(1:4) == 'sect') then
        varid = cmor_variable( &
          table=trim(tablepath), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/saxid, taxid/), &
          original_name=trim(ivnm), &
          missing_value=1e20, &
          comment=trim(vcomment), &
          positive=trim(vpositive))
      else if (index(special, 'glbave') > 0 &
        .or. index(special, '2zos') > 0.) then
        varid = cmor_variable( &
          table=trim(tablepath), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/taxid/), &
          original_name=trim(ivnm), &
          missing_value=1e20, &
          comment=trim(vcomment), &
          positive=trim(vpositive))
      else
        varid = cmor_variable( &
          table=trim(tablepath), &
          table_entry=trim(ovnm), &
          units=trim(vunits), &
          axis_ids=(/grdid, kaxid, taxid/), &
          original_name=trim(ivnm), &
          missing_value=1e20, &
          comment=trim(vcomment), &
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

    status = cmor_close(varid, fnmo, 1)
    if (status /= 0) stop 'problem closing CMOR output file'

  end subroutine close_ofile

  ! -----------------------------------------------------------------

  subroutine read_field

    implicit none

    !real                    :: fac1, fac2, fac3, fac4, fac5, fac6
    integer                 :: i, j, k
    !character(len=slenmax)  :: coord, ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6
    character(len=slenmax)  :: coord

    ! Open input file
    status = nf90_open(fnm, nf90_nowrite, ncid)
    call handle_ncerror(status)

    ! Read data
    if (index(special, 'glbave3d') > 0) then
      fld = 0.
      s1 = 'dp'
      call add_fixed(s1, 1., ncid)
      do k = 1, kk
        do j = 1, jj
          do i = 1, ii
            if (fld(i, j, k) == 1e20) then
              dp(i, j, k) = 0.
            else
              dp(i, j, k) = fld(i, j, k)
            end if
          end do
        end do
      end do
    end if

    !call resolve_vnm(slenmax, ivnm, ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6, &
      !fac1, fac2, fac3, fac4, fac5, fac6)
    fld = 0.

    if (allocated(vars)) then
      do k =1, size(vars)
        call add_fixed(vars(k), facs(k), ncid)
      end do
    else
      call add_fixed(ivnm, 1.0, ncid)
    end if

    if (index(special, 'volcello') > 0) then
      fld2 = fld
      fld = 0.
    end if
    !call add_fixed(ivnm2, fac2, ncid)
    !call add_fixed(ivnm3, fac3, ncid)
    !call add_fixed(ivnm4, fac4, ncid)
    !call add_fixed(ivnm5, fac5, ncid)
    !call add_fixed(ivnm6, fac6, ncid)

    status = nf90_close(ncid)
    call handle_ncerror(status)

  end subroutine read_field

  ! -----------------------------------------------------------------

  subroutine read_tslice(rec, badrec, fname)

    implicit none

    !real                            :: fac1, fac2, fac3, fac4, fac5, fac6
    integer, intent(in)             :: rec
    logical, intent(out)            :: badrec
    character(len=*), intent(in), optional  :: fname
    integer, save                   :: fid
    integer                         :: i, j, k, rec1
    !character(len=slenmax)          :: ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6

    ! Exception for fill day
    rec1 = max(rec, 1)

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
      status = nf90_inq_varid(fid, 'time', rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find time variable'
        stop
      end if
      status = nf90_get_var(fid, rhid, tval, (/rec1/), (/1/))
      call handle_ncerror(status)
      if (rec == 0) tval = tval - 1

      tbnds(1, 1) = max(0., tbnds(1, 1))
      tval = 0.5 * (tbnds(1, 1) + tbnds(2, 1))
    end if

    ! Read data
    if (index(special, 'glbave3d') > 0 .or. &
      (index(special, '2rho') > 0 .and. vtype == 'layer') &
      .or. index(special, '2zos') > 0) then
      fld = 0.
      s1 = 'dp'
      call add_tslice(s1, 1., rec1, fid)
      do k = 1, kk
        do j = 1, jj
          do i = 1, ii
            if (fld(i, j, k) == 1e20) then
              dp(i, j, k) = 0.
            else
              dp(i, j, k) = fld(i, j, k)
            end if
          end do
        end do
      end do
      ! Compute rescaled dp
      if (index(special, '2zoss') > 0) then
        do j = 1, jdm
          do i = 1, idm
            if (pmask(i, j) > 0.5) &
              dp(i, j, :) = dp(i, j, :) * bpini(i, j) / sum(dp(i, j, :))
          end do
        end do
      else if (index(special, '2zost') > 0) then
        do j = 1, jdm
          do i = 1, idm
            if (pmask(i, j) > 0.5) &
              dp(i, j, :) = dp(i, j, :) * bpinit(i, j) / sum(dp(i, j, :))
          end do
        end do
      end if
    end if

    !call resolve_vnm(slenmax, ivnm, ivnm1, ivnm2, ivnm3, ivnm4, ivnm5, ivnm6, &
      !fac1, fac2, fac3, fac4, fac5, fac6)
    if (index(special, '2rho') > 0 .or. index(special, '2zoss') > 0 .or. &
      index(special, 'strmf') > 0 .or. index(special, 'Xfield2') > 0 .or. &
      index(special, 'Dfield2') > 0 .or. index(special, 'dpint') > 0 .or. &
      index(special, 'dpavg') > 0) then
      fld = 0.
      !call add_tslice(ivnm2, fac2, rec1, fid)
      call add_tslice(vars(2), facs(2), rec1, fid)
      fld2 = fld
      fld = 0.
      !call add_tslice(ivnm1, fac1, rec1, fid)
      call add_tslice(vars(1), facs(1), rec1, fid)
      if (index(special, 'strmf') > 0) then
        fldtmp = 0.
        call strmf_eval(idm, jdm, kdm, fld, fld2, fldtmp)
        fld = fldtmp
      end if
    else
      fld = 0.
      if (allocated(vars)) then
        do k = 1, size(vars)
          call add_tslice(vars(k), facs(k), rec1, fid)
        end do
      else
          call add_tslice(ivnm, 1.0, rec1, fid)
      end if
        !call add_tslice(ivnm1, fac1, rec1, fid)
      !call add_tslice(ivnm2, fac2, rec1, fid)
      !call add_tslice(ivnm3, fac3, rec1, fid)
      !call add_tslice(ivnm4, fac4, rec1, fid)
      !call add_tslice(ivnm5, fac5, rec1, fid)
      !call add_tslice(ivnm6, fac6, rec1, fid)
    end if

    ! Read sea level height if necessary
    if (index(special, 'dz2') > 0) then
      status = nf90_inq_varid(fid, 'sealv', rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable sealv '
        stop
      end if
      status = nf90_get_var(fid, rhid, sealv, (/1, 1, rec1/), (/idm, jdm, 1/))
      call handle_ncerror(status)
      status = nf90_get_att(fid, rhid, 'scale_factor', sfac)
      if (status /= nf90_noerr) sfac = 1.
      status = nf90_get_att(fid, rhid, 'add_offset', offs)
      if (status /= nf90_noerr) offs = 0.
      status = nf90_get_att(fid, rhid, '_FillValue', fill)
      do j = 1, jj
        do i = 1, ii
          if (sealv(i, j) == fill) then
            sealv(i, j) = 1e20
          else
            sealv(i, j) = sealv(i, j) * sfac + offs
          end if
        end do
      end do
    end if

    ! Read bottom pressure if necessary
    if (index(special, 'dzavg') > 0) then
      status = nf90_inq_varid(fid, 'pbot', rhid)
      if (status /= nf90_noerr) then
        write(*, *) 'cannot find input variable pbot '
        stop
      end if
      status = nf90_get_var(fid, rhid, pbot, (/1, 1, rec1/), (/idm, jdm, 1/))
      call handle_ncerror(status)
      status = nf90_get_att(fid, rhid, 'scale_factor', sfac)
      if (status /= nf90_noerr) sfac = 1.
      status = nf90_get_att(fid, rhid, 'add_offset', offs)
      if (status /= nf90_noerr) offs = 0.
      status = nf90_get_att(fid, rhid, '_FillValue', fill)
      do j = 1, jj
        do i = 1, ii
          if (pbot(i, j) == fill) then
            pbot(i, j) = 1e20
          else
            pbot(i, j) = pbot(i, j) * sfac + offs
          end if
        end do
      end do
    end if

    status = nf90_close(fid)
    call handle_ncerror(status)

  end subroutine read_tslice

  ! -----------------------------------------------------------------

  subroutine add_tslice(vnm, fac, rec, fid)

    ! Description: add one time slice to output variable fld

    implicit none

    character(len=slenmax), intent(in)  :: vnm
    real, intent(in)                    :: fac
    integer, intent(in)                 :: rec, fid
    integer                             :: i, j, k

    ! Return if variable name is empty
    if (len(trim(vnm)) == 0) return

    ! Read time slice
    status = nf90_inq_varid(fid, trim(vnm), rhid)
    if (status /= nf90_noerr) then
      write(*, *) 'cannot find input variable ', trim(vnm)
      stop
    end if
    if (trim(vtype) == '2d') then
      status = nf90_get_var(fid, rhid, fldtmp, (/1, 1, rec/), (/idm, jdm, 1/))
    else if (trim(vtype) == '1d') then
      status = nf90_get_var(fid, rhid, fldtmp, (/1, 1, rec/), (/1, 1, 1/))
    else if (trim(vtype) == 'layer') then
      status = nf90_get_var(fid, rhid, fldtmp, (/1, 1, 1, rec/), (/idm, jdm, kdm, 1/))
    else if (trim(vtype) == 'level') then
      status = nf90_get_var(fid, rhid, fldtmp, (/1, 1, 1, rec/), (/idm, jdm, ddm, 1/))
    else if (trim(vtype) == 'merk') then
      status = nf90_get_var(fid, rhid, fldtmp, (/1, 1, 1, rec/), (/ldm, kdm, rdm, 1/))
    else if (trim(vtype) == 'merd') then
      status = nf90_get_var(fid, rhid, fldtmp, (/1, 1, 1, rec/), (/ldm, ddm, rdm, 1/))
    else if (trim(vtype) == 'mert') then
      status = nf90_get_var(fid, rhid, fldtmp, (/1, 1, rec/), (/ldm, rdm, 1/))
    else if (trim(vtype) == 'sect') then
      write(*, *) 'read sections ', secdm
      status = nf90_get_var(fid, rhid, fldtmp, (/1, rec/), (/secdm, 1/))
    end if
    call handle_ncerror(status)
    status = nf90_get_att(fid, rhid, 'scale_factor', sfac)
    if (status /= nf90_noerr) sfac = 1.
    status = nf90_get_att(fid, rhid, 'add_offset', offs)
    if (status /= nf90_noerr) offs = 0.
    status = nf90_get_att(fid, rhid, '_FillValue', fill)
    if (status /= nf90_noerr) then
      do k = 1, kk
        do j = 1, jj
          do i = 1, ii
            fld(i, j, k) = fld(i, j, k) + (fldtmp(i, j, k) * sfac + offs) * fac
          end do
        end do
      end do
    else
      do k = 1, kk
        do j = 1, jj
          do i = 1, ii
            if (fldtmp(i, j, k) == fill) then
              fld(i, j, k) = 1e20
            else
              fld(i, j, k) = fld(i, j, k) + (fldtmp(i, j, k) * sfac + offs) * fac
            end if
          end do
        end do
      end do
    end if

  end subroutine add_tslice

  ! -----------------------------------------------------------------

  subroutine add_fixed(vnm, fac, fid)

    ! Description: add one time slice to output variable fld

    implicit none

    character(len=slenmax), intent(in)  :: vnm
    real, intent(in)                    :: fac
    integer, intent(in)                 :: fid
    integer                             :: i, j, k

    ! Return if variable name is empty
    if (len(trim(vnm)) == 0) return

    ! Read time slice
    status = nf90_inq_varid(fid, trim(vnm), rhid)
    if (status /= nf90_noerr) then
      write(*, *) 'cannot find input variable ', trim(vnm)
      stop
    end if
    status = nf90_get_var(fid, rhid, fldtmp)
    call handle_ncerror(status)
    status = nf90_get_att(fid, rhid, 'scale_factor', sfac)
    if (status /= nf90_noerr) sfac = 1.
    status = nf90_get_att(fid, rhid, 'add_offset', offs)
    if (status /= nf90_noerr) offs = 0.
    status = nf90_get_att(fid, rhid, '_FillValue', fill)
    if (status /= nf90_noerr) then
      do k = 1, kk
        do j = 1, jj
          do i = 1, ii
            fld(i, j, k) = fld(i, j, k) + (fldtmp(i, j, k) * sfac + offs) * fac
          end do
        end do
      end do
    else
      do k = 1, kk
        do j = 1, jj
          do i = 1, ii
            if (fldtmp(i, j, k) == fill) then
              fld(i, j, k) = 1e20
            else
              fld(i, j, k) = fld(i, j, k) + (fldtmp(i, j, k) * sfac + offs) * fac
            end if
          end do
        end do
      end do
    end if

  end subroutine add_fixed

  ! -----------------------------------------------------------------

  subroutine write_field

    implicit none

    integer :: i, j, k

    ! Set zero on ocean grid cells
    do k = 1, kk
      do j = 1, jj
        do i = 1, ii
          if (abs(fld(i, j, k)) > 2e20) fld(i, j, k) = 0.
        end do
      end do
    end do

    ! Store variable
    if (index(special, 'glbave') > 0 .or. index(special, '2zos') > 0.) then
      error_flag = cmor_write( &
        var_id=varid, &
        data=(/fld(1, 1, 1)/))
    else if (vtype == '2d') then
      error_flag = cmor_write( &
        var_id=varid, &
        data=reshape(fld, (/idm, jdm/)))
    else
      error_flag = cmor_write( &
        var_id=varid, &
        data=fld)
    end if

  end subroutine write_field

  ! -----------------------------------------------------------------

  subroutine write_tslice

    implicit none

    integer :: i, j, k

    ! Populate field defined at interface level
    if (index(special, 'zhalf') > 0) then
      fldhalf(:, :, 1) = sealv
      fldhalf(:, :, 2:kdm + 1) = fld
      do k = 1, kk + 1
        do j = 1, jj
          do i = 1, ii
            if (abs(fldhalf(i, j, k)) > 2e20) fldhalf(i, j, k) = 0.
          end do
        end do
      end do
    else if (index(special, 'halfl') > 0) then
      fldhalf(:, :, 2:kdm + 1) = fld
      do j = 1, jj
        do i = 1, ii
          if (abs(fldhalf(i, j, 2)) >= 1e20) then
            fldhalf(i, j, 1) = 1e20
          else
            fldhalf(i, j, 1) = 0.
          end if
        end do
      end do
    end if

    ! Set missing on land grid cells
    if (index(special, 'glbave') <= 0) then
      do k = 1, kk
        do j = 1, jj
          do i = 1, ii
            if (abs(fld(i, j, k)) > 1e20) fld(i, j, k) = 1e20
          end do
        end do
      end do
    end if

    ! Store variable
    if (index(special, 'half') > 0) then
      if (trim(tcoord) /= 'time1') then
        error_flag = cmor_write( &
          var_id=varid, &
          data=fldhalf, &
          ntimes_passed=1, &
          time_vals=tval, &
          time_bnds=tbnds)
      else
        error_flag = cmor_write( &
          var_id=varid, &
          data=fldhalf, &
          ntimes_passed=1, &
          time_vals=tval)
      end if
    else
      if (trim(tcoord) == 'time1') then
        error_flag = cmor_write( &
          var_id=varid, &
          data=fld, &
          ntimes_passed=1, &
          time_vals=tval)
      else
        if ((lsumz .or. index(special, 'level1') > 0) .and. &
          .not. index(special, 'glbave') > 0 &
          .or. index(special, 'lvl2srf') > 0 &
          .or. index(special, 'dpint') > 0 &
          .or. index(special, 'dpavg') > 0 &
          .or. index(special, 'locmin') > 0 &
          .or. index(special, 'omega2z') > 0) then
          error_flag = cmor_write( &
            var_id=varid, &
            data=fld(:, :, 1), &
            ntimes_passed=1, &
            time_vals=tval, &
            time_bnds=tbnds)
        else if (index(special, 'glbave') > 0 &
          .or. index(special, '2zos') > 0.) then
          error_flag = cmor_write( &
            var_id=varid, &
            data=(/fld(1, 1, 1)/), &
            ntimes_passed=1, &
            time_vals=tval, &
            time_bnds=tbnds)
        else if (index(special, 'dzavg') > 0) then
          error_flag = cmor_write( &
            var_id=varid, &
            data=(reshape(fld(:, :, 1), (/idm, jdm, 1/))), &
            ntimes_passed=1, &
            time_vals=tval, &
            time_bnds=tbnds)
        else if (vtype(1:4) == 'sect') then
          error_flag = cmor_write( &
            var_id=varid, &
            data=fld(:, 1, 1), &
            ntimes_passed=1, &
            time_vals=tval, &
            time_bnds=tbnds)
        else
          !write(*,*) 'write tslice'
          !write(*,*) 'tval:',tval
          !write(*,*) 'tnbds:',tbnds
          !write(*,*) 'shape fld:',shape(fld)
          error_flag = cmor_write( &
            var_id=varid, &
            data=fld, &
            ntimes_passed=1, &
            time_vals=tval, &
            time_bnds=tbnds)
        end if
      end if
    end if

  end subroutine write_tslice

end module m_modelsocn
