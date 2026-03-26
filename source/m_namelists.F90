module m_namelists

  implicit none

  ! Namelist limits
  !integer, parameter :: rowmax = 200, colmax = 3, lenmax = 200
  !integer, parameter :: slenmax = 1024, smax = 10
  integer, parameter :: rowmax = 200, slenmax = 1024

  ! System namelist
  character(len=slenmax), save  :: ibasedir, obasedir, tabledir, griddata
  logical, save                 :: createsubdirs, forcefilescan, verbose
  namelist /sys/ ibasedir, obasedir, tabledir, griddata, createsubdirs, &
    forcefilescan, verbose

  ! Model namelist
  character(len=slenmax), save                  :: model_id, institute_id
  !character(len=slenmax), dimension(smax), save :: institution, source, &
    !references, contact
  character(len=slenmax), save :: institution, source, references, contact
  !character(len=(slenmax+1)*smax), save :: institution1, source1, references1, contact1
  character(len=slenmax), save :: tagoyr, tagoyrbgc, tagomon, tagomonbgc, tagoday, &
    tagodaybgc, tagimon, tagiday, tagamon, tagaday, taga6hr, taga6hri, taga3hr, &
    taga3hri, taglmon, taglday, tagl3hr, tagl3hri, taglyr
  character(len=slenmax), save :: secindexfile, ocngridfile, ocninitfile, ocnmertfile, &
    rhotablesuff, atmgridfile, ocnregnfile
  !logical, save                 :: linebreaks
  character(len=slenmax), save  :: parent_source_id, coordtable, namelist_file_json, &
    atmgrid, atmgrid_label, atmgrid_resolution, ocngrid, ocngrid_label, ocngrid_resolution, &
    icegrid, icegrid_label, icegrid_resolution, lndgrid, lndgrid_label, lndgrid_resolution

  namelist /model/  model_id, institute_id, &
                    institution, source, references, contact, &
                    tagoyr, tagoyrbgc, tagomon, tagomonbgc, tagoday, tagodaybgc, &
                    tagimon, tagiday, tagamon, tagaday, taga6hr, taga6hri, taga3hr, &
                    taga3hri, taglmon, taglday, tagl3hr, tagl3hri, taglyr, &
                    secindexfile, ocngridfile, ocninitfile, ocnmertfile, &
                    rhotablesuff, atmgridfile, ocnregnfile, &
                    !linebreaks, &
                    parent_source_id, coordtable, namelist_file_json, &
                    atmgrid, atmgrid_label, atmgrid_resolution, &
                    ocngrid, ocngrid_label, ocngrid_resolution, &
                    icegrid, icegrid_label, icegrid_resolution, &
                    lndgrid, lndgrid_label, lndgrid_resolution

  ! Experiment namelist
  character(len=slenmax), save :: casename, experiment_id, parent_experiment_id, &
    parent_experiment_rip, isubdir, osubdir, membertag
  !character(len=slenmax), dimension(smax), save :: history, comment, forcing
  character(len=slenmax), save :: history, comment, forcing
  integer, save                                 :: realization, exprefyear, year1, yearn, month1, monthn
  real(kind=8), save                            :: branch_time
  logical, save :: dry_run, plevdummy, readdummy, add_fill_day, scanallfiles
  integer, save :: physics_version = 1, initialization_method = 1
  character(len=slenmax), save  :: activity_id, parent_variant_label, &
    parent_mip_era, mip_era, sub_experiment_id, parent_sub_experiment, &
    parent_activity_id, branch_method, parent_time_units, tracking_prefix, &
    variant_label, source_type
  real(kind=8), save            :: branch_time_in_child, branch_time_in_parent
  character(len=slenmax), save  :: forcing_index, physics_index, realization_index, &
    initialization_index
  namelist /experiment/ casename, experiment_id, parent_experiment_id, &
                        parent_experiment_rip, isubdir, osubdir, membertag, &
                        history, comment, forcing, &
                        realization, exprefyear, year1, yearn, month1, monthn, &
                        branch_time, &
                        dry_run, plevdummy, readdummy, add_fill_day, scanallfiles, &
                        activity_id, parent_variant_label, parent_mip_era, mip_era, &
                        sub_experiment_id, parent_sub_experiment, parent_activity_id, branch_method, &
                        parent_time_units, tracking_prefix, variant_label, source_type, &
                        branch_time_in_child, branch_time_in_parent, &
                        forcing_index, physics_index, realization_index, initialization_index

  ! Variables
  character(len=slenmax), save  :: pomon
  integer, save                 :: n_variables
  character(len=slenmax), dimension(rowmax), save :: compound_names, branded_names, &
    realms, frequencies, regions
  character(len=slenmax) :: realm, frequency

  namelist /variables/ compound_names

  ! Misc
  integer :: istatus, funit
  character(len=slenmax), save :: fnm, fnm2, itag, nmlfp
  character(len=slenmax), save :: nmlfpsys, nmlfpmod, nmlfpexp, nmlfpvar

  ! Time related variables
  logical, save :: linstant
  integer, save :: year, month, rec
  real(kind=8), save :: tval(1), tval2(2), tbnd(2), mbnd(2), tbnds(2,1), mbnds(2,1)
  character(len=slenmax), save :: calendar = 'noleap', &
    calunits = 'days since 1850-01-01 00:00:00'

contains

  ! -----------------------------------------------------------------

  subroutine read_namelists

    implicit none

    integer :: n, pos, nstr
    logical :: fexist
    character(len=slenmax)  :: substr, tmpstr

    ! Initialise namelist variables
    ibasedir      = ' '
    obasedir      = ' '
    tabledir      = ' '
    griddata      = ' '
    tagoyr        = ' '
    tagoyrbgc     = ' '
    tagomon       = ' '
    tagomonbgc    = ' '
    tagoday       = ' '
    tagodaybgc    = ' '
    tagimon       = ' '
    tagiday       = ' '
    tagamon       = ' '
    tagaday       = ' '
    taga6hr       = ' '
    taga6hri      = ' '
    taga3hr       = ' '
    taga3hri      = ' '
    taglyr        = ' '
    taglmon       = ' '
    taglday       = ' '
    tagl3hr       = ' '
    tagl3hri      = ' '
    atmgridfile   = ' '
    ocngridfile   = ' '
    ocninitfile   = ' '
    ocnmertfile   = ' '
    ocnregnfile   = ' '
    secindexfile  = ' '
    rhotablesuff  = 'OnRho'
    coordtable    = 'CMIP7_coordinate.json'
    year1         = 0
    month1        = 1
    yearn         = 0
    monthn        = 12
    createsubdirs = .true.
    forcefilescan = .true.
    verbose       = .true.
    dry_run       = .false.
    plevdummy     = .false.
    readdummy     = .false.
    add_fill_day  = .false.
    !newcolumnorder= .true.
    scanallfiles  = .true.

    casename      = ' '
    experiment_id = ' '
    institute_id  = ' '
    institution   = ' '
    source        = ' '
    contact       = ' '
    history       = ' '
    comment       = ' '
    references    = ' '
    model_id      = ' '
    forcing       = ' '
    realization   = 1
    branch_time   = 0.0
    parent_experiment_id = ' '
    parent_experiment_rip = ' '
    isubdir       = ' '
    osubdir       = ' '
    membertag     = ' '

    compound_names = ''

    pomon = ''


    ! Read namelists
    if (iargc() /= 4) then
      write(*, *) 'Usage: noresm2cmor <system nml-file> <model nml-file>' // &
        '<exp nml-file> <variable nml-file>'
      stop
    end if
    !else if (iargc() == 4 .or. iargc() == 5) then
    nmlfpsys = ' '
    call getarg(1, nmlfpsys)
    inquire(file=trim(nmlfpsys), exist=fexist)
    if (.not. fexist) then
      write(*, *) 'cannot find namelist file: '//trim(nmlfpsys)
      stop
    end if
    nmlfpmod = ' '
    call getarg(2, nmlfpmod)
    inquire(file=trim(nmlfpmod), exist=fexist)
    if (.not. fexist) then
      write(*, *) 'cannot find namelist file: '//trim(nmlfpmod)
      stop
    end if
    nmlfpexp = ' '
    call getarg(3, nmlfpexp)
    inquire(file=trim(nmlfpexp), exist=fexist)
    if (.not. fexist) then
      write(*, *) 'cannot find namelist file: '//trim(nmlfpexp)
      stop
    end if
    nmlfpvar = ' '
    call getarg(4, nmlfpvar)
    inquire(file=trim(nmlfpvar), exist=fexist)
    if (.not. fexist) then
      write(*, *) 'cannot find namelist file: '//trim(nmlfpvar)
      stop
    end if
      !if (iargc() == 5) then
        !call getarg(5, vsingle)
      !end if
    !end if

    funit = get_free_unit()
    open(funit, file=trim(nmlfpsys), status='old', action='read', recl=200)
    read(funit, nml=sys, iostat=istatus)
    close(funit)
    if (istatus /= 0) then
      write(*, *) 'Problem reading namelist system in file '//trim(nmlfpsys)
      stop
    end if
    griddata = trim(griddata)//'/'
    tabledir = trim(tabledir)//'/'

    open(funit, file=trim(nmlfpmod), status='old', action='read', recl=200)
    read(funit, nml=model, iostat=istatus)
    close(funit)
    if (istatus /= 0) then
      write(*, *) 'Problem reading namelist model in file '//trim(nmlfpmod)
      stop
    end if

    open(funit, file=trim(nmlfpexp), status='old', action='read', recl=200)
    read(funit, nml=experiment, iostat=istatus)
    close(funit)
    if (istatus /= 0) then
      write(*, *) 'Problem reading namelist experiment in file '//trim(nmlfpexp)
      stop
    end if

    open(funit, file=trim(nmlfpvar), status='old', action='read', recl=200)
    read(funit, nml=variables, iostat=istatus)
    !rewind(funit)
    if (istatus /= 0) write(*, *) &
      'WARNING: compound_names not in namelist file.'
    close(funit)

    n_variables = 0
    do n = 1, rowmax
      nstr = len_trim(compound_names(n))
      if (nstr /= 0) then
        substr = compound_names(n)
        pos = index(substr,'.')
        realms(n) = substr(1:pos-1)
        substr = substr(pos+1:nstr)
        pos = index(substr,'.')
        tmpstr = substr(1:pos-1)
        substr = substr(pos+1:nstr)
        pos = index(substr,'.')
        branded_names(n) = trim(tmpstr)//'_'//substr(1:pos-1)
        substr = substr(pos+1:nstr)
        pos = index(substr,'.')
        frequencies(n) = substr(1:pos-1)
        regions(n) = substr(pos+1:nstr)
        n_variables = n_variables + 1
      else
        exit
      end if
    end do

    ! Extend input path if necessary
    if (len_trim(isubdir) > 0) &
      ibasedir = trim(ibasedir)//'/'//trim(isubdir)

    ! Modify output path and create output folder
    obasedir = trim(obasedir)//'/'//trim(osubdir)
    write(*,*) 'obsedir:',trim(obasedir)
    call system('mkdir -p '//trim(obasedir))

  end subroutine read_namelists

  ! -----------------------------------------------------------------

  subroutine print_namelists

    implicit none

    integer :: n

    write(*, *)
    write(*, *) 'System namelist:'
    write(*, *) ' input directory  = ', trim(ibasedir)
    write(*, *) ' output directory = ', trim(obasedir)
    write(*, *) ' table directory  = ', trim(tabledir)
    write(*, *) ' grid data dir.   = ', trim(griddata)
    write(*, *) ' create sub-dirs  = ', createsubdirs
    write(*, *) ' verbose          = ', verbose
    write(*, *)
    write(*, *) 'Model namelist:'
    write(*, *) ' institution      = ', trim(institution)
    write(*, *) ' model id         = ', trim(model_id)
    write(*, *) ' source           = ', trim(source)
    write(*, *) ' references       = ', trim(references)
    write(*, *) ' contact          = ', trim(contact)
    write(*, *) ' tag annual ocn   = ', trim(tagoyr)
    write(*, *) ' tag annual bgc   = ', trim(tagoyrbgc)
    write(*, *) ' tag monthly ocn  = ', trim(tagomon)
    write(*, *) ' tag monthly bgc  = ', trim(tagomonbgc)
    write(*, *) ' tag daily ocn    = ', trim(tagoday)
    write(*, *) ' tag daily bgc    = ', trim(tagodaybgc)
    write(*, *) ' tag monthly ice  = ', trim(tagimon)
    write(*, *) ' tag daily ice    = ', trim(tagiday)
    write(*, *) ' tag monthly atm  = ', trim(tagamon)
    write(*, *) ' tag daily atm    = ', trim(tagaday)
    write(*, *) ' tag 6hourly atm  = ', trim(taga6hr)
    write(*, *) ' tag 6hourly insa = ', trim(taga6hri)
    write(*, *) ' tag 3hourly atm  = ', trim(taga3hr)
    write(*, *) ' tag 3hourly insa = ', trim(taga3hri)
    write(*, *) ' tag yearly lnd   = ', trim(taglyr)
    write(*, *) ' tag monthly lnd  = ', trim(taglmon)
    write(*, *) ' tag daily lnd    = ', trim(taglday)
    write(*, *) ' tag 3hourly lnd  = ', trim(tagl3hr)
    write(*, *) ' tag 3hourly insl = ', trim(tagl3hri)
    write(*, *) ' atmos grid file  = ', trim(atmgridfile)
    write(*, *) ' ocean grid file  = ', trim(ocngridfile)
    write(*, *) ' ocean ini file   = ', trim(ocninitfile)
    write(*, *) ' ocean sec file   = ', trim(secindexfile)
    write(*, *) ' ocean moc file   = ', trim(ocnmertfile)
    write(*, *) ' ocean reg file   = ', trim(ocnregnfile)
    !write(*, *) ' allow line break = ', linebreaks

    write(*, *)
    write(*, *) 'Experiment namelist:'
    write(*, *) ' case name        = ', trim(casename)
    write(*, *) ' experiment id    = ', trim(experiment_id)
    write(*, *) ' history          = ', trim(history)
    write(*, *) ' comment          = ', trim(comment)
    write(*, *) ' forcing          = ', trim(forcing)
    write(*, *) ' realization      = ', realization
    write(*, *) ' start year       = ', year1
    write(*, *) ' end year         = ', yearn
    write(*, *) ' start month      = ', month1
    write(*, *) ' end month        = ', monthn
    write(*, *) ' add dummy day    = ', add_fill_day
    write(*, *) ' dry run          = ', dry_run

    write(*, *)
    print *, 'Variable list:'
    do n = 1, n_variables
      print *, trim(compound_names(n))
    end do

  end subroutine print_namelists

  ! -----------------------------------------------------------------

  subroutine merge_strarr(slen, sdm, strin, strout, lb)

    implicit none

    integer, intent(in) :: sdm, slen
    character(len=slen), dimension(sdm), intent(in) :: strin
    character(len=(slen+1)*sdm), intent(out) :: strout
    logical, intent(in) :: lb

    integer :: n, pos

    strout = ' '
    pos = 0
    do n = 1, sdm
      if (len_trim(strin(n)) > 0) then
        if (pos /= 0) then
          pos = pos + 1
          if (lb) then
            strout(pos:pos) = achar(10)
          else
            strout(pos:pos) = ' '
          end if
        end if
        strout(pos+1:pos+len_trim(strin(n))) = trim(strin(n))
        pos = pos + len_trim(strin(n))
      end if
    end do
    pos = 1
    do
      if (index(strout(pos:), '\n') > 0) then
        pos = pos + index(strout(pos:), '\n')
        strout(pos-1:pos) = ' '//achar(10)
        pos = pos + 2
        if (pos >= len(strout)) exit
      else
        exit
      end if
    end do

  end subroutine merge_strarr

  ! -----------------------------------------------------------------

  function get_free_unit() result(free_unit)

    implicit none

    integer :: free_unit
    logical :: in_use

    do free_unit = 10, 99
      inquire(free_unit, opened=in_use)
      if (.not. in_use) exit
    end do

  end function get_free_unit

  ! -----------------------------------------------------------------

  subroutine write_namelist_json(grid, grid_label, grid_resolution, varname)

    use json_module

    implicit none

#ifdef MPI
    include 'mpif.h'
    integer :: mpirank, mpisize, mpierror
#endif

    character(len=*), intent(in) :: grid, grid_label, grid_resolution, varname
    character :: yyyymm1*6, yyyymm2*6, c2*2, r3*3
    type(json_core) :: json
    type(json_value), pointer :: p


    !write(*,*) 'varname:',trim(varname)
    call json%initialize()
    call json%create_object(p, '')

    ! mapped from cmor2
    call json%add(p, 'outpath', trim(obasedir))
    call json%add(p, 'experiment_id', trim(experiment_id))
    call json%add(p, 'institution_id', trim(institute_id))
    call json%add(p, 'institution', trim(institution))
    call json%add(p, 'source_id', trim(model_id))
    call json%add(p, 'source', trim(source))
    call json%add(p, 'calendar', trim(calendar))
    call json%add(p, 'realization_index', trim(realization_index))
    call json%add(p, 'physics_index', trim(physics_index))
    call json%add(p, 'initialization_index', trim(initialization_index))
    call json%add(p, 'contact', trim(contact))
    call json%add(p, 'history', trim(history))
    call json%add(p, 'comment', trim(comment))
    call json%add(p, 'references', trim(references))
    call json%add(p, 'model_id', trim(model_id))
    call json%add(p, 'run_variant', trim(forcing))
    call json%add(p, 'branch_time', branch_time)
    call json%add(p, 'parent_experiment_id', trim(parent_experiment_id))
    ! new for cmor3
    call json%add(p, 'forcing_index', trim(forcing_index))
    call json%add(p, 'parent_variant_label', trim(parent_variant_label))
    call json%add(p, '_controlled_vocabulary_file', 'cmor-cvs.json')
    call json%add(p, '_AXIS_ENTRY_FILE', 'CMIP7_coordinate.json')
    call json%add(p, '_FORMULA_VAR_FILE', 'CMIP7_formula_terms.json')
    call json%add(p, '_cmip7_option', 1)

    ! required global atrtributes
    call json%add(p, 'activity_id', trim(activity_id))
    !!call json%add(p, 'area_label', trim(area_label))
    call json%add(p, 'source_type', trim(source_type))
    call json%add(p, 'sub_experiment_id', trim(sub_experiment_id))
    call json%add(p, 'parent_sub_experiment_id', trim(parent_sub_experiment))
    call json%add(p, 'parent_mip_era', trim(parent_mip_era))
    call json%add(p, 'mip_era', trim(mip_era))
    call json%add(p, 'parent_activity_id', trim(parent_activity_id))
    call json%add(p, 'parent_source_id', trim(parent_source_id))
    call json%add(p, 'grid', trim(grid))
    call json%add(p, 'grid_label', trim(grid_label))
    call json%add(p, 'nominal_resolution', trim(grid_resolution))
    call json%add(p, 'branch_method', trim(branch_method))
    call json%add(p, 'branch_time_in_child', branch_time_in_child)
    call json%add(p, 'branch_time_in_parent', branch_time_in_parent)
    call json%add(p, 'parent_time_units', trim(parent_time_units))
    call json%add(p, 'tracking_prefix', trim(tracking_prefix))
    call json%add(p, 'output_path_template', &
      '<activity_id><institution_id><source_id><experiment_id><_member_id><table><variable_id><grid_label><version>')
    call json%add(p, 'output_file_template', &
      '<variable_id><table><source_id><experiment_id><_member_id><grid_label>')
    call json%add(p, 'license_id', 'CC-BY-4.0')
    call json%add(p, 'archive_id', 'WCRP')

    call json%add(p, 'frequency', trim(frequency))
    call json%add(p, 'region', 'glb')
    call json%add(p, 'branded_variable', trim(varname))
    call json%add(p, 'drs_specs', 'MIP-DRS7')

    write(yyyymm1, '(I4.4,I2.2)') year1, month1
    write(yyyymm2, '(I4.4,I2.2)') yearn, monthn
    write(r3, '(I3.3)') realization
#ifdef MPI
    call mpi_comm_rank(mpi_comm_world, mpirank, mpierror)
    write(c2, '(I2.2)') mpirank
#else
    c2 = '00'
#endif

    namelist_file_json = 'namelist_'//trim(casename)//'_'// &
      trim(varname)//'_'//yyyymm1//'-'//yyyymm2//'_'//r3//'_'//c2//'.json'
    call json%print(p, namelist_file_json)
    call json%destroy(p)            ! cleanup
    !write(*, *) 'json-namelist file: ', trim(namelist_file_json)

  end subroutine write_namelist_json

end module m_namelists
