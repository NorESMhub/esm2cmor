module m_jsons

  use json_module
  !use m_namelists, only : slenmax, verbose
  use m_namelists

  implicit none

contains

  ! -----------------------------------------------------------------
  subroutine json_write_attributes(grid, grid_label, grid_resolution, varname)

    !use json_module

    !implicit none

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

    json_file_attributes = 'attributes_'//trim(casename)//'_'// &
      trim(varname)//'_'//yyyymm1//'-'//yyyymm2//'_'//r3//'_'//c2//'.json'
    call json%print(p, json_file_attributes)
    call json%destroy(p)            ! cleanup
    !write(*, *) 'json-namelist file: ', trim(json_file_attributes)

  end subroutine json_write_attributes

  ! -----------------------------------------------------------------
  subroutine json_get_keys(fnm, path, keys, separator, lfound)

    character(len=*),                                  intent(in)  :: fnm, path
    character(len=slenmax), dimension(:), allocatable, intent(out) :: keys
    character(len=*), intent(in), optional :: separator
    logical, intent(out), optional :: lfound

    type(json_file) :: jsonf
    type(json_core) :: jsonc
    type(json_value), pointer :: entry_ptr => null()
    type(json_value), pointer :: child_ptr => null()

    integer :: i, num_children
    logical :: found

    character(len=:), allocatable :: key_name

    if (present(separator)) then
      call jsonf%initialize(path_separator=trim(separator))
    else
      call jsonf%initialize()
    end if

    call jsonf%load_file(filename=trim(fnm))
    if (jsonf%failed()) stop "Error: Could not load JSON file: "//fnm
    call jsonf%get(path, entry_ptr,found)

    if (associated(entry_ptr)) then
      call jsonc%info(entry_ptr, n_children=num_children)
      allocate(keys(num_children))
      do i = 1, num_children
          call jsonc%get_child(entry_ptr, i, child_ptr)
          call jsonc%info(child_ptr, name=key_name)
          keys(i) = key_name
      end do
!   else
!     if (verbose) write(*,*) "Warning: '"//trim(path)//"' not found in JSON."
    end if

    call jsonf%destroy()

    !if (allocated(keys)) deallocate(keys)
    if (allocated(key_name)) deallocate(key_name)

    if (present(lfound)) lfound = found

  end subroutine json_get_keys

  ! -----------------------------------------------------------------
  subroutine json_get_val_str(fnm, path, str, separator, lfound)

    character(len=*),              intent(in)  :: fnm, path
    character(len=*), intent(out) :: str
    logical, intent(out), optional :: lfound
    character(len=*), intent(in), optional :: separator

    type(json_file) :: json

    character(len=:), allocatable :: cval
    logical :: found

    if (present(separator)) then
      call json%initialize(path_separator=trim(separator))
    else
      call json%initialize()
    end if

    call json%load_file(filename=trim(fnm))
    if (json%failed()) stop "Error: Could not load JSON file: "//trim(fnm)
    call json%get(path, cval, found)
    call json%destroy()
    !if (.not. found) write(*,*) "Warning: value of "//trim(path)//" in "//trim(fnm)//" not found"

    if (found) then
      str=cval
      deallocate(cval)
    end if

    if (present(lfound)) lfound = found

  end subroutine json_get_val_str

  ! -----------------------------------------------------------------
  subroutine json_get_array_string(fnm,path,array,separator,lfound)

    character(len=*),                                  intent(in)  :: fnm, path
    character(len=slenmax), dimension(:), allocatable, intent(out) :: array
    logical, intent(out), optional :: lfound
    character(len=*), intent(in), optional :: separator

    type(json_file) :: json

    logical :: found

    if (present(separator)) then
      call json%initialize(path_separator=trim(separator))
    else
      call json%initialize()
    end if
    call json%load_file(filename=trim(fnm))
    if (json%failed()) stop "Error: Could not load JSON file: "//fnm
    call json%get(path, array, found)
    call json%destroy()
    !if (.not. found) write(*,*) "Warning: value of "//trim(path)//" in "//trim(fnm)//" not found"

    !if (found) deallocate(array)
    if (present(lfound)) lfound = found

  end subroutine json_get_array_string

  ! -----------------------------------------------------------------
  subroutine json_get_array_real(fnm,path,array,separator,lfound)

    character(len=*),                        intent(in)  :: fnm, path
    real(kind=8), dimension(:), allocatable, intent(out) :: array
    logical, intent(out), optional :: lfound
    character(len=*), intent(in), optional :: separator

    type(json_file) :: json

    logical :: found

    if (present(separator)) then
      call json%initialize(path_separator=trim(separator))
    else
      call json%initialize()
    end if
    call json%load_file(filename=trim(fnm))
    if (json%failed()) stop "Error: Could not load JSON file: "//fnm
    call json%get(path, array, found)
    call json%destroy()
    !if (.not. found) write(*,*) "Warning: value of "//trim(path)//" in "//trim(fnm)//" not found"

    !if (found) deallocate(array)
    if (present(lfound)) lfound = found

  end subroutine json_get_array_real

  ! -----------------------------------------------------------------
  subroutine json_get_vertcoord(tnm, vnm, str, separator, lfound)

    character(len=*), intent(in) :: tnm, vnm
    character(len=*), intent(out) :: str
    logical, intent(out), optional :: lfound
    character(len=*), intent(in), optional :: separator

    type(json_file) :: json
    integer :: k, pdm
    character(len=slenmax), dimension(:), allocatable :: cval
    logical :: found

    if (present(separator)) then
      call json%initialize(path_separator=trim(separator))
    else
      call json%initialize()
    end if
    call json%load_file(filename=trim(tnm))
    call json%get('variable_entry.' // trim(vnm) // '.dimensions', cval, found)
    call json%destroy()
      

    if (.not. found) then
      if (present(lfound)) lfound = found
      return
    end if

    found = .false.
    pdm = size(cval)
    do k = 1, pdm
      if (cval(k) == 'alevel' .or. cval(k)(1:4) == 'plev' .or. &
          cval(k) == 'olevel') then
        str = cval(k)
        found = .true.
        exit
      end if
    end do

    if (present(lfound)) lfound = found

    if (allocated(cval)) deallocate(cval)

  end subroutine json_get_vertcoord

  ! -----------------------------------------------------------------
  subroutine json_get_timecoord(tnm, vnm, str, separator, lfound)

    character(len=*), intent(in) :: tnm, vnm
    character(len=*), intent(out) :: str
    logical, intent(out), optional :: lfound
    character(len=*), intent(in), optional :: separator

    type(json_file) :: json
    integer :: k, pdm
    character(len=slenmax), dimension(:), allocatable :: cval
    logical         :: found

    if (present(separator)) then
      call json%initialize(path_separator=trim(separator))
    else
      call json%initialize()
    end if
    call json%load_file(filename=trim(tnm))
    call json%get('variable_entry.' // trim(vnm) // '.dimensions', cval, found)
    call json%destroy()
    
    if (.not. found) then
      if (present(lfound)) lfound = found
      return
    end if

    found = .false.
    pdm = size(cval)
    do k = 1, pdm
      if (cval(k)(1:4) == 'time') then
        str = cval(k)
        found = .true.
        exit
      end if
    end do

    if(present(lfound)) lfound = found
    if(allocated(cval)) deallocate(cval)

  end subroutine json_get_timecoord

  ! -----------------------------------------------------------------
  subroutine json_get_original_name(tnm, vnm, str, lfound)

    character(len=*),              intent(in)  :: tnm, vnm
    character(len=*), intent(out) :: str
    logical, intent(out), optional :: lfound
    !character(len=*), intent(in), optional :: separator

    call json_get_val_str(trim(tnm), 'variable_entry:'//trim(vnm)// &
            ':original_name', str, separator=':', lfound=lfound)

  end subroutine json_get_original_name

  ! -----------------------------------------------------------------
  subroutine json_get_units(tnm, vnm, str, lfound)

    character(len=*),              intent(in)  :: tnm, vnm
    character(len=*), intent(out) :: str
    logical, intent(out), optional :: lfound
    !character(len=*), intent(in), optional :: separator

    call json_get_val_str(trim(tnm), 'variable_entry.'//trim(vnm)// &
            '.units', str, separator='.', lfound=lfound)

  end subroutine json_get_units

! -----------------------------------------------------------------
  subroutine json_get_vars(fnm,vnm,array,lfound)

    character(len=*),                                  intent(in)  :: fnm, vnm
    character(len=slenmax), dimension(:), allocatable, intent(out) :: array
    logical, intent(out), optional :: lfound

    call json_get_array_string(trim(fnm), 'variable_entry:'//trim(vnm)// &
            ':sources:vars', array, separator=':', lfound=lfound)

  end subroutine json_get_vars

! -----------------------------------------------------------------
  subroutine json_get_facs(fnm,vnm,array,lfound)

    character(len=*),                        intent(in)  :: fnm, vnm
    !character(len=slenmax), dimension(:), allocatable, intent(out) :: array
    real(kind=8), dimension(:), allocatable, intent(out) :: array
    logical, intent(out), optional :: lfound

    call json_get_array_real(trim(fnm), 'variable_entry:'//trim(vnm)// &
            ':sources:facs', array, separator=':', lfound=lfound)

  end subroutine json_get_facs

! -----------------------------------------------------------------
  subroutine json_get_preproc_keys(fnm, vnm, keys, lfound)

    character(len=*),                                  intent(in)  :: fnm, vnm
    character(len=slenmax), dimension(:), allocatable, intent(out) :: keys
    logical, intent(out), optional :: lfound

    call json_get_keys(trim(fnm),'variable_entry:'//trim(vnm)//':preprocs', keys, &
        separator=':',lfound=lfound)

  end subroutine json_get_preproc_keys

  ! -----------------------------------------------------------------
  subroutine json_get_preproc_val(tnm, vnm, key, val, lfound)

    character(len=*), intent(in)  :: tnm, vnm, key
    character(len=*), intent(out) :: val
    logical, intent(out), optional :: lfound

    call json_get_val_str(tnm,'variable_entry:'//trim(vnm)//':preprocs:' &
        //trim(key), val, separator=':', lfound=lfound)

  end subroutine json_get_preproc_val

! -----------------------------------------------------------------
  subroutine json_get_postproc_keys(fnm, vnm, keys, lfound)

    character(len=*),                                  intent(in)  :: fnm, vnm
    character(len=slenmax), dimension(:), allocatable, intent(out) :: keys
    logical, intent(out), optional :: lfound

    call json_get_keys(trim(fnm),'variable_entry:'//trim(vnm)//':postprocs', keys, &
        separator=':',lfound=lfound)

  end subroutine json_get_postproc_keys

  ! -----------------------------------------------------------------
  subroutine json_get_postproc_val(tnm, vnm, key, val, lfound)

    character(len=*), intent(in)  :: tnm, vnm, key
    character(len=*), intent(out) :: val
    logical, intent(out), optional :: lfound

    call json_get_val_str(tnm,'variable_entry:'//trim(vnm)//':postprocs:' &
        //trim(key), val, separator=':', lfound=lfound)

  end subroutine json_get_postproc_val

! ! -----------------------------------------------------------------

! subroutine json_get_variable()
!   real :: r
!   integer :: ind

! end subroutine json_get_variable

! ! -----------------------------------------------------------------

! subroutine json_get_cmor_table()
!   real :: r
!   integer :: ind

! end subroutine json_get_cmor_table

! ! -----------------------------------------------------------------

! subroutine json_get_cmor_coordinate()
!   real :: r
!   integer :: ind

! end subroutine json_get_cmor_coordinate

! ! -----------------------------------------------------------------

! subroutine json_get_cmor_grid()
!   real :: r
!   integer :: ind

! end subroutine json_get_cmor_grid

! ! -----------------------------------------------------------------

! subroutine json_get_cmor_cv()
!   real :: r
!   integer :: ind

! end subroutine json_get_cmor_cv

  ! -----------------------------------------------------------------

end module m_jsons
