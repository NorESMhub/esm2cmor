module m_jsons

  use json_module
  use m_namelists, only : slenmax, verbose

  implicit none

contains

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

    call json_get_keys(trim(fnm),'variable_entry:'//trim(vnm)//':preproc', keys, &
        separator=':',lfound=lfound)

  end subroutine json_get_preproc_keys

  ! -----------------------------------------------------------------
  subroutine json_get_preproc_val(tnm, vnm, key, val, lfound)

    character(len=*), intent(in)  :: tnm, vnm, key
    character(len=*), intent(out) :: val
    logical, intent(out), optional :: lfound

    call json_get_val_str(tnm,'variable_entry:'//trim(vnm)//':preproc:'//trim(key), &
        val, separator=':', lfound=lfound)

  end subroutine json_get_preproc_val

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
