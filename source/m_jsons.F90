module m_jsons

  use json_module
  use m_namelists, only : slenmax

  implicit none

contains

  ! -----------------------------------------------------------------

  subroutine json_get_keys(fnm,path,keys,separator,lfound)

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
    !write(*,*) 'found:',found

    if (associated(entry_ptr)) then
      ! get the number of keys (children)
      call jsonc%info(entry_ptr, n_children=num_children)
      allocate(keys(num_children))
      ! get key names
      do i = 1, num_children
          call jsonc%get_child(entry_ptr, i, child_ptr)
          call jsonc%info(child_ptr, name=key_name)
          keys(i) = key_name
      end do
    !else
        !print *, "Error: key path: '"//trim(path)//"' not found in JSON."
    end if

    call jsonf%destroy()

    if (allocated(keys)) deallocate(keys)
    if (allocated(key_name)) deallocate(key_name)

    if (present(lfound)) lfound = found

  end subroutine json_get_keys

  ! -----------------------------------------------------------------

  subroutine json_get_value_string(fnm, path, str, separator, lfound)

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
    !if (.not. found) write(*,*) "Warning: value of "//trim(path)//" in "//trim(fnm)//" not found"
    call json%destroy()

    if (found) then
      str=cval
      deallocate(cval)
    end if

    if (present(lfound)) lfound = found

  end subroutine json_get_value_string
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
    !if (.not. found) write(*,*) "Warning: value of "//trim(path)//" in "//trim(fnm)//" not found"
    call json%destroy()
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
    !if (.not. found) write(*,*) "Warning: value of "//trim(path)//" in "//trim(fnm)//" not found"
    call json%destroy()
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

    !cnm = ''

    !write(*,*) 'subroutine: get_vertcoord'
    !write(*,*) 'tnm:', trim(tnm)
    !write(*,*) 'vnm:', trim(vnm)
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
    if (allocated(cval)) deallocate(cval)

    if (present(lfound)) lfound = found

  end subroutine json_get_vertcoord

  ! -----------------------------------------------------------------

  subroutine json_get_timecoord(tnm, vnm, str, separator, lfound)

    !use json_module
    !implicit none

    character(len=*), intent(in) :: tnm, vnm
    character(len=*), intent(out) :: str
    logical, intent(out), optional :: lfound
    character(len=*), intent(in), optional :: separator

    type(json_file) :: json
    integer :: k, pdm
    character(len=slenmax), dimension(:), allocatable :: cval
    logical         :: found

    !write(*, *) 'l643, vnm:', trim(vnm)
    !write(*, *) 'l643, tnm:', trim(tnm)

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

  subroutine json_get_units(tnm, vnm, units, separator, lfound)

    !use json_module
    !implicit none

    character(len=*), intent(in) :: tnm, vnm
    character(len=*), intent(out) :: units
    logical, intent(out), optional :: lfound
    character(len=*), intent(in), optional :: separator

    type(json_file) :: json

    character(len=:), allocatable :: cval
    logical   :: found

    if (present(separator)) then
      call json%initialize(path_separator=trim(separator))
    else
      call json%initialize()
    end if
    call json%load_file(filename=trim(tnm))
    call json%get('variable_entry.' // trim(vnm) // '.units', cval, found)
    call json%destroy()
    if (found) then
      units = cval
      deallocate(cval)
    end if

    if (present(lfound)) lfound = found

  end subroutine json_get_units


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
