module m_jsons

  use json_module
  use m_namelists, only : slenmax

  implicit none

contains

  ! -----------------------------------------------------------------

  subroutine json_get_keys(fnm,path,keys,found)

    character(len=*),                                  intent(in)  :: fnm, path
    character(len=slenmax), dimension(:), allocatable, intent(out) :: keys
    logical, intent(out), optional :: found

    type(json_file) :: jsonf
    type(json_core) :: jsonc
    type(json_value), pointer :: entry_ptr => null()
    type(json_value), pointer :: child_ptr => null()

    integer :: i, num_children

    character(len=:), allocatable :: key_name

    call jsonf%initialize(path_separator=':')
    call jsonf%load_file(filename=trim(fnm))
    if (jsonf%failed()) stop "Error: Could not load JSON file: "//fnm
    call jsonf%get(path, entry_ptr,found)
    write(*,*) 'found:',found

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

  end subroutine json_get_keys

  ! -----------------------------------------------------------------

  subroutine json_get_value(fnm,path,str,found)

    character(len=*),              intent(in)  :: fnm, path
    character(len=*), intent(out) :: str
    logical, intent(out), optional :: found

    type(json_file) :: jsonf

    character(len=:), allocatable :: cval


    call jsonf%initialize(path_separator=':')
    call jsonf%load_file(filename=trim(fnm))
    if (jsonf%failed()) stop "Error: Could not load JSON file: "//trim(fnm)
    call jsonf%get(path, cval, found)
    !if (.not. found) write(*,*) "Warning: value of "//trim(path)//" in "//trim(fnm)//" not found"
    call jsonf%destroy()
    if (found) then
      str=cval
      deallocate(cval)
    end if

  end subroutine json_get_value
  ! -----------------------------------------------------------------

  subroutine json_get_array_string(fnm,path,array,found)

    character(len=*),                                  intent(in)  :: fnm, path
    character(len=slenmax), dimension(:), allocatable, intent(out) :: array
    logical, intent(out), optional :: found

    type(json_file) :: jsonf

    !logical :: found

    call jsonf%initialize(path_separator=':')
    call jsonf%load_file(filename=trim(fnm))
    if (jsonf%failed()) stop "Error: Could not load JSON file: "//fnm
    call jsonf%get(path, array, found)
    !if (.not. found) write(*,*) "Warning: value of "//trim(path)//" in "//trim(fnm)//" not found"
    call jsonf%destroy()
    if (found) deallocate(array)

  end subroutine json_get_array_string

  ! -----------------------------------------------------------------

  subroutine json_get_array_real(fnm,path,array,found)

    character(len=*),                        intent(in)  :: fnm, path
    real(kind=8), dimension(:), allocatable, intent(out) :: array
    logical, intent(out), optional :: found

    type(json_file) :: jsonf

    !logical :: found

    call jsonf%initialize(path_separator=':')
    call jsonf%load_file(filename=trim(fnm))
    if (jsonf%failed()) stop "Error: Could not load JSON file: "//fnm
    call jsonf%get(path, array, found)
    !if (.not. found) write(*,*) "Warning: value of "//trim(path)//" in "//trim(fnm)//" not found"
    call jsonf%destroy()
    if (found) deallocate(array)

  end subroutine json_get_array_real

  ! -----------------------------------------------------------------

  subroutine json_get_vertcoord(tnm, vnm, str, found)

    use json_module
    implicit none

    character(len=*), intent(in) :: tnm, vnm
    character(len=*), intent(out) :: str
    logical, intent(out), optional :: found

    type(json_file) :: json
    integer :: k, pdm
    character(len=slenmax), dimension(:), allocatable :: cval

    !cnm = ''

    !write(*,*) 'subroutine: get_vertcoord'
    !write(*,*) 'tnm:', trim(tnm)
    !write(*,*) 'vnm:', trim(vnm)
    call json%initialize()
    call json%load_file(filename=trim(tnm))
    call json%get('variable_entry.' // trim(vnm) // '.dimensions', cval, found)
    call json%destroy()

    if (.not. found) return

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
    deallocate(cval)

  end subroutine json_get_vertcoord

  ! -----------------------------------------------------------------

  subroutine json_get_timecoord(tnm, vnm, str, found)

    use json_module
    implicit none

    character(len=*), intent(in) :: tnm, vnm
    character(len=*), intent(out) :: str
    logical, intent(out), optional :: found

    type(json_file) :: json
    integer :: k, pdm
    character(len=slenmax), dimension(:), allocatable :: cval

    !write(*, *) 'l643, vnm:', trim(vnm)
    !write(*, *) 'l643, tnm:', trim(tnm)

    call json%initialize()
    call json%load_file(filename=trim(tnm))
    call json%get('variable_entry.' // trim(vnm) // '.dimensions', cval, found)
    call json%destroy()

    if (.not. found) return
    pdm = size(cval)
    do k = 1, pdm
      if (cval(k)(1:4) == 'time') then
        str = cval(k)
        exit
      end if
    end do

    deallocate(cval)

  end subroutine json_get_timecoord

  ! -----------------------------------------------------------------

  subroutine json_get_vunits(tnm, vnm, units, found)

    use json_module
    implicit none

    character(len=*), intent(in) :: tnm, vnm
    character(len=*), intent(out) :: units
    logical, intent(out), optional :: found

    type(json_file) :: json

    character(len=:), allocatable :: cval

    call json%initialize()
    call json%load_file(filename=trim(tnm))
    call json%get('variable_entry.' // trim(vnm) // '.units', cval, found)
    call json%destroy()
    if (found) then
      units = cval
      deallocate(cval)
    end if

  end subroutine json_get_vunits


  ! -----------------------------------------------------------------

  subroutine json_get_variable()
    real :: r
    integer :: ind

  end subroutine json_get_variable

  ! -----------------------------------------------------------------

  subroutine json_get_cmor_table()
    real :: r
    integer :: ind

  end subroutine json_get_cmor_table

  ! -----------------------------------------------------------------

  subroutine json_get_cmor_coordinate()
    real :: r
    integer :: ind

  end subroutine json_get_cmor_coordinate

  ! -----------------------------------------------------------------

  subroutine json_get_cmor_grid()
    real :: r
    integer :: ind

  end subroutine json_get_cmor_grid

  ! -----------------------------------------------------------------

  subroutine json_get_cmor_cv()
    real :: r
    integer :: ind

  end subroutine json_get_cmor_cv

  ! -----------------------------------------------------------------

end module m_jsons
