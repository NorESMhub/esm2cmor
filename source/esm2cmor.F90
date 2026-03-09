program main
   
  use m_namelists, only: verbose, read_namelists, print_namelists, &
    ibasedir, casename, forcefilescan, membertag
  use m_modelsatm, only: atm2cmor
  use m_modelslnd, only: lnd2cmor
  use m_modelsice, only: ice2cmor
  use m_modelsocn, only: ocn2cmor
   
  implicit none
   
  logical :: fileexists
#ifdef MPI
  include 'mpif.h'
  integer :: mpierror, mpirank
   
  ! --- Initialise mpi
  call MPI_INIT(mpierror)
#endif
   
  ! --- Read namelists
  call read_namelists
   
  ! --- Create file list if it does not exist
#ifdef MPI
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpierror)
  if (mpirank .eq. 0) then
#endif
    if (verbose) call print_namelists
    inquire(file='filelist_'//trim(casename)//trim(membertag), &
      exist=fileexists)
    if (.not. fileexists .or. forcefilescan) then
      write(*,*) 'get_file_info: create new file list ' &
        //trim('filelist_'//casename)//trim(membertag)
      if (len_trim(membertag) .gt. 0) then
        call SYSTEM('find '//trim(ibasedir)//'/'//trim(casename) &
          //' -path "*/hist/*" -name "*.*_'//trim(membertag) &
          //'.h*.nc" | sort > '//trim('filelist_'//casename) &
          //trim(membertag))
      else
        call SYSTEM('find '//trim(ibasedir)//'/'//trim(casename) &
          //'/{atm,ice,lnd,ocn,rof}' &
          //' \( -path "*/hist/*" -or ' &
          //'    -path "*/hist_true/*" \)' &
          //' -name "*.nc"' &
          //' | sort > '//trim('filelist_'//casename))
      end if
    else
      write(*,*) 'get_file_info: read existing file list ' &
        //trim('filelist_'//casename)//trim(membertag)
    end if
#ifdef MPI
  end if
  call MPI_BARRIER(MPI_COMM_WORLD, mpierror)
#endif
   
  ! --- Run cmor processing for individual components
  ! call atm2cmor
  ! call lnd2cmor
  call ice2cmor
  call ocn2cmor
  ! call glc2cmor
   
#ifdef MPI
  ! --- Finalise mpi
  call MPI_FINALIZE(mpierror)
   
#endif
  write(*,*)
  write(*,*) '===================='
  write(*,*) '   ALL JOBS DONE'
  write(*,*) '===================='
  write(*,*)
end program main
