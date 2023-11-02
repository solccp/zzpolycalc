module options_module
  use ISO_FORTRAN_ENV 
  implicit none
  save
  integer :: nbuckets
  integer(int64) :: writemark 
  integer :: maxrecords 
  logical :: is_adjacencyfile
  logical :: print_bondlevel 
contains
  subroutine initialize_options()
    nbuckets = 2097152
    writemark =   200000000
    maxrecords =  12000000
    is_adjacencyfile = .false.
    print_bondlevel = .false. 
  end subroutine
end module
