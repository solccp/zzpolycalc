module options_module
  use ISO_FORTRAN_ENV 
  implicit none
  save
  integer :: nbuckets
  integer(int64) :: writemark 
  integer :: maxrecords 
  logical :: is_adjacencyfile
  logical :: has_read_cache_file
  character(len=100) :: read_cache_fname
  logical :: has_write_cache_file
  character(len=100) :: write_cache_fname
  logical :: print_bondlevel 
contains
  subroutine initialize_options()
    nbuckets = 2097152
    writemark =   200000000
    maxrecords =  52000000
    is_adjacencyfile = .false.
    print_bondlevel = .false. 
    has_read_cache_file = .false.
    has_write_cache_file = .false. 
  end subroutine
end module
