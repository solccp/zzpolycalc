module options_module
  use ISO_FORTRAN_ENV 
  implicit none
  save
  integer :: nbuckets
  integer :: chunksize
  integer(int64) :: writemark 
  integer :: maxrecords 
  logical :: is_adjacencyfile
  logical :: has_read_cache_file
  character(len=100) :: read_cache_fname
  logical :: has_write_cache_file
  character(len=100) :: write_cache_fname
  logical :: print_bondlevel 
  logical :: print_XML
  logical :: verbose
contains
  subroutine initialize_options()
    chunksize = 1
    nbuckets = 2097152
    writemark =   200000000
    maxrecords =  52000000
    is_adjacencyfile = .false.
    print_bondlevel = .false.
    print_XML = .false. 
    verbose = .false.
    has_read_cache_file = .false.
    has_write_cache_file = .false. 
  end subroutine
  subroutine print_options()
    use types_module
!    write(*,'(a,i)')'chunksize',chunksize
    write(*,'(a,i)')'nbuckets',nbuckets
    write(*,'(a,i)')'maxrecords',maxrecords
    if (has_read_cache_file) write(*,'(a,a)')'Read cache from:  ',trim(read_cache_fname)
    if (has_write_cache_file) write(*,'(a,a)')'Write cache to:  ',trim(write_cache_fname)
    write(*,'(a,i)')'vlongmax',vlongmax
  end subroutine print_options   
    
    
end module
