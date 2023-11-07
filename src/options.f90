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
  logical :: unsorted_geometry
  logical :: read_connection_table
contains
  subroutine initialize_options()
    chunksize = 1
    nbuckets = 1000000
    writemark =   200000000
    maxrecords =  20000000
    is_adjacencyfile = .false.
    print_bondlevel = .false.
    print_XML = .false. 
    verbose = .false.
    unsorted_geometry = .false. 
    has_read_cache_file = .false.
    has_write_cache_file = .false. 
    read_connection_table = .false. 
  end subroutine
  subroutine print_options()
    use types_module
    write(*,'(a,X,i0)')'nbuckets',nbuckets
    write(*,'(a,X,i0)')'maxrecords',maxrecords
    if (has_read_cache_file) write(*,'(a,a)')'Read cache from:  ',trim(read_cache_fname)
    if (has_write_cache_file) write(*,'(a,a)')'Write cache to:  ',trim(write_cache_fname)
    write(*,'(a,X,i0)')'vlongmax',vlongmax
    if (is_adjacencyfile) write(*,*)'Input file is adjacency matrix'
    if (unsorted_geometry) write(*,*)'Unmodified geometry used' 
    if (read_connection_table) write(*,*)'Connection table read from the bottom of the xyz file'
  end subroutine print_options   
    
    
end module
