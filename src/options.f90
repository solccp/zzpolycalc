module options_module
  use ISO_FORTRAN_ENV 
  implicit none
  save
  integer :: nbuckets
  integer(int64) :: writemark 
  integer :: maxrecords 

contains
  subroutine initialize_options()
    nbuckets = 209715
    writemark =   200000000
    maxrecords =  12000000
  end subroutine
end module
