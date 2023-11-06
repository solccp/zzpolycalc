! Simple insertion sort of an XYZ array
subroutine sort(a,n)
  use types_module
  implicit none
  
  integer, intent(in) :: n
  real(kreal), intent(inout) :: a(3,n)
  real(kreal) :: tmp(3)
  integer :: i,j

  do i=2,n
    tmp=a(:,i)
    j=i-1
    do while (j>=1) ! equivalent to sort -k 3,3g -k 2,2g -k 1,1g
      if (a(3,j) < tmp(3) .or.(a(3,j) == tmp(3) .and. a(2,j) < tmp(2)) .or. &
             ( a(3,j) == tmp(3) .and. a(2,j) == tmp(2) .and. a(1,j) < tmp(1))) exit 
      a(:,j+1)=a(:,j)
      j=j-1
    end do
    a(:,j+1)=tmp
  end do
end subroutine
