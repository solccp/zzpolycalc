!################################ module types_module ###############################
!####################################################################################
module types_module
use ISO_FORTRAN_ENV

!  integer, parameter :: kint = kind(0)
  integer, parameter :: kint = 4
  integer, parameter :: kreal = kind(0.0d0)
  integer, parameter :: maxatoms = 1000
  integer, parameter :: vlongmax = 5
  integer, parameter :: vbase = 1000000000

  type,public :: vlonginteger
    integer,public :: leadpow
    integer(kind=4),private :: tabl(vlongmax)
  end type vlonginteger

  type,public :: structure
    integer(kint) :: nat
    integer(kint),allocatable,dimension(:) :: initiallabel
    integer(kint),allocatable,dimension(:) :: neighbornumber
    integer(kint),allocatable,dimension(:,:) :: neighborlist
    integer(kint) :: order
    type(vlonginteger),allocatable,dimension(:) :: polynomial
    integer(kint), allocatable,dimension(:,:) :: bondlist
    integer(kint) :: nbondlistentries
  end type structure

contains
  subroutine cpvli(a,b)
  type(vlonginteger),intent(in) :: a
  type(vlonginteger),intent(out) :: b
  b=a
  end subroutine cpvli
  


  function setvli(a) result(c)   
  implicit none
  integer(kint),intent(in) :: a
  integer(kint) :: b,i
  type(vlonginteger) :: c

  c%tabl=0
  b=a
  do i=1,vlongmax
    c%tabl(i)=mod(b,vbase)
    b=(b-c%tabl(i))/vbase
    if (b == 0) then
      c%leadpow=i
      exit
    end if
  end do

  end function setvli


  function addvli(a,b) result(c)
  implicit none
  type(vlonginteger),intent(in) :: a,b
  type(vlonginteger) :: c
  integer(kint) :: i,val

  c%tabl=0
  c%leadpow=max(a%leadpow,b%leadpow)
  do i=1,c%leadpow
    c%tabl(i)=c%tabl(i)+a%tabl(i)+b%tabl(i)
  end do
  do i=1,min(vlongmax-1,c%leadpow)
    if (c%tabl(i) >= vbase) then
      val=mod(c%tabl(i),vbase)
      c%tabl(i+1)=c%tabl(i+1)+(c%tabl(i)-val)/vbase
      c%tabl(i)=val
    end if
  end do

  if (c%tabl(vlongmax) >= vbase) then
    print*,"overflow in addvli, enlarge size of tabl"
    stop
  end if
 
  if (c%leadpow < vlongmax) then  
    if (c%tabl(c%leadpow+1) /= 0) c%leadpow=c%leadpow+1
  end if

  end function addvli
 
  function multvli(a,b) result(c)
  implicit none
  type(vlonginteger),intent(in) :: a,b
  type(vlonginteger) :: c
  integer(int64) :: tmp(vlongmax)
  integer(kint) :: i,j,val

  tmp=0
  do i=1,a%leadpow
    do j=1,min(vlongmax-i,b%leadpow)
      tmp(i+j-1)=tmp(i+j-1)+a%tabl(i)*b%tabl(j)
    end do
  end do
  do i=1,min(vlongmax-1,a%leadpow+b%leadpow)
    if (tmp(i) >= vbase) then
      val=mod(tmp(i),vbase)
      tmp(i+1)=tmp(i+1)+(tmp(i)-val)/vbase
      tmp(i)=val
    end if
  end do

  if (tmp(vlongmax) >= vbase) then
    print*,"overflow in multvli, enlarge size of tabl"
    stop
  end if

  c%tabl=tmp

  c%leadpow=0
  do i=min(vlongmax,a%leadpow+b%leadpow+1),1,-1
    if (c%tabl(i) /=0) then
      c%leadpow=i
      exit
    end if
  end do
 
  end function multvli
 
  subroutine printvli(a)
  implicit none
  type(vlonginteger),intent(in) :: a
  integer(kint) :: i,j,val

  if (a%leadpow == 0) then
    write(*,*)"0"
  else
    write(*,'(1X,I4,40I4.4)')(a%tabl(i),i=a%leadpow,1,-1)
  end if

  end subroutine printvli

  subroutine printvlinoadv(a)
  implicit none
  type(vlonginteger),intent(in) :: a
  integer(kint) :: i,j,val

  if (a%leadpow == 0) then
    write(*,*)"0"
  else
    write(*,'(1X,I4,40I4.4)',advance='no')(a%tabl(i),i=a%leadpow,1,-1)
  end if

  end subroutine printvlinoadv



  subroutine print_vli_in_string(pos,string,val)
!   prints integer val in the string at position pos
    implicit none
    integer(kint) :: i
    integer(kint), intent(inout) :: pos
    type(vlonginteger) :: val
    character(len=500) :: string

    do i=val%leadpow,val%leadpow
      select case (val%tabl(i))
        case (0:9)
          write(string(pos:),'(i1)')val%tabl(i)
          pos=pos+1
        case (10:99)
          write(string(pos:),'(i2)')val%tabl(i)
          pos=pos+2
        case (100:999)
          write(string(pos:),'(i3)')val%tabl(i)
          pos=pos+3
        case (1000:9999)
          write(string(pos:),'(i4)')val%tabl(i)
          pos=pos+4
        case (10000:99999)
          write(string(pos:),'(i5)')val%tabl(i)
          pos=pos+5
        case (100000:999999)
          write(string(pos:),'(i6)')val%tabl(i)
          pos=pos+6
        case (1000000:9999999)
          write(string(pos:),'(i7)')val%tabl(i)
          pos=pos+7
        case (10000000:99999999)
          write(string(pos:),'(i8)')val%tabl(i)
          pos=pos+8
        case default
          write(string(pos:),'(i9)')val%tabl(i)
          pos=pos+9
      end select
    end do
    do i=val%leadpow-1,1,-1
      write(string(pos:),'(i9.9)')val%tabl(i)
      pos=pos+9
    end do
    return

  end subroutine print_vli_in_string

end module types_module

!####################################################################################
!############################ end of module types_module ############################

module lookup_module_list
use types_module
use ISO_FORTRAN_ENV
  integer, parameter :: maxat = 1000
  integer, parameter :: packlen = 1024

  integer :: nstruct = 0
  type,public :: neigh
     integer(kint) :: nat
     integer(int64), pointer :: nlist(:)     
     integer(kint) :: order
     type(vlonginteger),pointer :: polynomial(:)
  end type neigh

  type ptrneigh
  type(neigh), pointer :: p
  type(ptrneigh), pointer :: next => null()
  end type ptrneigh

  type(ptrneigh) :: x(maxat)
contains 
subroutine add_neigh(nat,a,order,poly)
  implicit none
  integer(kint),intent(in) :: nat
  integer(kint), intent(in) :: a(3,nat)
  integer(kint), intent(in) :: order
  type(vlonginteger), intent(in) :: poly(order)
  type(vlonginteger) :: polytmp
  type(ptrneigh),pointer :: curr
  type(ptrneigh), target :: xtmp
  logical :: first

  integer :: i,j
  nstruct=nstruct+1
  if (nat>packlen) then 
  print*,"overflow in packlen"
    stop
  end if

  if (nat>maxat) then 
  print*,"overflow in maxat"
    stop
  end if

! find last
  xtmp=x(nat)
  curr=>xtmp
  if (associated(curr%p)) then 
     first=.false.
  else
     first=.true.
  end if
  j=0
  do while(associated(curr%next))
    curr=>curr%next
    j=j+1
  end do
!  write(*,*)first,j
  if (.not. first) then
    allocate(curr%next)
    curr=>curr%next
  end if
  allocate(curr%p)
  allocate(curr%p%nlist(nat))
  do i=1,nat
      curr%p%nlist(i)=a(1,i)+a(2,i)*packlen+a(3,i)*packlen**2
  end do
  allocate(curr%p%polynomial(order+1))
  curr%p%order=order
  do i=1,order+1
    call cpvli(poly(i),curr%p%polynomial(i))
  end do
  if (first) then
    x(nat)%p=>curr%p
  end if
  if (.not. first .and. j.eq.0) then ! second in chain
    x(nat)%next=>curr
  end if
!  write(*,*)nat,curr%p%order,%loc(curr),%loc(x(nat)%next),%loc(x(nat)%p)
!  write(*,*)nstruct
!   write(*,'(A)',advance='no')"X "   
!   do i=1,nat
!     write(*,'(3I3)', advance='no')(a(j,i),j=1,3)
!   end do
!   write(*,*)



end subroutine add_neigh

function check_seen(nat,a,order,poly) result(seen)
  implicit none
  integer(kint),intent(in) :: nat
  integer(kint), intent(in) :: a(3,nat)
  integer(kint), intent(out) :: order
  type(vlonginteger), allocatable,intent(out) :: poly(:)  
  logical :: seen
  integer :: i,j,k
  type(ptrneigh),pointer :: curr
  type(neigh), pointer :: match
  type(ptrneigh), target :: xtmp


  seen = .false.
!  write(*,*)nat
  xtmp=x(nat)
  curr=>xtmp
  structloop:   do while(associated(curr))
   if (.not. associated(curr%p) ) exit structloop
    do j=1,nat
!        write(*,*)j,curr%p%nlist(j),a(1,j)+a(2,j)*packlen+a(3,j)*packlen**2,a(1,j),a(2,j),a(3,j)
        if (curr%p%nlist(j).ne. a(1,j)+a(2,j)*packlen+a(3,j)*packlen**2) then 
          curr=>curr%next
          cycle structloop
        end if
      if (j.eq.nat) then 
        seen=.true.
!        write(*,*)'match'
        match=>curr%p
        exit structloop
      end if
    end do
    curr=>curr%next
 end do structloop
 if (seen) then
    order=match%order
    allocate(poly(order+1))
     do i=1,order+1
      call cpvli(match%polynomial(i),poly(i))
    end do
 end if
end function check_seen


end module lookup_module_list




module lookup_module_crc
use types_module
use ISO_FORTRAN_ENV
  integer, parameter :: maxtab = 2097152
  integer(int32), parameter :: packshift = 10
  integer(int32), parameter :: packlen = 2**packshift

  integer :: nstruct = 0
  type,public :: neigh
     integer(int32), pointer :: nlist(:)     
     integer(kint) :: order
     type(vlonginteger),pointer :: polynomial(:)
  end type neigh

  type ptrneigh
  type(neigh), pointer :: p(:)
  end type ptrneigh

  type(ptrneigh) :: x(maxtab)
  integer(kint) :: xlen(maxtab)
  integer(kint) :: xuse(maxtab,100)
  interface
  function crc32_hash(a,cont) result(crc64)
    use,intrinsic :: ISO_FORTRAN_ENV, only : int32,int64
    integer(int64)               :: crc64
    logical,intent(in),optional  :: cont
    character(len=1),intent(in)  :: a(:)
  end function crc32_hash
  end interface


contains 
subroutine add_neigh(nat,a,order,poly)
  implicit none
  integer(kint),intent(in) :: nat
  integer(kint), intent(in) :: a(3,nat)
  integer(kint), intent(in) :: order
  type(vlonginteger), intent(in) :: poly(order+1)
  type(vlonginteger) :: polytmp
  type(ptrneigh),pointer :: curr
  type(ptrneigh), target :: xtmp
  logical :: first
  type(neigh), pointer :: ptmp(:)
  integer(int64) :: ipack,idx2l,idx1l
  character(len=1) :: buf (3*6*kint),buf2
  integer :: i,j,ilen,idx2,idx1


!  if (nat.le.10 ) return

  nstruct=nstruct+1
  if (mod(nstruct,10000).eq.0) write(*,*)'nstruct',nstruct ! print progress
  if (nat>packlen) then 
  print*,"overflow in packlen"
    stop
  end if

!  if (nat>maxat) then 
!  print*,"overflow in maxat"
!    stop
!  end if


  buf=transfer(a(1:3,1:6),buf)
!first neighbor is rather wastful. Add natoms
!  buf(1:kint)=transfer(nat,buf2)
  idx1l=crc32_hash(buf)
  if (nat.ge.12) then
    buf=transfer(a(1:3,nat/2:nat/2+6),buf)
    idx1l=idx1l+crc32_hash(buf)
  end if

 if (nat.ge.18) then
    buf=transfer(a(1:3,2*nat/3:2*nat/3+6),buf)
    idx1l=idx1l+crc32_hash(buf)
  end if

  idx1=mod(idx1l,maxtab)
  idx1=iabs(idx1)+1


  xtmp=x(idx1)

  

  curr=>xtmp
  if (associated(curr%p)) then 
     first=.false.
  else
     first=.true.
     xlen(idx1)=0
  end if
  
  ilen=xlen(idx1)
!  write(*,*)nat,ilen
  allocate(ptmp(ilen+1))

!  xuse(idx1,ilen)=0
   
  do i=1,ilen
!      allocate(ptmp(i)%nlist(nat))
!      write(*,*)'alloc',i,%loc(ptmp(i)%nlist)
!      allocate(ptmp(i)%polynomial(curr%p(i)%order+1))
!I have no idea why but this assign pointers and not values
      ptmp(i)%order=x(idx1)%p(i)%order
      ptmp(i)%nlist=>x(idx1)%p(i)%nlist
      ptmp(i)%polynomial=>x(idx1)%p(i)%polynomial
  end do

  allocate(ptmp(ilen+1)%nlist(nat))
!  write(*,*)'alloc ilen+1',ilen+1,%loc(ptmp(ilen+1)%nlist)


  do i=1,nat
      ptmp(ilen+1)%nlist(i)=a(1,i)+ishft(a(2,i),packshift)+ishft(a(3,i),2*packshift)
  end do
  allocate(ptmp(ilen+1)%polynomial(order+1))
  ptmp(ilen+1)%order=order
  do i=1,order+1
    call cpvli(poly(i),ptmp(ilen+1)%polynomial(i))
  end do

!  do i=1,ilen
!    write(*,*)'Dealloc nlist',i,%loc(x(nat)%p(i)%nlist),%loc(x(nat)%p(i)),%loc(ptmp(i)%nlist),%loc(ptmp(i))
!    deallocate(x(nat)%p(i)%nlist)
!    write(*,*)'Dealloc polynomial',i
!    deallocate(x(nat)%p(i)%polynomial)
!  end do

  if (.not.first)  then 
!     write(*,*)'Dealloc x(nat)%p ptmp',%loc(x(nat)%p),%loc(ptmp)
   deallocate(x(idx1)%p)
  end if
  
!  write(*,*)'ptmp',%loc(ptmp)
  x(idx1)%p=>ptmp
  xlen(idx1)=xlen(idx1)+1

!   write(*,'(A,I5)',advance='no')"P ",nat
!   do i=1,order+1
!     call printvlinoadv(poly(i))
!   end do
!   write(*,*)



!  write(*,*)nat,curr%p%order,%loc(curr),%loc(x(nat)%next),%loc(x(nat)%p)
!  write(*,*)nstruct
!   write(*,'(A)',advance='no')"X "   
!   write(*,'(2I3)'),nat,iabs(a(1,nat)-a(1,nat/2))+iabs(a(3,nat/2)-a(3,1))
!   do i=1,nat
!     write(*,'(3I3)', advance='no')(a(j,i),j=1,3)
!   end do
!   write(*,*)

end subroutine add_neigh

function check_seen(nat,a,order,poly) result(seen)
  implicit none
  integer(kint),intent(in) :: nat
  integer(kint), intent(in) :: a(3,nat)
  integer(kint), intent(out) :: order
  type(vlonginteger), allocatable,intent(out) :: poly(:)  
  logical :: seen
  integer :: i,j,k,idx2,idx1
  type(ptrneigh),pointer :: curr
  type(neigh), pointer :: match
  type(ptrneigh), target :: xtmp
  integer(int64) :: idx2l,idx1l
  character(len=1) :: buf (3*6*kint),buf2

!  if (nat.le.10 ) return
  buf=transfer(a(1:3,1:6),buf)
!first neighbor is rather wastful. Add natoms
!  buf(1:kint)=transfer(nat,buf2)



!  buf=transfer(a,buf)
!  write(*,*)'size buf',size(buf)
  idx1l=crc32_hash(buf)
  if (nat.ge.12) then
    buf=transfer(a(1:3,nat/2:nat/2+6),buf)
    idx1l=idx1l+crc32_hash(buf)
  end if
  if (nat.ge.18) then
    buf=transfer(a(1:3,2*nat/3:2*nat/3+6),buf)
    idx1l=idx1l+crc32_hash(buf)
  end if


  idx1=mod(idx1l,maxtab)
  idx1=iabs(idx1)+1


  seen = .false.
!  write(*,*)nat,xlen(nat)
  xtmp=x(idx1)
  curr=>xtmp
  structloop:   do i=1,xlen(idx1)
    do j=1,nat
!        write(*,*)i,j,curr%p(i)%nlist(j),a(1,j)+a(2,j)*packlen+a(3,j)*packlen**2,a(1,j),a(2,j),a(3,j)
        if (curr%p(i)%nlist(j).ne. a(1,j)+ishft(a(2,j),packshift)+ishft(a(3,j),2*packshift)) then 
          cycle structloop
        end if
      if (j.eq.nat) then 
        seen=.true.
!        write(*,*)'match'
        match=>curr%p(i)
        exit structloop
      end if
    end do
 end do structloop
 if (seen) then 
!    xuse(idx1,i)=xuse(idx1,i)+1
    order=match%order
    allocate(poly(order+1))
     do i=1,order+1
      call cpvli(match%polynomial(i),poly(i))
    end do
 end if
end function check_seen


end module lookup_module_crc


module lookup_module_md5
use types_module
use ISO_FORTRAN_ENV
use, intrinsic :: iso_c_binding
  integer, parameter :: maxtab = 2097152
!  integer, parameter :: highmark =  50000000
integer, parameter :: highmark =  5000000
!  integer, parameter :: maxtab = 100
  integer(int32), parameter :: packshift = 10
  integer(int32), parameter :: packlen = 2**packshift
  logical :: firstrun = .true.
  integer :: nstruct = 0
  integer :: nstructall = 0
  type,public :: neigh
     integer(C_signed_char), pointer :: nlist(:)     
     integer(kint) :: order
     integer(kint) :: iseen
     type(vlonginteger),pointer :: polynomial(:)
  end type neigh

  type ptrneigh
  type(neigh), pointer :: p(:)
  end type ptrneigh

  type(ptrneigh),allocatable :: x(:)
  integer(kint),allocatable :: xlen(:)
  integer(kint),allocatable :: irepl(:)
  interface
  function crc32_hash(a,cont) result(crc64)
    use,intrinsic :: ISO_FORTRAN_ENV, only : int32,int64
    integer(int64)               :: crc64
    logical,intent(in),optional  :: cont
    character(len=1),intent(in)  :: a(:)
  end function crc32_hash
    subroutine MD5(dat,size,result) bind(C, name='MD5')
    use, intrinsic :: iso_c_binding
    integer(C_signed_char) :: result(*)
    integer(C_long), value :: size
    character(kind=c_char) :: dat(size)
  end subroutine MD5
  end interface


contains 
subroutine add_neigh(nat,a,order,poly)
  implicit none
  integer(kint),intent(in) :: nat
  integer(kint), intent(in) :: a(3,nat)
  integer(int16) :: asmall(3,nat)
  integer(kint), intent(in) :: order
  type(vlonginteger), intent(in) :: poly(order+1)
  type(vlonginteger) :: polytmp
  type(ptrneigh),pointer :: curr
  type(ptrneigh), target :: xtmp
  logical :: first
  type(neigh), pointer :: ptmp(:)
  integer(int64) :: ipack,idx2l,idx1l
  character(len=1) :: buf (3*nat*2),buf2
  integer :: i,j,ilen,idx2,idx1
  integer(C_signed_char) :: md5sum(16)
  logical :: replace
  integer :: mem

!  if (nat.le.10 ) return

  if (nat>packlen) then 
  print*,"overflow in packlen"
    stop
  end if

!  if (nat>maxat) then 
!  print*,"overflow in maxat"
!    stop
!  end if

  nstructall=nstructall+1
  if (mod(nstructall,1000000).eq.0)  then
    call system_mem_usage(mem)
    write(*,*)'nstruct',nstructall,nstruct,'mem',mem
  end if 

  asmall=a
  buf=transfer(asmall(1:3,1:nat),buf)
  call MD5(buf,size(buf,1,C_long),md5sum)
  idx1l=transfer(md5sum,idx1l)

!  idx1l=crc32_hash(buf)
!  if (nat.ge.12) then
!    buf=transfer(a(1:3,nat/2:nat/2+6),buf)
!    idx1l=idx1l+crc32_hash(buf)
!  end if

! if (nat.ge.18) then
!    buf=transfer(a(1:3,2*nat/3:2*nat/3+6),buf)
!    idx1l=idx1l+crc32_hash(buf)
!  end if

  idx1=mod(idx1l,maxtab)
  idx1=iabs(idx1)+1


  xtmp=x(idx1)

  

  curr=>xtmp
  if (associated(curr%p)) then 
     first=.false.
  else
     first=.true.
     xlen(idx1)=0
     irepl(idx1)=1
  end if
  
  ilen=xlen(idx1)
  if (nstruct.gt.highmark .and. .not. first) then
    replace=.true.
  else
    replace=.false.
  end if

  if (.not. replace) then

!  write(*,*)nat,ilen
  allocate(ptmp(ilen))

!   write(*,*)ilen,sizeof(ptmp)
!  xuse(idx1,ilen)=0
   
  do i=1,ilen
      allocate(ptmp(i)%nlist(16))
      allocate(ptmp(i)%polynomial(x(idx1)%p(i)%order+1))
    
      ptmp(i)%order=x(idx1)%p(i)%order
      ptmp(i)%iseen=x(idx1)%p(i)%iseen
      ptmp(i)%nlist=x(idx1)%p(i)%nlist
      ptmp(i)%polynomial=x(idx1)%p(i)%polynomial
      deallocate(x(idx1)%p(i)%nlist,x(idx1)%p(i)%polynomial)
  end do

  if (.not. first) deallocate(x(idx1)%p)
  
  allocate(x(idx1)%p(ilen+1))
  do i=1,ilen
      allocate(x(idx1)%p(i)%nlist(16))
      allocate(x(idx1)%p(i)%polynomial(ptmp(i)%order+1))
    
      x(idx1)%p(i)%order=ptmp(i)%order
      x(idx1)%p(i)%iseen=ptmp(i)%iseen
      x(idx1)%p(i)%nlist=ptmp(i)%nlist
      x(idx1)%p(i)%polynomial=ptmp(i)%polynomial
      deallocate(ptmp(i)%nlist,ptmp(i)%polynomial)
  end do
  deallocate(ptmp)

  allocate(x(idx1)%p(ilen+1)%nlist(16))
  allocate(x(idx1)%p(ilen+1)%polynomial(order+1))


  do i=1,16
      x(idx1)%p(ilen+1)%nlist(i)=md5sum(i)
  end do
  x(idx1)%p(ilen+1)%order=order
  x(idx1)%p(ilen+1)%iseen=0
  do i=1,order+1
    call cpvli(poly(i),x(idx1)%p(ilen+1)%polynomial(i))
  end do

!  do i=1,ilen
!    write(*,*)'Dealloc nlist',i,%loc(x(nat)%p(i)%nlist),%loc(x(nat)%p(i)),%loc(ptmp(i)%nlist),%loc(ptmp(i))
!    deallocate(x(nat)%p(i)%nlist)
!    write(*,*)'Dealloc polynomial',i
!    deallocate(x(nat)%p(i)%polynomial)
!  end do

 ! if (.not.first)  then 
!     write(*,*)'Dealloc x(nat)%p ptmp',%loc(x(nat)%p),%loc(ptmp)
!    deallocate(x(idx1)%p)
!  end if
  
!  write(*,*)'ptmp',%loc(ptmp)
!  x(idx1)%p=>ptmp
  xlen(idx1)=xlen(idx1)+1
  nstruct=nstruct+1
  if (nstruct.eq.highmark) write(*,*)'High mark',nstruct,' achieved' 
else ! replace

  j=irepl(idx1) 
  do i=1,16
      x(idx1)%p(j)%nlist(i)=md5sum(i)
  end do
 
  if (x(idx1)%p(j)%order .ne. order) then
     deallocate(x(idx1)%p(j)%polynomial) 
     allocate(x(idx1)%p(j)%polynomial(order+1))
  end if
 

  x(idx1)%p(j)%order=order
  x(idx1)%p(j)%iseen=0
 
  do i=1,order+1
    call cpvli(poly(i),x(idx1)%p(j)%polynomial(i))
  end do

  j=j+1
  if (j.gt.xlen(idx1)) j=1
  irepl(idx1)=j
end if

!   write(*,'(A,I5)',advance='no')"P ",nat
!   do i=1,order+1
!     call printvlinoadv(poly(i))
!   end do
!   write(*,*)



!  write(*,*)nat,curr%p%order,%loc(curr),%loc(x(nat)%next),%loc(x(nat)%p)
!  write(*,*)nstruct
!   write(*,'(A)',advance='no')"X "   
!   write(*,'(2I3)'),nat,iabs(a(1,nat)-a(1,nat/2))+iabs(a(3,nat/2)-a(3,1))
!   do i=1,nat
!     write(*,'(3I3)', advance='no')(a(j,i),j=1,3)
!   end do
!   write(*,*)

end subroutine add_neigh

function check_seen(nat,a,order,poly) result(seen)
  implicit none
  integer(kint),intent(in) :: nat
  integer(kint), intent(in) :: a(3,nat)
  integer(int16) :: asmall(3,nat)
  integer(kint), intent(out) :: order
  type(vlonginteger), allocatable,intent(out) :: poly(:)  
  logical :: seen
  integer :: i,j,k,idx2,idx1
  type(ptrneigh),pointer :: curr
  type(neigh), pointer :: match
  type(ptrneigh), target :: xtmp
  integer(int64) :: idx2l,idx1l
  character(len=1) :: buf (3*nat*2),buf2
  integer(C_signed_char) :: md5sum(16)

  type(neigh) :: temp
  integer(int64), pointer :: temp2


  if (firstrun) then
    allocate(x(maxtab),xlen(maxtab),irepl(maxtab))
    write(*,*)sizeof(x),sizeof(xlen),sizeof(curr),sizeof(temp),sizeof(temp2)
    firstrun=.false.
  end if

!  if (nat.le.10 ) return
  asmall=a
  buf=transfer(asmall(1:3,1:nat),buf)

  call MD5(buf,size(buf,1,C_long),md5sum)
  idx1l=transfer(md5sum,idx1l)



!  idx1l=crc32_hash(buf)
!  if (nat.ge.12) then
!    buf=transfer(a(1:3,nat/2:nat/2+6),buf)
!    idx1l=idx1l+crc32_hash(buf)
!  end if
!  if (nat.ge.18) then
!    buf=transfer(a(1:3,2*nat/3:2*nat/3+6),buf)
!    idx1l=idx1l+crc32_hash(buf)
!  end if


  idx1=mod(idx1l,maxtab)
  idx1=iabs(idx1)+1


  seen = .false.
!  write(*,*)nat,xlen(nat)
  xtmp=x(idx1)
  curr=>xtmp
  structloop:   do i=1,xlen(idx1)
    do j=1,16
        if (curr%p(i)%nlist(j).ne. md5sum(j)) then 
          cycle structloop
        end if
      if (j.eq.16) then 
        seen=.true.
!        write(*,*)'match'
        match=>curr%p(i)
        exit structloop
      end if
    end do
 end do structloop
 if (seen) then 
!    xuse(idx1,i)=xuse(idx1,i)+1
    order=match%order
    x(idx1)%p(i)%iseen=x(idx1)%p(i)%iseen+1
    allocate(poly(order+1))
     do i=1,order+1
      call cpvli(match%polynomial(i),poly(i))
    end do
 end if
end function check_seen


!from stackexchange
subroutine system_mem_usage(valueRSS)
use ifport !if on intel compiler
implicit none
integer, intent(out) :: valueRSS

character(len=200):: filename=' '
character(len=80) :: line
character(len=8)  :: pid_char=' '
integer :: pid
logical :: ifxst

valueRSS=-1    ! return negative number if not found

!--- get process ID

pid=getpid()
write(pid_char,'(I8)') pid
filename='/proc/'//trim(adjustl(pid_char))//'/status'

!--- read system file

inquire (file=filename,exist=ifxst)
if (.not.ifxst) then
  write (*,*) 'system file does not exist'
  return
endif

open(unit=100, file=filename, action='read')
do
  read (100,'(a)',end=120) line
  if (line(1:6).eq.'VmRSS:') then
     read (line(7:),*) valueRSS
     exit
  endif
enddo
120 continue
close(100)

return
end subroutine system_mem_usage


end module lookup_module_md5

