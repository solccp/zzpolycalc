!################################ module hash_module ###############################
!####################################################################################

module hash_module
#ifdef USE_SHA256
  integer, parameter :: hashsize = 32
#else
  integer, parameter :: hashsize = 16
#endif
  interface 
#ifdef USE_XXHASH
  subroutine hash(dat,size,result) bind(C, name= 'xxhash')
#else
#ifdef USE_SHA256
  subroutine hash(dat,size,result) bind(C, name= 'SHA256')
#else
  subroutine hash(dat,size,result) bind(C, name= 'MD5')
#endif
#endif
    use, intrinsic :: iso_c_binding
    integer(C_signed_char) :: result(*)
    integer(C_long), value :: size
    character(kind=c_char) :: dat(size)
  end subroutine hash
  end interface
end module hash_module

!################################ module lookup_module_hash #########################
!####################################################################################


module lookup_module_hash
use types_module
use hash_module
use ISO_FORTRAN_ENV
use, intrinsic :: iso_c_binding
!  integer, parameter :: maxtab = 2097152
  integer, parameter :: maxtab = 209715
!  integer, parameter :: highmark =   180 000 000
  integer(int64), parameter :: writemark =   200000000
!integer(int64), parameter :: writemark =   1000000
  integer, parameter :: highmark =  12000000
!  integer, parameter :: maxtab = 100
  integer(int32), parameter :: packshift = 10
  integer(int32), parameter :: packlen = 2**packshift
  logical :: firstrun = .true.
  integer(int64) :: nstruct = 0
  integer(int64) :: nstructall = 0
  type,public :: neigh
     integer(C_signed_char), pointer :: nlist(:)     
     integer(int16) :: order
     integer(int16) :: nat
     integer(kint) :: iseen
     integer(int64) :: lastseen
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
  end interface


contains 
subroutine add_neigh(nat,nbnum,a,order,poly,hashseen,iseen,lastseen,duringread)
!  USE IFPORT ! for rename
  use hash_module
  implicit none
  integer(kint),intent(in) :: nat
  integer(kint), intent(in) :: nbnum(nat)
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
  real(8)        :: highscore,score
  character(len=1) :: buf (3*nat*2),buf2
  integer :: i,j,ilen,idx2,idx1,ihighscore
  integer(C_signed_char) :: hashsum(hashsize)
  integer(C_signed_char),OPTIONAL :: hashseen(hashsize)
  integer(kint),OPTIONAL :: iseen
  integer(int64),OPTIONAL :: lastseen
  logical,OPTIONAL :: duringread

  logical :: replace
  integer :: mem,ires

!  if (nat.le.10 ) return

!  if (nat>packlen) then 
!  print*,"overflow in packlen"
!    stop
!  end if

!  if (nat>maxat) then 
!  print*,"overflow in maxat"
!    stop
!  end if

  if (.not. present(duringread)) then 
    nstructall=nstructall+1
    if (mod(nstructall,1000000).eq.0)  then
      call system_mem_usage(mem)
      write(*,*)'nstruct',nstructall,nstruct,'mem',mem
    end if 
  end if


  
  asmall=0 ! fills outside of neighbornumber with zeroes

  if(.not. present(hashseen)) then
    do i=1,nat
     do j=1,nbnum(i)
       asmall(j,i)=a(j,i)
     end do
    end do
    buf=transfer(asmall(1:3,1:nat),buf)
    call hash(buf,size(buf,1,C_long),hashsum)
  else
    hashsum=hashseen
  end if
  idx1l=transfer(hashsum,idx1l)

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
      allocate(ptmp(i)%nlist(hashsize))
      allocate(ptmp(i)%polynomial(x(idx1)%p(i)%order+1))
    
      ptmp(i)%order=x(idx1)%p(i)%order
      ptmp(i)%iseen=x(idx1)%p(i)%iseen
      ptmp(i)%nat=x(idx1)%p(i)%nat
      ptmp(i)%lastseen=x(idx1)%p(i)%lastseen
      ptmp(i)%nlist=x(idx1)%p(i)%nlist
      ptmp(i)%polynomial=x(idx1)%p(i)%polynomial
      deallocate(x(idx1)%p(i)%nlist,x(idx1)%p(i)%polynomial)
  end do

  if (.not. first) deallocate(x(idx1)%p)
  
  allocate(x(idx1)%p(ilen+1))
  do i=1,ilen
      allocate(x(idx1)%p(i)%nlist(hashsize))
      allocate(x(idx1)%p(i)%polynomial(ptmp(i)%order+1))
    
      x(idx1)%p(i)%order=ptmp(i)%order
      x(idx1)%p(i)%iseen=ptmp(i)%iseen
      x(idx1)%p(i)%nat=ptmp(i)%nat
      x(idx1)%p(i)%lastseen=ptmp(i)%lastseen
      x(idx1)%p(i)%nlist=ptmp(i)%nlist
      x(idx1)%p(i)%polynomial=ptmp(i)%polynomial
      deallocate(ptmp(i)%nlist,ptmp(i)%polynomial)
  end do
  deallocate(ptmp)

  allocate(x(idx1)%p(ilen+1)%nlist(hashsize))
  allocate(x(idx1)%p(ilen+1)%polynomial(order+1))


  do i=1,hashsize
      x(idx1)%p(ilen+1)%nlist(i)=hashsum(i)
  end do
  x(idx1)%p(ilen+1)%order=order
  if (.not.present(iseen)) then 
    x(idx1)%p(ilen+1)%iseen=0
  else
    x(idx1)%p(ilen+1)%iseen=iseen
  end if
  x(idx1)%p(ilen+1)%nat=nat
  if (.not.present(lastseen)) then
    x(idx1)%p(ilen+1)%lastseen=nstructall
  else
    x(idx1)%p(ilen+1)%lastseen=lastseen
  end if
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

! find best candidate to replace

highscore=0
ihighscore=0
do j=1,xlen(idx1)
  score=nstructall-x(idx1)%p(j)%lastseen
  score=score/(x(idx1)%p(j)%iseen+1)
  score=score/(x(idx1)%p(j)%nat**2) 
!  score=score/(x(idx1)%p(j)%nat) 
!  write(*,*)'score',j,score,x(idx1)%p(j)%lastseen,x(idx1)%p(j)%iseen
  if (score.gt.highscore .or. j.eq.1) then
     highscore=score
     ihighscore=j
  end if
end do
  j=ihighscore
!  j=irepl(idx1)
!  write(*,*)'highscore',j,highscore,xlen(idx1)
  do i=1,hashsize
      x(idx1)%p(j)%nlist(i)=hashsum(i)
  end do
 
  if (x(idx1)%p(j)%order .ne. order) then
     deallocate(x(idx1)%p(j)%polynomial) 
     allocate(x(idx1)%p(j)%polynomial(order+1))
  end if
 

  x(idx1)%p(j)%order=order
  x(idx1)%p(j)%iseen=0
  x(idx1)%p(j)%nat=nat
  x(idx1)%p(j)%lastseen=nstructall
 
  do i=1,order+1
    call cpvli(poly(i),x(idx1)%p(j)%polynomial(i))
  end do

!  j=j+1
!  if (j.gt.xlen(idx1)) j=1
  irepl(idx1)=j
end if

  if (.not.present(hashseen) .and. .not. present(duringread)) then
   if (mod(nstructall,writemark).eq.0)  then
    write(*,*)'Saving cache'
    call execute_command_line ("mv cache.bin cache.bin.bak")
!     ires=rename('cache.bin','cache.bin.bak')
     call writetodisk
!    call execute_command_line ("rm cache.bin.bak")
     open(unit=23, iostat=ires, file='cache.bin.bak', status='old')
     if (ires .eq. 0) close(23, status='delete')
   end if
  end if


!   write(*,'(A,I5)',advance='no')"P ",nat
!   do i=1,order+1
!     call printvlinoadv(poly(i))
!   end do
!   write(*,*)



!  write(*,*)nat,curr%p%order,%loc(curr),%loc(x(nat)%next),%loc(x(nat)%p)
!  write(*,*)nstruct
!   write(*,'(A,I7)',advance='no')"X ",nat
!   do i=1,nat
!     write(*,'(3I5)', advance='no')(asmall(j,i),j=1,3)
!   end do
!   write(*,*)

end subroutine add_neigh


subroutine writetodisk
use hash_module
implicit none
  type(ptrneigh),pointer :: curr
  type(ptrneigh), target :: xtmp
  logical :: first
  integer(int64) :: ipack,idx2l,idx1l
  real(8)        :: highscore,score
  integer :: i,j,ilen,idx2,idx1,ihighscore
  integer(C_signed_char) :: hashsum(hashsize)
  integer(kint) :: leadpowmax=0
!  call MD5(buf,size(buf,1,C_long),hashsum)
!  idx1l=transfer(hashsum,idx1l)

!  idx1l=crc32_hash(buf)
!  if (nat.ge.12) then
!    buf=transfer(a(1:3,nat/2:nat/2+6),buf)
!    idx1l=idx1l+crc32_hash(buf)
!  end if

! if (nat.ge.18) then
!    buf=transfer(a(1:3,2*nat/3:2*nat/3+6),buf)
!    idx1l=idx1l+crc32_hash(buf)
!  end if

  open(23,file='cache.bin',FORM='UNFORMATTED')
  write(23)vlongmax,nstructall

  do idx1=1,maxtab
    xtmp=x(idx1)
    curr=>xtmp
    if (associated(curr%p)) then 
       first=.false.
       ilen=xlen(idx1)
    else
       first=.true.
       ilen=0
    end if
    do i=1,ilen
      write(23)x(idx1)%p(i)%order,x(idx1)%p(i)%iseen,x(idx1)%p(i)%nat,x(idx1)%p(i)%lastseen,x(idx1)%p(i)%nlist
      call writetofile(23,x(idx1)%p(i)%polynomial,x(idx1)%p(i)%order+1)
      do j=1,x(idx1)%p(i)%order+1
        if (x(idx1)%p(i)%polynomial(j)%leadpow.gt.leadpowmax) leadpowmax=x(idx1)%p(i)%polynomial(j)%leadpow
!        write(*,*)leadpowmax,j,x(idx1)%p(i)%order
      end do
    end do
  end do
  write(*,*)'Max large integer size: ',leadpowmax
  close(23)
end subroutine writetodisk


subroutine readfromdisk
use hash_module
implicit none
  integer(C_signed_char) :: hashsum(hashsize)
  integer :: vlong,ires
  integer(kint) :: nbnum(maxatoms)
  integer(kint) ::  a(3,maxatoms)
     integer(int16) :: order
     integer(int16) :: nat
     integer(kint) :: order32
     integer(kint) :: nat32
     integer(kint) :: iseen
     integer(int64) :: lastseen
     type(vlonginteger),allocatable :: poly(:)

   allocate(x(maxtab),xlen(maxtab),irepl(maxtab))
   firstrun=.false.
   xlen=0
  write(*,*)'Reading cache.bin'
  open(23,file='cache.bin',FORM='UNFORMATTED')
  read(23)vlong,nstructall
  if (vlong.gt.vlongmax) then
    write(*,*)'WARNING!'
    write(*,*)'Size of vlong on file is larger than in the current code'
    write(*,*)'The readfromfile procedure may not catch corruption'
    write(*,*)'Inspect the results carefully'
  end if
  ires=0
  do while (ires.eq.0)
  read(23,IOSTAT=ires)order,iseen,nat,lastseen,hashsum
  if (ires.eq.0) then
   allocate(poly(order+1))
   order32=order
   nat32=nat
!   write(*,*)order,iseen,nat,lastseen,hashsum
   call readfromfile(23,poly,order32+1)
   call add_neigh(nat32,nbnum,a,order32,poly,hashsum,iseen,lastseen,.true.)
   deallocate(poly)
  end if
  end do
  close(23)  
  write(*,*)'cache.bin processed'

end subroutine readfromdisk


function check_seen(nat,nbnum,a,order,poly) result(seen)
  use hash_module
  implicit none
  integer(kint),intent(in) :: nat
  integer(kint), intent(in) :: nbnum(nat)
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
  integer(C_signed_char) :: hashsum(hashsize)

  type(neigh) :: temp
  integer(int64), pointer :: temp2


  if (firstrun) then
    allocate(x(maxtab),xlen(maxtab),irepl(maxtab))
!    write(*,*)sizeof(x),sizeof(xlen),sizeof(curr),sizeof(temp),sizeof(temp2)
    firstrun=.false.
    xlen=0
  end if

!  if (nat.le.10 ) return

  asmall=0 ! fills outside of neighbornumber with zeroes
  do i=1,nat
   do j=1,nbnum(i)
     asmall(j,i)=a(j,i)
   end do
  end do

  buf=transfer(asmall(1:3,1:nat),buf)

  call hash(buf,size(buf,1,C_long),hashsum)
  idx1l=transfer(hashsum,idx1l)



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
    do j=1,hashsize
        if (curr%p(i)%nlist(j).ne. hashsum(j)) then 
          cycle structloop
        end if
      if (j.eq.hashsize) then 
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
    x(idx1)%p(i)%lastseen=nstructall
    allocate(poly(order+1))
     do i=1,order+1
      call cpvli(match%polynomial(i),poly(i))
    end do
 end if
!   write(*,'(A,I5,L)',advance='no')"F ",nat,seen   
!   do i=1,nat
!     write(*,'(3I4)', advance='no')(a(j,i),j=1,3)
!   end do
!   write(*,*)


end function check_seen


!from stackexchange
subroutine system_mem_usage(valueRSS)
#ifdef __INTEL_COMPILER
use ifport !if on intel compiler
#endif

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


end module lookup_module_hash

