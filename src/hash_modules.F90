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
use options_module
use ISO_FORTRAN_ENV
use, intrinsic :: iso_c_binding
  logical :: firstrun = .true.
  integer(int64) :: nstruct = 0
  integer(int64) :: nstructall = 0
  integer(int64) :: nstructseen = 0
  integer(int64) :: ncachebytes = 0
  type,public :: cachedstructure
     integer(C_signed_char), pointer :: hash(:)     
     integer(int16) :: order
     integer(int16) :: nat
     integer(kint) :: iseen
     integer(int64) :: lastseen
     integer(kint) :: mpacksize
     integer(kint),pointer :: packedpolynomial(:)
  end type cachedstructure

  type ptrcachedstructure
  type(cachedstructure), allocatable :: p(:)
  end type ptrcachedstructure

  type(ptrcachedstructure),allocatable :: x(:)
  integer(kint),allocatable :: xlen(:)
  integer(kint),allocatable :: irepl(:)
  interface
  end interface


contains 
subroutine add_cachedstructure(nat,nbnum,a,order,poly,hashseen,iseen,lastseen,duringread)
  use hash_module
  use options_module
  use types_module
  implicit none
  integer(kint),intent(in) :: nat
  integer(kint), intent(in) :: nbnum(nat)
  integer(kint), intent(in) :: a(3,nat)
  integer(int16) :: asmall(3,nat)
  integer(kint), intent(in) :: order
  type(vlonginteger), intent(in) :: poly(order+1)
  type(ptrcachedstructure), target :: xtmp
  logical :: first
  type(cachedstructure), allocatable :: ptmp(:)
  integer(int64) :: idx1l
  real(8)        :: highscore,score
  character(len=1) :: buf (3*nat*2)
  integer :: i,j,ilen,idx1,ihighscore
  integer(C_signed_char) :: hashsum(hashsize)
  integer(C_signed_char),OPTIONAL :: hashseen(hashsize)
  integer(kint),OPTIONAL :: iseen
  integer(int64),OPTIONAL :: lastseen
  logical,OPTIONAL :: duringread

  logical :: replace
  integer :: mpacksize

  if (.not. present(duringread)) then ! only print in main calculation not when reading cache.bin
    nstructall=nstructall+1
    if (mod(nstructall,100000_int64).eq.0)  then
      if (verbose) write(*,'(A,I0,1X,I0,A,I0)')'nstruct ',nstructall,nstruct,' filesize ',ncachebytes
    end if 
  end if


  
  asmall=0 ! not used connections filled with zeroes

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

  idx1=mod(int(idx1l,kind=kint),nbuckets)
  idx1=iabs(idx1)+1


  xtmp=x(idx1)

  

  if (xlen(idx1)>0) then ! xlen is zeroed at init of check_seen
     first=.false.
  else
     first=.true.
     irepl(idx1)=1
  end if
  
  ilen=xlen(idx1)
  if (nstruct.gt.maxrecords .and. .not. first) then
    replace=.true.
  else
    replace=.false.
  end if

  if (.not. replace) then

    if (ilen>0) then 
      if ( ilen == size(x(idx1)%p)) then ! no room in current array, needs enlarging by chunksize
        allocate(ptmp(ilen+chunksize))
        ptmp(1:size(x(idx1)%p)) = x(idx1)%p
        call move_alloc(ptmp,x(idx1)%p)
      end if
    else  
      allocate(x(idx1)%p(ilen+chunksize))
    end if
   
    allocate(x(idx1)%p(ilen+1)%hash(hashsize))
    
    

    do i=1,hashsize
      x(idx1)%p(ilen+1)%hash(i)=hashsum(i)
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

    mpacksize=getpackedsize(poly,order+1)
    x(idx1)%p(ilen+1)%mpacksize = mpacksize
    allocate(x(idx1)%p(ilen+1)%packedpolynomial(mpacksize))
    
    call packvliarray2(poly,x(idx1)%p(ilen+1)%packedpolynomial,order+1,mpacksize)
  
    xlen(idx1)=xlen(idx1)+1
    nstruct=nstruct+1
    ncachebytes=ncachebytes+8*2+2*2+4+8+hashsize+mpacksize*kint

    if (nstruct.eq.maxrecords .and. verbose) write(*,*)'Max records',nstruct,' achieved' 
else ! replace

! find best candidate to replace

highscore=0
ihighscore=0
do j=1,xlen(idx1)
  score=nstructall-x(idx1)%p(j)%lastseen
  score=score/(x(idx1)%p(j)%iseen+1)
  score=score/(x(idx1)%p(j)%nat**2) 
  if (score.gt.highscore .or. j.eq.1) then
     highscore=score
     ihighscore=j
  end if
end do
  j=ihighscore
  do i=1,hashsize
      x(idx1)%p(j)%hash(i)=hashsum(i)
  end do
 
  mpacksize=getpackedsize(poly,order+1)
  ncachebytes=ncachebytes+(mpacksize-x(idx1)%p(j)%mpacksize)*kint

  if (x(idx1)%p(j)%order .ne. order .or. x(idx1)%p(j)%mpacksize .ne. mpacksize) then
     deallocate(x(idx1)%p(j)%packedpolynomial) 
     allocate(x(idx1)%p(j)%packedpolynomial(mpacksize))
  end if
   x(idx1)%p(j)%mpacksize = mpacksize
  call packvliarray2(poly,x(idx1)%p(j)%packedpolynomial,order+1,mpacksize)
 

  x(idx1)%p(j)%order=order
  x(idx1)%p(j)%iseen=0
  x(idx1)%p(j)%nat=nat
  x(idx1)%p(j)%lastseen=nstructall
 
  irepl(idx1)=j
end if

  if (.not.present(hashseen) .and. .not. present(duringread)) then
   if (mod(nstructall,writemark).eq.0)  then
    if (has_write_cache_file) then   
      if (verbose) write(*,*)'Saving cache',trim(write_cache_fname)
      call writetodisk(write_cache_fname)
    end if 
   end if
  end if

end subroutine add_cachedstructure


subroutine writetodisk(fname)
use hash_module
use options_module
implicit none
  type(ptrcachedstructure), target :: xtmp
  logical :: first
  integer :: i,ilen,idx1,ibuf
  character(len=*), intent(in) :: fname

  open(23,file=trim(fname),FORM='UNFORMATTED')
  write(23)vlongmax,nstructall

  do idx1=1,nbuckets
    xtmp=x(idx1)
    if (allocated(x(idx1)%p)) then 
       first=.false.
       ilen=xlen(idx1)
    else
       first=.true.
       ilen=0
    end if
    do i=1,ilen
      write(23)x(idx1)%p(i)%order,x(idx1)%p(i)%iseen,x(idx1)%p(i)%nat,x(idx1)%p(i)%lastseen,&
               (x(idx1)%p(i)%hash(ibuf),ibuf=1,hashsize) !&
      write(23)(x(idx1)%p(i)%packedpolynomial(ibuf),ibuf=1,x(idx1)%p(i)%mpacksize)

    end do
  end do
  if (verbose) write (*,*)'cache saved',ncachebytes+20,'bytes'
  close(23)
end subroutine writetodisk


subroutine readfromdisk(fname)
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
  character(len=*), intent(in) :: fname
  integer :: ibuf
  
   allocate(x(nbuckets),xlen(nbuckets),irepl(nbuckets))
   firstrun=.false.
   xlen=0
  if (verbose) write(*,*)'Reading cache'
  open(23,file=trim(fname),FORM='UNFORMATTED')
  read(23)vlong,nstructall
  if (vlong.gt.vlongmax) then
    write(*,*)'WARNING!'
    write(*,*)'Size of vlong on file is larger than in the current code'
    write(*,*)'The readfromfile procedure may not catch corruption'
    write(*,*)'Inspect the results carefully'
  end if
  ires=0
  do while (ires.eq.0)
  read(23,IOSTAT=ires)order,iseen,nat,lastseen,(hashsum(ibuf),ibuf=1,hashsize)
  if (ires.eq.0) then
   allocate(poly(order+1))
   order32=order
   nat32=nat
   call readfromfile(23,poly,order32+1)
   call add_cachedstructure(nat32,nbnum,a,order32,poly,hashsum,iseen,lastseen,.true.)
   deallocate(poly)
  end if
  end do
  close(23)  
  if (verbose) write(*,*)'cache processed'

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
  integer :: i,j,idx1
  type(ptrcachedstructure),pointer :: curr
  type(cachedstructure), pointer :: match
  type(ptrcachedstructure), target :: xtmp
  integer(int64) :: idx1l
  character(len=1) :: buf (3*nat*2)
  integer(C_signed_char) :: hashsum(hashsize)

  nstructseen = nstructseen + 1

  if (firstrun) then
    allocate(x(nbuckets),xlen(nbuckets),irepl(nbuckets))
    firstrun=.false.
    xlen=0
  end if

  
  asmall=0 ! zeroes parts unbonded parts
  do i=1,nat
   do j=1,nbnum(i)
     asmall(j,i)=a(j,i)
   end do
  end do

  buf=transfer(asmall(1:3,1:nat),buf)

  call hash(buf,size(buf,1,C_long),hashsum)
  idx1l=transfer(hashsum,idx1l)


  idx1=mod(int(idx1l, kind=kint),nbuckets)

  idx1=iabs(idx1)+1


  seen = .false.
  xtmp=x(idx1)
  curr=>xtmp
  structloop:   do i=1,xlen(idx1)
    do j=1,hashsize
        if (curr%p(i)%hash(j).ne. hashsum(j)) then 
          cycle structloop
        end if
      if (j.eq.hashsize) then 
        seen=.true.
        match=>curr%p(i)
        exit structloop
      end if
    end do
 end do structloop
 if (seen) then 
    order=match%order
    x(idx1)%p(i)%iseen=x(idx1)%p(i)%iseen+1
    x(idx1)%p(i)%lastseen=nstructall
    allocate(poly(order+1))
    call unpackvliarray2(match%packedpolynomial,poly,order+1,match%mpacksize)

 end if


end function check_seen

end module lookup_module_hash
