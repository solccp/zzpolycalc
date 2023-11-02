!######################### program zhang_polynomial #################################
!####################################################################################
program zhang_polynomial
!
! This program calculates the Zhang-Zhang polynomial for benzenoid structures;
! the ZZ polynomial contains beside other quantities also the Clar number,
! Clar count, and Kekule number of a given structure
!   
! Reference:
!   I. Gutman, B. Furtula, and A. T. Balaban
!   Polycyclic Aromatic Compounds 26 pp.17-35, 2006
!
  use types_module
  use lookup_module_hash
  use options_module
  use getopt_module
  implicit none
  integer(kint) :: i,nhex,level=0,j,notused
  integer(kint),allocatable,dimension(:,:) :: lista
  type(structure) :: pah
  logical :: cacheexists
  character(len=200) :: input_fname


#ifdef USE_XXHASH
  write(*,*)'xxhash used'
#else
#ifdef USE_SHA256
  write(*,*)'SHA256 used'
#else
  write(*,*)'MD5 used'
#endif
#endif


! ############################################################
! # read initial geometry data and create topological matrix #
! ############################################################

  call initialize_options() 
  call read_options(input_fname)
  call read_input(input_fname,pah)

  
  
! read cache from disk 
  inquire(file='cache.bin',exist=cacheexists)
  if (cacheexists) then
    call readfromdisk
  end if

! #############################################################
! # find recursively the ZZ polynomial of the given structure #
! #############################################################
  call find_ZZ_polynomial(pah,level,0)

! ###########################
! # print the ZZ polynomial #
! ###########################
  call print_ZZ_polynomial(pah)
  write(*,*)'Unique',nstruct
  notused=0
  do i=1,nbuckets
    do j=1,xlen(i)
!      write(*,*)'iseen',i,j,x(i)%p(j)%iseen
       if (x(i)%p(j)%iseen .gt. 0) notused=notused+1
    end do
  end do

  write(*,*)'Not used:',notused
!  call writetodisk

end
!####################################################################################
!###################### end of program zhang_polynomial #############################
