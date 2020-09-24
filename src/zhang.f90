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
  use lookup_module_md5
  implicit none
  integer(kint) :: i,nhex,level=0,j,notused
  integer(kint),allocatable,dimension(:,:) :: lista
  type(structure) :: pah

! ############################################################
! # read initial geometry data and create topological matrix #
! ############################################################
  call read_input(pah)

! #############################################################
! # find recursively the ZZ polynomial of the given structure #
! #############################################################
  call find_ZZ_polynomial(pah,level)

! ###########################
! # print the ZZ polynomial #
! ###########################
  call print_ZZ_polynomial(pah)
  write(*,*)'Unique',nstruct
  notused=0
  do i=1,maxtab
    do j=1,xlen(i)
!      write(*,*)'iseen',i,j,x(i)%p(j)%iseen
       if (x(i)%p(j)%iseen .gt. 0) notused=notused+1
    end do
  end do

  write(*,*)'Not used:',notused

end
!####################################################################################
!###################### end of program zhang_polynomial #############################
