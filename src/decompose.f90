!############################ subroutine decompose ##################################
!####################################################################################
recursive subroutine decompose(pah,level)
!
! decompose the original parent structure (pah) into three daughter structures:
!   1. with one edge bond deleted (between atom1 and atom2)
!   2. with two edge atoms (atom1 & atom2) deleted
!   3. with one ring containing atom1 & atom2 deleted
! 
  use types_module
  implicit none
  integer(kint) :: atom1,atom2,atom3,i,level,nelim,j
  integer(kint),dimension(6) :: sextet1,sextet2
  integer(kint),dimension(2) :: atoms
  type(structure),intent(inout) :: pah
  type(structure) :: bond,corners,ring1,ring2
  logical :: ring1_exists,ring2_exists,are_neighbors

! #############################################
! # select an edge bond between atom1 & atom2 #
! #############################################
  atom1=1
  atom2=pah%neighborlist(1,1)

! ######################################################
! # check if the selected edge belongs to any hexagons #
! ######################################################
  ring1_exists=.false.
  ring2_exists=.false.
  if (pah%neighbornumber(atom1) == 2) then
    atom3=pah%neighborlist(2,1)
    call find_aromatic_sextet(pah,sextet1,atom1,atom2,atom3,ring1_exists)
  else if (pah%neighbornumber(atom1) == 3) then
    atom3=pah%neighborlist(2,1)
    call find_aromatic_sextet(pah,sextet1,atom1,atom2,atom3,ring1_exists)
    atom3=pah%neighborlist(3,1)
    call find_aromatic_sextet(pah,sextet2,atom1,atom2,atom3,ring2_exists)
  end if
!  if (ring1_exists) write(*,*)'Ring 1'
!  if (ring2_exists) write(*,*)'Ring 2'

! ##################################
! # create the daughter structures #
! ##################################
  call create_nobond_daughter(pah,bond,atom1,atom2)
  nelim=2
  atoms(1)=atom1
  atoms(2)=atom2
  call create_noatoms_daughter(pah,corners,nelim,atoms)
  if (ring1_exists) then
    nelim=6
    call create_noatoms_daughter(pah,ring1,nelim,sextet1)
  end if
  if (ring2_exists) then
    nelim=6
    call create_noatoms_daughter(pah,ring2,nelim,sextet2)
  end if

! ###################################################
! # eliminate dangling bonds in daughter structures #
! ###################################################
  call cut_dangling_bonds(bond)
  call cut_dangling_bonds(corners)
  if (ring1_exists) then
    call cut_dangling_bonds(ring1)
  end if
  if (ring2_exists) then
    call cut_dangling_bonds(ring2)
  end if

! ###############################################
! # find ZZ polynomials for daughter structures #
! ###############################################
  call find_ZZ_polynomial(bond,level+1)
  call find_ZZ_polynomial(corners,level+1)
  if (ring1_exists) then
    call find_ZZ_polynomial(ring1,level+1)
  end if
  if (ring2_exists) then
    call find_ZZ_polynomial(ring2,level+1)
  end if

! ###############################################
! # find ZZ polynomial for the parent structure #
! ###############################################
   call sum_polynomials(pah,bond,corners,ring1,ring1_exists,ring2,ring2_exists)
!   write(*,*)'Level',level,pah%nat,pah%order

!   write(*,'(A)',advance='no')"N "   
!   do i=1,pah%nat
!     write(*,'(3I3)', advance='no')(pah%neighborlist(j,i),j=1,3)
!   end do
!   write(*,*)
!  call print_ZZ_polynomial(pah)  

! ##################################
! # deallocate daughter structures #
! ##################################
  deallocate(bond%neighbornumber)
  deallocate(corners%neighbornumber)
  deallocate(bond%neighborlist)
  deallocate(corners%neighborlist)
  deallocate(bond%polynomial)
  deallocate(corners%polynomial)
  if (ring1_exists) then
    deallocate(ring1%neighbornumber)
    deallocate(ring1%neighborlist)
    deallocate(ring1%polynomial)
    if (ring1%nbondlistentries > 0) then
      deallocate(ring1%bondlist)
    end if
  end if
  if (ring2_exists) then
    deallocate(ring2%neighbornumber)
    deallocate(ring2%neighborlist)
    deallocate(ring2%polynomial)
    if (ring2%nbondlistentries > 0) then
      deallocate(ring2%bondlist)
    end if
  end if
  if (bond%nbondlistentries > 0) then
    deallocate(bond%bondlist)
  end if
  if (corners%nbondlistentries > 0) then
    deallocate(corners%bondlist)
  end if

end subroutine decompose
!####################################################################################
!######################## end of subroutine decompose ###############################
