!##################### subroutine create_nobond_daughter ############################
!####################################################################################
recursive subroutine create_nobond_daughter(pah,bond,atom1,atom2)
!
! creates a daughter structure (bond) from a parent structure (pah)
! by deleting a bond between atoms: atom1 and atom2
! 
  use types_module
  implicit none
  integer(kint) :: atom1,atom2,j,i
  type(structure) :: pah,bond

! #####################################
! # initialize the daughter structure #
! #####################################
  bond%nat=pah%nat
  bond%order=0
  bond%nbondlistentries=pah%nbondlistentries
  allocate(bond%neighbornumber(bond%nat))
  allocate(bond%neighborlist(3,bond%nat))
  if (bond%nbondlistentries > 0) then
    allocate(bond%bondlist(2,bond%nbondlistentries))
    bond%bondlist=pah%bondlist
  end if
  bond%neighbornumber=0
  bond%neighborlist=0

! ###############################
! # fill the daughter structure #
! ###############################
  bond%neighbornumber=pah%neighbornumber
  bond%neighbornumber(atom1)=bond%neighbornumber(atom1)-1
  bond%neighbornumber(atom2)=bond%neighbornumber(atom2)-1
  bond%neighborlist=pah%neighborlist
  bond%neighborlist(1:3,atom1)=0
  bond%neighborlist(1:3,atom2)=0
  j=0
  do i=1,3
    if (pah%neighborlist(i,atom1) /= atom2) then
      j=j+1
      bond%neighborlist(j,atom1)=pah%neighborlist(i,atom1)
    end if
  end do
  j=0
  do i=1,3
    if (pah%neighborlist(i,atom2) /= atom1) then
      j=j+1
      bond%neighborlist(j,atom2)=pah%neighborlist(i,atom2)
    end if
  end do
  call clean_bond_list(bond)
  return

end subroutine create_nobond_daughter
!####################################################################################
!################## end of subroutine create_nobond_daughter ########################




!######################## subroutine create_noatoms_daughter ########################
!####################################################################################
recursive subroutine create_noatoms_daughter(pah,pah1,nelim,delatoms)
!
! creates a daughter structure (pah1) from a parent structure (pah)
! by deleting nelim atoms: delatom(1),...,delatom(nelim)
! 
  use types_module
  implicit none
  integer(kint):: delatoms(nelim)
  integer(kint) :: j,i,k,l,m,nelim
  integer(kint),dimension(maxatoms) :: map
  type(structure) :: pah,pah1
  logical,dimension(maxatoms) :: offlist

  map(1:pah%nat)=0
  offlist(1:pah%nat)=.true.

! ########################################  
! # logical vector of atoms for deleting #
! ########################################  
  forall (i=1:nelim)
    offlist(delatoms(i))=.false.
  end forall

! #####################################
! # initialize the daughter structure #
! #####################################
  pah1%nat=pah%nat-nelim
  pah1%order=0
  pah1%nbondlistentries=pah%nbondlistentries
  allocate(pah1%neighbornumber(pah1%nat))
  allocate(pah1%neighborlist(3,pah1%nat))
  if (pah1%nbondlistentries > 0) then
    allocate(pah1%bondlist(2,pah1%nbondlistentries))
  end if
  pah1%neighbornumber=0
  pah1%neighborlist=0

! #############################
! # create the transition map #
! #############################
  j=0
  do i=1,pah%nat
    if (offlist(i)) then
      j=j+1
      map(i)=j
    end if
  end do

! ###########################################
! # create the structure with deleted atoms #
! ###########################################
  do i=1,pah%nat
    if (map(i) /= 0) then 
      k=0
      do l=1,pah%neighbornumber(i)
        m=map(pah%neighborlist(l,i))
        if (m /= 0) then
          k=k+1
          pah1%neighborlist(k,map(i))=m
        end if
      end do
      pah1%neighbornumber(map(i))=k
    end if
  end do

! #######################################
! # translate to the new atom numbering #
! #######################################
  forall (i=1:pah1%nbondlistentries, j=1:2)
      pah1%bondlist(j,i)=map(pah%bondlist(j,i))
  end forall
  call clean_bond_list(pah1)
  return

end subroutine create_noatoms_daughter
!####################################################################################
!#################### end of subroutine create_noatoms_daughter ######################
