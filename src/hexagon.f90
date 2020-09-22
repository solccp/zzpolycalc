!####################### subroutine find_all_hexagons ###############################
!####################################################################################
subroutine find_all_hexagons(nat,pah,nhex,lista)
! 
! find a list of all hexagons in a polycyclic structure pah
!
  use types_module
  implicit none
  integer(kint) :: i,j,k,l,atom2,atom3,nat
  integer(kint),dimension(6) :: sextet
  integer(kint) :: nhex ! number of hexagons in structure pah
  integer(kint),dimension(6,nat) :: lista
  type(structure) :: pah
  logical :: ring_exists

  nhex=0
! #######################
! # loop over all atoms #
! #######################
  atomloop: do i=1,pah%nat

!   ##################################
!   # if atom i has only 2 neighbors #
!   ##################################
    if (pah%neighbornumber(i) == 2) then
      if (pah%neighborlist(1,i) < i) cycle atomloop
      if (pah%neighborlist(2,i) < i) cycle atomloop
      call find_aromatic_sextet(pah,sextet,i,pah%neighborlist(2,i),pah%neighborlist(1,i),ring_exists)
      if (ring_exists) then
        do j=4,6
          if (sextet(j) < i) cycle atomloop
        end do
        nhex=nhex+1
        do l=1,6
          lista(l,nhex)=sextet(l)
        end do
      end if

!   #############################
!   # if atom i has 3 neighbors #
!   #############################
    else if (pah%neighbornumber(i) == 3) then

!     ###################################################
!     # loop over all 2-combinations of three neighbors #
!     ###################################################
      innerloop: do j=1,3
        atom2=pah%neighborlist(mod(j,3)+1,i)
        if (atom2 < i) cycle innerloop
        atom3=pah%neighborlist(mod(j+1,3)+1,i)
        if (atom3 < i) cycle innerloop
        call find_aromatic_sextet(pah,sextet,i,atom2,atom3,ring_exists)
        if (ring_exists) then
          do k=4,6
            if (sextet(k) < i) cycle innerloop
          end do
          nhex=nhex+1
          do l=1,6
            lista(l,nhex)=sextet(l)
          end do
        end if
      end do innerloop
    end if
  end do atomloop
  return

end subroutine find_all_hexagons
!####################################################################################
!#################### end of subroutine find_all_hexagons ###########################


!####################### subroutine find_all_pentagons ##############################
!####################################################################################
subroutine find_all_pentagons(nat,pah,npent,lista)
! 
! find a list of all pentagons in a polycyclic structure pah
!
  use types_module
  implicit none
  integer(kint) :: i,j,k,l,atom2,atom3,nat
  integer(kint),dimension(5) :: pentagon
  integer(kint) :: npent ! number of pentagons in structure pah
  integer(kint),dimension(5,nat) :: lista
  type(structure) :: pah
  logical :: ring_exists

  npent=0
! #######################
! # loop over all atoms #
! #######################
  atomloop: do i=1,pah%nat

!   ##################################
!   # if atom i has only 2 neighbors #
!   ##################################
    if (pah%neighbornumber(i) == 2) then
      if (pah%neighborlist(1,i) < i) cycle atomloop
      if (pah%neighborlist(2,i) < i) cycle atomloop
      call find_pentagon(pah,pentagon,i,pah%neighborlist(2,i),pah%neighborlist(1,i),ring_exists)
      if (ring_exists) then
        do j=4,5
          if (pentagon(j) < i) cycle atomloop
        end do
        npent=npent+1
        do l=1,5
          lista(l,npent)=pentagon(l)
        end do
      end if

!   #############################
!   # if atom i has 3 neighbors #
!   #############################
    else if (pah%neighbornumber(i) == 3) then

!     ###################################################
!     # loop over all 2-combinations of three neighbors #
!     ###################################################
      innerloop: do j=1,3
        atom2=pah%neighborlist(mod(j,3)+1,i)
        if (atom2 < i) cycle innerloop
        atom3=pah%neighborlist(mod(j+1,3)+1,i)
        if (atom3 < i) cycle innerloop
        call find_pentagon(pah,pentagon,i,atom2,atom3,ring_exists)
        if (ring_exists) then
          do k=4,5
            if (pentagon(k) < i) cycle innerloop
          end do
          npent=npent+1
          do l=1,5
            lista(l,npent)=pentagon(l)
          end do
        end if
      end do innerloop
    end if
  end do atomloop
  return

end subroutine find_all_pentagons
!####################################################################################
!#################### end of subroutine find_all_pentagons ##########################
