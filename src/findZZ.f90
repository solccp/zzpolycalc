!###################### subroutine find_ZZ_polynomial ###############################
!####################################################################################
recursive subroutine find_ZZ_polynomial(pah,level,path)
! 
! find recursively the ZZ polynomial for the structure pah
!
  use types_module
  use lookup_module_hash
  implicit none
  integer(kint) :: medat,level,path
  type(structure) :: pah
  logical :: seen

!   write(*,*)'Level',level,pah%nat,pah%order


! ###########################
! # if pah contains 0 atoms #
! ###########################
  if (pah%nat == 0) then
    pah%order=0
    allocate(pah%polynomial(pah%order+1))
    pah%polynomial(1)=setvli(1)

! ##########################################
! # if pah contains an odd number of atoms #
! ##########################################
  else if (mod(pah%nat,2) == 1) then
    pah%order=0
    allocate(pah%polynomial(pah%order+1))
    pah%polynomial(1)=setvli(0)

! ##########################################
! # check if the graph of pah is connected #
! ##########################################
  else

    call check_if_connected(pah,medat)

!   #############################################
!   # if pah is connected, decompose it further #
!   #############################################
    if (medat == 0) then

    seen=check_seen(pah%nat,pah%neighbornumber,pah%neighborlist,pah%order,pah%polynomial)

    if (.not. seen) then
      call decompose(pah,level,path)

      call add_neigh(pah%nat,pah%neighbornumber,pah%neighborlist,pah%order,pah%polynomial)
    end if


!   ########################################################################
!   # if pah is disconnected, split it and decompose the fragments further #
!   ########################################################################
    else
      call split_and_decompose(pah,medat,level,path)
    end if

  end if
  return

end subroutine find_ZZ_polynomial
!####################################################################################
!################## end of subroutine find_ZZ_polynomial ############################
