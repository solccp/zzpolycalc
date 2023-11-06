!######################## subroutine schlegel_diagram ###############################
!####################################################################################
subroutine schlegel_diagram(nat,pah,geom)
!
! construct a planar graph for a fullerene containing at least one pentagon;
! if the fullerene contains more pentagons, all Schlegel diagrams are constructed
! the construction is based on a stereographic projection of the fullerene cage
! onto the 2D surface, assuming that the nadir is located in the center of the selected pentagon
!
  use types_module
  implicit none
  integer(kint) :: nat,i,j,k,npent,l
  integer(kint),allocatable,dimension(:,:) :: lista5
  real(kreal) :: cpnorm,umost,dmost,rmost,lmost,xscal,yscal
  real(kreal),dimension(3,nat) :: geom
  real(kreal),allocatable,dimension(:,:) :: geom1,xyst
  real(kreal),allocatable,dimension(:) :: theta,phi
  real(kreal),dimension(3) :: tpmap,xv,yv,cp,cm
  type(structure) :: pah
  character(len = 6) :: ftitle

! ############################
! # find center of fullerene #
! ############################
  do j=1,3
    cm(j)=0.0d0
    do i=1,pah%nat
      cm(j)=cm(j)+geom(j,i)
    end do
    cm(j)=cm(j)/pah%nat
  end do

! ###################################
! # unit vector to each carbon atom #
! ###################################
  allocate(theta(pah%nat))
  allocate(phi(pah%nat))
  allocate(xyst(2,pah%nat))
  allocate(geom1(3,pah%nat))
  do i=1,pah%nat
    cpnorm=0.0d0
    do j=1,3
      cpnorm=cpnorm+(geom(j,i)-cm(j))**2
    end do
    cpnorm=dsqrt(cpnorm)
    do j=1,3
      geom1(j,i)=(geom(j,i)-cm(j))/cpnorm
    end do
  end do
  
! ######################
! # find all pentagons #
! ######################
  allocate(lista5(5,pah%nat))
  call find_all_pentagons(pah%nat,pah,npent,lista5)

! #######################
! # loop over pentagons #
! #######################
  do i=1,npent
    
!   ###########################
!   # find center of pentagon #
!   ###########################
    do j=1,3
      cp(j)=0.0d0
      do k=1,5
        cp(j)=cp(j)+geom(j,lista5(k,i))
      end do
      cp(j)=cp(j)/5.0d0
    end do

!   #########################################
!   # unit vector to the center of pentagon #
!   #########################################
    cpnorm=0.0d0
    do j=1,3
      cpnorm=cpnorm+(cp(j)-cm(j))**2
    end do
    cpnorm=dsqrt(cpnorm)
    do j=1,3
      cp(j)=(cp(j)-cm(j))/cpnorm
    end do

!   ##########################################################################
!   # unit vector perpendicular to the unit vector to the center of pentagon #
!   ##########################################################################
    j=1
    call radial_neighbor(nat,pah,lista5,i,j,k)
    j=2
    call radial_neighbor(nat,pah,lista5,i,j,l)
    do j=1,3
      xv(j)=geom1(j,k)-geom1(j,l)
    end do
    cpnorm=dsqrt(xv(1)*xv(1)+xv(2)*xv(2)+xv(3)*xv(3))
    do j=1,3
      xv(j)=xv(j)/cpnorm
    end do
    cpnorm=xv(1)*cp(1)+xv(2)*cp(2)+xv(3)*cp(3)
    do j=1,3
      xv(j)=xv(j)-cp(j)*cpnorm
    end do
!    proj=dsqrt(geom1(1,k)*geom1(1,k)+geom1(2,k)*geom1(2,k)+geom1(3,k)*geom1(3,k))    
!    cpnorm=cp(1)*geom1(1,k)+cp(2)*geom1(2,k)+cp(3)*geom1(3,k)
!    xv(1)=geom1(1,k)-cpnorm/proj*cp(1)    
!    xv(2)=geom1(2,k)-cpnorm/proj*cp(2)    
!    xv(3)=geom1(3,k)-cpnorm/proj*cp(3)    
    cpnorm=dsqrt(xv(1)*xv(1)+xv(2)*xv(2)+xv(3)*xv(3))
    do j=1,3
      yv(j)=xv(j)/cpnorm
    end do

!   #####################
!   # third unit vector #
!   #####################
    xv(1)=cp(2)*yv(3)-cp(3)*yv(2)
    xv(2)=cp(3)*yv(1)-cp(1)*yv(3)
    xv(3)=cp(1)*yv(2)-cp(2)*yv(1)

!   #################################################
!   # theta and phi coordinates of each carbon atom #
!   # and stereographic projection on XY plane      #
!   #################################################
    do j=1,nat
      tpmap(1)=cp(1)*geom1(1,j)+cp(2)*geom1(2,j)+cp(3)*geom1(3,j)
      tpmap(2)=xv(1)*geom1(1,j)+xv(2)*geom1(2,j)+xv(3)*geom1(3,j)
      tpmap(3)=yv(1)*geom1(1,j)+yv(2)*geom1(2,j)+yv(3)*geom1(3,j)
      xyst(1,j)=tpmap(3)/(1.049-tpmap(1))
      xyst(2,j)=tpmap(2)/(1.049-tpmap(1))
    end do

!   #################################
!   # find extremal points of graph #
!   #################################
    umost=0.0
    dmost=0.0
    lmost=0.0
    rmost=0.0
    do j=1,5
      call radial_neighbor(nat,pah,lista5,i,j,k)
      if (xyst(1,k) .gt. rmost) rmost=xyst(1,k)
      if (xyst(1,k) .lt. lmost) lmost=xyst(1,k)
      if (xyst(2,k) .gt. umost) umost=xyst(2,k)
      if (xyst(2,k) .lt. dmost) dmost=xyst(2,k)
    end do 
    yscal=3*(1+0.809)/(umost-dmost)
    xscal=6*(0.951)/(rmost-lmost)


!   #####################################
!   # rescale all graph point and shift #
!   #####################################
    do j=1,nat
      xyst(1,j)=(xyst(1,j)-rmost)*xscal+3*0.951
      xyst(2,j)=(xyst(2,j)-umost)*yscal+3
    end do


!   ###########################################
!   # update geometry of the leading pentagon #
!   ###########################################
!    xyst(1,lista5(4,i))=0.0
!    xyst(2,lista5(4,i))=3.0
!    xyst(1,lista5(5,i))=-3*0.951
!    xyst(2,lista5(5,i))=3*0.309
!    xyst(1,lista5(2,i))=-3*0.588
!    xyst(2,lista5(2,i))=-3*0.809
!    xyst(1,lista5(1,i))=3*0.588
!    xyst(2,lista5(1,i))=-3*0.809
!    xyst(1,lista5(3,i))=3*0.951
!    xyst(2,lista5(3,i))=3*0.309

    do j=1,5
      call radial_neighbor(nat,pah,lista5,i,j,k)
      xyst(1,lista5(j,i))=xyst(1,k)/(dble(nat-5))*(dble(nat))*1.1
      xyst(2,lista5(j,i))=xyst(2,k)/(dble(nat-5))*(dble(nat))*1.1
    end do




!   ##################################
!   # write schlegel diagram to file #
!   ##################################
    write(ftitle(1:4),'(a4)')'pent'
    write(ftitle(5:6),'(i2.2)')i
    open(21,file=ftitle)

!   ##################
!   # write vertices #
!   ##################
    do j=1,nat
      write(21,'(2f12.6)')xyst(1,j),xyst(2,j)
    end do
    write(21,*)" "

!   ###############
!   # write edges #
!   ###############
    do j=1,nat
      do k=1,3
        if (pah%neighborlist(k,j) .gt. j) then
          write(21,'(2f12.6)')xyst(1,j),xyst(2,j)
          write(21,'(2f12.6)')xyst(1,pah%neighborlist(k,j)),xyst(2,pah%neighborlist(k,j))
          write(21,*)" "
        end if
      end do
    end do
    close(21)

  
!  print*,cp
!  print*,xv
!  print*,yv


  end do

  return

end subroutine schlegel_diagram
!####################################################################################
!################### end of subroutine schlegel_diagram #############################


!########################  subroutine radial neighbor ###############################
!####################################################################################
subroutine radial_neighbor(nat,pah,lista5,pn,cn,rn)
!
! find the neigbor of the atom cn in pentagon pn and return the results as rn
! all atoms of pentagon pn are stored in lista5(1..5,pn)
! particularly, the index of the atom cn in pentagon pn is lista5(cn,pn)
! rn is the only neighbor of lista5(1..5,pn) not located in pentagon pn
!
  use types_module
  implicit none
  integer(kint) :: nat,pn,cn,rn,i,j,istat
  integer(kint),dimension(5,nat) :: lista5
  type(structure) :: pah

  do i=1,3
    istat=1
    do j=1,5
      if (pah%neighborlist(i,lista5(cn,pn)) .eq. lista5(j,pn)) istat=0
    end do
    if (istat .eq. 1) then
      rn=pah%neighborlist(i,lista5(cn,pn))
    end if
  end do

  return

end subroutine radial_neighbor

!####################################################################################
!##################### end of subroutine radial neighbor ############################
