!############################ subroutine read_input #################################
!####################################################################################
subroutine read_input(pah)
!
! read the geometry of a polycyclic benzenoid structure pah from the file 'geometry',
! filter out only the carbon atoms, and create the topological matrix for it
!
  use types_module
  implicit none
  integer(kind=4) :: info
  integer(kint) :: cnat=0,status,bnat,i,j,k,nhex,l,m,errorcode,a1,a2
  integer(kint),allocatable,dimension(:,:) :: lista
  integer(kint),allocatable,dimension(:,:) :: localbondlist
  character(len=2) :: atname
  real(kreal),parameter :: ccdist=1.7d0
  real(kreal) :: inertia(3,3),eival(3),work(100)
  real(kreal),allocatable,dimension(:,:) :: geom
  integer(kint),allocatable,dimension(:) :: map
  real(kreal),dimension(3) :: x
  real(kreal) :: dist
  type(structure) :: pah
  logical :: inlist,bondfileexists,adjexists

  inquire(file='adjlist',exist=adjexists)
  if (adjexists) then
    open(20,file='adjlist')
    read(20,*)cnat
  else 
  

! ######################
! # read geometry file #
! ######################
    open(20,file='geometry')
    read(20,*)bnat
    read(20,*)
    allocate(geom(3,bnat))
    allocate(map(bnat))
    map=0
    do i=1,bnat
      read(20,*)atname,(x(j),j=1,3)
      if (atname=='C' .or. atname=='c') then
        cnat=cnat+1
        map(i)=cnat
        do j=1,3
         geom(j,cnat)=x(j)
        end do
      end if
    end do
    close(20)
  end if

! ######################################
! # verify the maximal number of atoms #
! ######################################
  if (cnat > maxatoms) then
    write(*,*)"Error"
    write(*,*)"Recompile your code with larger value of the maxatoms in types_module.f90"
    write(*,*)"Minimal value of maxatoms you need is:",cnat
    stop
  end if

! ######################
! # allocate structure #
! ######################
  pah%nat=cnat
  pah%order=0
  allocate(pah%initiallabel(pah%nat))
  allocate(pah%neighbornumber(pah%nat))
  allocate(pah%neighborlist(pah%nat,3))
  pah%neighbornumber=0

  if (adjexists) then
    do i=1,cnat-1
      read(20,*)atname,(pah%neighborlist(i,k),k=1,3)
      pah%neighbornumber(i)=3
      do k=1,3
        if (pah%neighborlist(i,k).eq.0) pah%neighbornumber(i)=pah%neighbornumber(i)-1
      end do
    end do
    pah%neighbornumber(cnat)=3

    j=1
! find connectivity of the last atom
    do i=1,cnat-1
      do k=1,3
        if (pah%neighborlist(i,k).eq.cnat) then
          pah%neighborlist(cnat,j)=i
          j=j+1
        end if
      end do
    end do
    do k=1,3
        if (pah%neighborlist(cnat,k).eq.0) pah%neighbornumber(cnat)=pah%neighbornumber(cnat)-1
    end do

  else
! #######################
! # find neighbor table #
! #######################
    do i=1,cnat
      pah%initiallabel(i)=i
      do j=i+1,cnat
        if (dist(cnat,i,j,geom) < ccdist) then
          pah%neighbornumber(i)=pah%neighbornumber(i)+1
          pah%neighborlist(i,pah%neighbornumber(i))=j
          pah%neighbornumber(j)=pah%neighbornumber(j)+1
          pah%neighborlist(j,pah%neighbornumber(j))=i
        end if
      end do
    end do
  end if
!  do i=1,cnat
!    write(*,'(4(I5,2x))'),i,(pah%neighborlist(i,k),k=1,3)
!  end do

  if (.not. adjexists) then
! ################################################
! # construct Schlegel diagram for the fullerene #
! ################################################
    call schlegel_diagram(pah%nat,pah,geom)

! ####################################################
! # read (if provided) the preferred partition order #
! ####################################################
    pah%nbondlistentries=0
    inquire(file='bondlist',exist=bondfileexists)
    if (bondfileexists) then
      allocate(localbondlist(2,cnat))
      open(21,file='bondlist')
      do
        read(21,*,iostat=errorcode)a1,a2
        if (errorcode == 0) then
          pah%nbondlistentries=pah%nbondlistentries+1
          localbondlist(1,pah%nbondlistentries)=map(a1)
          localbondlist(2,pah%nbondlistentries)=map(a2)
          if (a1 > bnat .or. a2 > bnat) then
            write(*,*)"Ooops, looks like your file: bondlist"
            write(*,*)"does not correspond to your input file"
            write(*,*)"cnat",cnat," a1,a2",a1,a2
            stop
          end if
        else if (errorcode == -1) then
          exit
        else
          write(*,*)"Ooops, reading error from the file: bondlist"
          write(*,*)"Verify if the file is not corrupted"
          stop
        end if
      end do
      allocate(pah%bondlist(2,pah%nbondlistentries))
      pah%bondlist=localbondlist(:,1:pah%nbondlistentries)
      deallocate(localbondlist)
      close(21)
      call clean_bond_list(pah)
    end if


! #########################################################
! # find all substructures with removed one aromatic ring #
! #########################################################
    allocate(lista(6,pah%nat))
    call find_all_hexagons(pah%nat,pah,nhex,lista)
    open(22,file='new.geometries')
    do i=1,nhex
      inertia=0.0d0
      do l=1,pah%nat
        inlist=.false.
        do m=1,6
          if (l == lista(m,i)) inlist=.true.
        end do
        if (inlist) cycle
        do j=1,3
          do k=1,3
            inertia(j,k)=inertia(j,k)+geom(j,l)*geom(k,l)
          end do
        end do
      end do
      write(22,'(i6)')pah%nat
      write(22,'(1x,4i6,3f15.3)')pah%nat,0,0,0,inertia(1,1)+inertia(2,2)+inertia(3,3)
      do j=1,pah%nat
        inlist=.false.
        do k=1,6
          if (j == lista(k,i)) inlist=.true.
        end do
        if (inlist) then
          write(22,'(a,3f20.10)')"B",(geom(k,j),k=1,3) 
        else
          write(22,'(a,3f20.10)')"C",(geom(k,j),k=1,3) 
        end if
      end do
    end do
    close(22)

  end if 
  return

end subroutine read_input
!####################################################################################
!######################### end of subroutine read_input #############################
