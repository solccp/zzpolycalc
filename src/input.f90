!############################ subroutine read_input #################################
!####################################################################################
subroutine read_input(input_fname,pah)
!
! read the geometry of a polycyclic benzenoid structure pah from the file 'geometry',
! filter out only the carbon atoms, and create the topological matrix for it
!
  use types_module
  use options_module
  implicit none
  integer(kint) :: cnat=0,bnat,i,j,k,errorcode,a1,a2
  integer(kint),allocatable,dimension(:,:) :: localbondlist
  character(len=200) :: atname
  real(kreal),parameter :: ccdist=1.7d0
  real(kreal),allocatable,dimension(:,:) :: geom
  integer(kint),allocatable,dimension(:) :: map
  real(kreal),dimension(3) :: x
  real(kreal) :: dist
  type(structure), intent(inout) :: pah
  logical :: bondfileexists
  character(len=*), intent(in) :: input_fname

!connection table
  integer :: nConnection
  integer :: con1, con2



  if (is_adjacencyfile) then
    open(20,file=trim(input_fname),status='old')
    read(20,*)cnat
  else 
! ######################
! # read geometry file #
! ######################
    open(20,file=trim(input_fname),status='old')
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
    if (.not. unsorted_geometry) call sort(geom,cnat)
  end if ! is_adjacency_file

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
  allocate(pah%neighborlist(3,pah%nat))
  pah%neighbornumber=0
  pah%neighborlist=0

  if (is_adjacencyfile) then
    do i=1,cnat-1
      read(20,*)atname,(pah%neighborlist(k,i),k=1,3)
      pah%neighbornumber(i)=3
      do k=1,3
        if (pah%neighborlist(k,i).eq.0) pah%neighbornumber(i)=pah%neighbornumber(i)-1
      end do
    end do
    pah%neighbornumber(cnat)=3

    j=1
! find connectivity of the last atom
    do i=1,cnat-1
      do k=1,3
        if (pah%neighborlist(k,i).eq.cnat) then
          pah%neighborlist(j,cnat)=i
          j=j+1
        end if
      end do
    end do
    do k=1,3
        if (pah%neighborlist(k,cnat).eq.0) pah%neighbornumber(cnat)=pah%neighbornumber(cnat)-1
    end do

  elseif (.not. read_connection_table ) then
! #######################
! # find neighbor table #
! #######################
      do i=1,cnat
        pah%initiallabel(i)=i
        do j=i+1,cnat
          if (dist(cnat,i,j,geom) < ccdist) then
            pah%neighbornumber(i)=pah%neighbornumber(i)+1
            pah%neighborlist(pah%neighbornumber(i),i)=j
            pah%neighbornumber(j)=pah%neighbornumber(j)+1
            pah%neighborlist(pah%neighbornumber(j),j)=i
          end if
        end do
      end do
    else
        nConnection = 0
        read (20, *) nConnection
        do i=1, nConnection
            read (20,*) con1, con2
            if (dist(cnat,con1,con2,geom) < ccdist) then
                pah%neighbornumber(con1)=pah%neighbornumber(con1)+1
                pah%neighborlist(pah%neighbornumber(con1),con1)=con2
                pah%neighbornumber(con2)=pah%neighbornumber(con2)+1
                pah%neighborlist(pah%neighbornumber(con2),con2)=con1
            end if
        end do
    end if
    
    close(20)
  
  if (.not. is_adjacencyfile) then
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

  end if 
  
  return

end subroutine read_input
!####################################################################################
!######################### end of subroutine read_input #############################
