!CRC32 from Fortran Wiki. Public domain
function crc32_hash(a,continue) result (crc_64)
use,intrinsic :: ISO_FORTRAN_ENV, only : int32,int64
implicit none
logical,intent(in),optional  :: continue
character(len=1),intent(in)  :: a(:)
integer(int64)               :: crc_64
integer(int32),save          :: crc
integer                      :: i
integer(int32),save          :: crc_table(0:255)
integer,save                 :: icalled=0
         integer :: j
         integer(int32) :: k

   if(present(continue))then
      if(continue .eqv. .false.)then
         crc=0_int32
      endif
   else
      crc=0_int32
   endif
   if(icalled.eq.0)then         ! on first call generate table and use table for speed
!      block
!         integer :: i, j
!         integer(int32) :: k
         do i = 0, 255
            k = i
            do j = 1, 8
               if (btest(k, 0)) then
                  k = ieor(shiftr(k, 1), -306674912_int32)
               else
                  k = shiftr(k, 1)
               endif
            enddo
            crc_table(i) = k
         enddo
!      end block
      icalled=1
   endif
   crc = not(crc)
   do i = 1, size(a)
      crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(a(i))), 255)))
   enddo
   crc = not(crc)
   crc_64=transfer([crc,0_int32],crc_64)
end function crc32_hash
