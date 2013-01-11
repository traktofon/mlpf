! vim: set ts=3 sw=3 :
module modeutil

   implicit none

   type :: modecomb_t
      integer :: imode ! first mode in the combination
      integer :: fmode ! last mode in the combination
   end type modecomb_t

   contains

   !--------------------------------------------------------------------
   subroutine vgn_shape(m,gdim,vd,gd,nd)
   !--------------------------------------------------------------------
      implicit none
      integer,intent(in)  :: m
      integer,intent(in)  :: gdim(:)
      integer,intent(out) :: vd,gd,nd
      integer             :: nmodes
      nmodes=size(gdim)
      vd=1
      nd=1
      gd=gdim(m)
      if(m>1)      vd=product(gdim(1:m-1))
      if(m<nmodes) nd=product(gdim(m+1:nmodes))
   end subroutine vgn_shape


   !--------------------------------------------------------------------
   subroutine make_mc_shape(gdim,mc,np,mdim)
   ! Given a shape (gdim) and a list of mode combinations to be
   ! performed (mc), returns the new number of modes (np) and the new
   ! shape (mdim).
   !--------------------------------------------------------------------
      implicit none
      integer,intent(in)          :: gdim(:)
      type(modecomb_t),intent(in) :: mc(:)
      integer,intent(out)         :: np
      integer,intent(out)         :: mdim(:)
      integer :: nmodes,m,m1,m2,nmc,c,p

      nmodes = size(gdim)
      nmc    = size(mc)
      c  = 1 ! counting mode combinations
      p  = 1 ! counting new modes
             ! m1:m2 = modes that will not be combined
      m1 = 1
   10 continue
      if (c <= nmc) then
         ! if there is another combination to be performed,
         ! m2 is the mode before the first mode in that combination
         m2 = mc(c)%imode - 1
      else
         ! otherwise m2 is the last mode
         m2 = nmodes
      endif
      ! Copy dimensions of uncombined modes.
      do m=m1,m2
         mdim(p) = gdim(m)
         p = p+1
      enddo
      ! If no more combinations are to be performed, exit.
      if (c > nmc) goto 20
      ! Otherwise find the dimension of the combined mode =
      ! the product of all dimensions of the modes that get combined
      mdim(p) = product(gdim(mc(c)%imode : mc(c)%fmode))
      p = p+1
      ! Now point m1 to the first mode after this combination.
      m1 = mc(c)%fmode + 1
      ! Move to the next combination.
      c = c+1
      ! Any old modes left? Then continue.
      if (m1 <= nmodes) goto 10
   20 continue
      ! All combinations done, set new number of modes.
      np = p-1
      return
   end subroutine make_mc_shape


end module modeutil
