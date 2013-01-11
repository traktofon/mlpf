! vim: set ts=3 sw=3 :
module modeutil

   implicit none

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

end module modeutil
