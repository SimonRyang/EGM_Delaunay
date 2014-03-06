module comfunc

use comtypes
use comparam
use comvar
use comfuncform

implicit none

contains

! --------------------------------------------------------------------------------------
function func_makegrid(x1,x2,n,c)

! builds grids, linear if scaling factor c=1, triple exponential for c=3

implicit none
real(prec)::grid(n),func_makegrid(n)
real(prec),intent(in)::x1,x2,c
integer,intent(in)::n

integer::i
real(prec)::scale

scale = x2-x1
grid(1) = x1
grid(n) = x2
do i = 2,n-1
	grid(i) = x1+scale*((i-1.0)/(n-1.0))**c
end do

func_makegrid = grid

end function func_makegrid
! --------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
subroutine sub_exopolation(xdash,hdash,cdash,kdash,vxdash,vhdash,vfdash)

! Interpolation for EXGM

implicit none

real(prec),intent(in)::xdash,hdash
real(prec),intent(out)::cdash,kdash,vxdash,vhdash,vfdash
real(prec)::valsx(2),valsh(2)
integer indsx(2),indsh(2)

real(prec):: a1,a2,b1,b2,c1,c2,d1,d2

!if (.not. marker(xc,hc)) then
!    if (hdash.gt.gridhex(nh)) then
!       print*,'h out of bounds, increase maximum human capital', xc,hc,it
!    end if
!    if (xdash .gt. gridxex(nx)) then
!       print*,'x out of bounds, increase maximum asset', xc,hc,it
!    end if
!    if (hdash .lt. gridhex(1)) then
!       print*,'h out of bounds, increase minimum human capital', xc,hc,it
!    end if
!end if 

call sub_basefun(gridxex,nx,xdash,valsx,indsx)
call sub_basefun(gridhex,nh,hdash,valsh,indsh)

cdash  = valsx(1)*valsh(1)*calt(indsx(1),indsh(1))+valsx(1)*valsh(2)*calt(indsx(1),indsh(2))+valsx(2)*valsh(1)*calt(indsx(2),indsh(1))+valsx(2)*valsh(2)*calt(indsx(2),indsh(2))
kdash  = valsx(1)*valsh(1)*kalt(indsx(1),indsh(1))+valsx(1)*valsh(2)*kalt(indsx(1),indsh(2))+valsx(2)*valsh(1)*kalt(indsx(2),indsh(1))+valsx(2)*valsh(2)*kalt(indsx(2),indsh(2))
vfdash = valsx(1)*valsh(1)*valt(indsx(1),indsh(1))+valsx(1)*valsh(2)*valt(indsx(1),indsh(2))+valsx(2)*valsh(1)*valt(indsx(2),indsh(1))+valsx(2)*valsh(2)*valt(indsx(2),indsh(2))

if (cdash .gt. 0.0) then
    vxdash = func_MUc(cdash)
    vhdash = (netw+1/func_MFk(kdash))*func_MUc(cdash)
else
    a1 = calt(indsx(1),indsh(1))
    b1 = calt(indsx(1),indsh(2))
    c1 = calt(indsx(2),indsh(1))
    d1 = calt(indsx(2),indsh(2))

    a2 = kalt(indsx(1),indsh(1))
    b2 = kalt(indsx(1),indsh(2))
    c2 = kalt(indsx(2),indsh(1))
    d2 = kalt(indsx(2),indsh(2))

    vxdash = valsx(1)*valsh(1)*func_MUc(a1)+valsx(1)*valsh(2)*func_MUc(b1)+valsx(2)*valsh(1)*func_MUc(c1)+valsx(2)*valsh(2)*func_MUc(d1)
    vhdash = valsx(1)*valsh(1)*(netw+1/func_MFk(a2))*func_MUc(a1)+valsx(1)*valsh(2)*(netw+1/func_MFk(b2))*func_MUc(b1)+valsx(2)*valsh(1)*(netw+1/func_MFk(c2))*func_MUc(c1)+valsx(2)*valsh(2)*(netw+1/func_MFk(d2))*func_MUc(d1)     
end if 
 
end subroutine sub_exopolation
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
subroutine sub_hypolation(xdash,hdash,cdash,kdash,vxdash,vhdash,vfdash)

! Interpolation for HEGM

implicit none

real(prec),intent(in)::xdash,hdash
real(prec),intent(out)::cdash,kdash,vxdash,vhdash,vfdash
real(prec)::vxtrans,vhtrans
real(prec)::valsx1(2),valsx2(2),valsh(2)
integer indsx1(2),indsx2(2),indsh(2)

call sub_basefun(gridhhy,nh,hdash,valsh,indsh)

call sub_basefun(gridxhy(:,indsh(1)),nx,xdash,valsx1,indsx1)

!if (xdash .gt. gridxhy(nx,indsh(2))) then
!   print*,'x out of bounds, increase maximum savings', xc,hc,it
!else if (xdash.lt.gridxhy(1,indsh(2))-epsi) then
!   print*,'x out of bounds, increase minimum savings', xc,hc,it
!else
   call sub_basefun(gridxhy(:,indsh(2)),nx,xdash,valsx2,indsx2)
!end if

cdash  = valsx1(1)*valsh(1)*calt(indsx1(1),indsh(1))+valsx1(2)*valsh(1)*calt(indsx1(2),indsh(1))+valsx2(1)*valsh(2)*calt(indsx2(1),indsh(2))+valsx2(2)*valsh(2)*calt(indsx2(2),indsh(2))
kdash  = valsx1(1)*valsh(1)*kalt(indsx1(1),indsh(1))+valsx1(2)*valsh(1)*kalt(indsx1(2),indsh(1))+valsx2(1)*valsh(2)*kalt(indsx2(1),indsh(2))+valsx2(2)*valsh(2)*kalt(indsx2(2),indsh(2))

vxdash = func_MUc(cdash)
vhdash = (netw+(1/gamma)*kdash**alpha)*func_MUc(cdash)

vfdash  = valsx1(1)*valsh(1)*valt(indsx1(1),indsh(1))+valsx1(2)*valsh(1)*valt(indsx1(2),indsh(1))+valsx2(1)*valsh(2)*valt(indsx2(1),indsh(2))+valsx2(2)*valsh(2)*valt(indsx2(2),indsh(2))

end subroutine sub_hypolation
!---------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------
subroutine sub_basefun (grid_x,npx,x,vals,inds) 

! interpolation: basis functions 

! this subroutine returns the values and the indices of the two basis
! functions that are positive on a given x in the grid_x

implicit none

real(prec),intent(in) :: x
integer , intent(in):: npx
real(prec), intent(in) :: grid_x (npx)
real(prec), intent(out) ::vals(2)
integer ,intent(out) ::inds(2)
integer :: i

call sub_lookup(i,x,grid_x,npx)

vals(2) = ( x-grid_x(i-1) )/(grid_x(i)-grid_x(i-1))
vals(1) = ( grid_x(i)-x )/(grid_x(i)-grid_x(i-1))
inds(2) = i
inds(1) = i-1

end subroutine sub_basefun
! --------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------
subroutine sub_lookup(i,x,grid_x,npx)

! lookup position of x in gridx

integer::ju,jl,jm
integer,intent(in)::npx
real(prec),intent(in)::x,grid_x(npx)
integer, intent(out)::i

jl = 1      
ju = npx   

do
	if (ju-jl .le. 1) exit
	jm = (ju+jl)/2
	if (x .ge. grid_x(jm)) then
		jl=jm
	else
		ju=jm
	endif
end do

i = jl+1

end subroutine sub_lookup
! --------------------------------------------------------------------------------------

end module comfunc