module comfuncform

use comtypes
use comparam
use comvar

implicit none

contains

! --------------------------------------------------------------------------------------
function func_U(c)

! Utility Function

implicit none
real(prec) func_U,c
real(prec)::epsi=0.000001

if (abs(tetta-1.0)<epsi) then
	func_U=log(c)
else
	func_U=1.0/(1.0-tetta)*(c**(1.0-tetta))
endif

end function func_U
! --------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------
function func_MUc(c)

! Marginal Utility of Consumption

implicit none
real(prec)::func_MUc,c
real(prec),parameter::penscale=100000000000000.0

if (c>0.0) then
	func_MUc=c**(-tetta)
else
	func_MUc=penscale+abs(100000*c)**3.0
endif

end function func_MUc
! --------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------
function func_invuc(marg)

! Inverse of the Marginal Utility of Consumption

implicit none
real(prec)::func_invuc,marg

if (marg<=0.0) then
	print*,'Trouble Inverting for Consumption ', xc,hc,it
	pause
else
	func_invuc=marg**(-1.0/tetta)
endif
   
end function func_invuc
! --------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------
function func_F(k)

! Health Production Function
implicit none
real(prec) func_F,k
real(prec)::epsi=0.000001

if (abs(alpha-1.0)<epsi) then
	func_F=log(k)
else
	func_F=gamma*1.0/(1.0-alpha)*(k**(1.0-alpha))
endif

end function func_F
! --------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------
function func_MFk(k)

! Marginal Health Production

implicit none
real(prec) func_MFk,k
real(prec),parameter::penscale=1000000000.0

if (k>0.0) then
	func_MFk=(gamma)*k**(-alpha)
else if (k>-epsi) then
    func_MFk=penscale+abs(k)**2.0
else
    print* , 'negative investment in human capital ', xc,hc,it   	
endif

end function func_MFk
! --------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------
function func_invfk(marg)

! Inverse of the Utility Function

implicit none
real(prec)::func_invfk,marg

if (marg<=0.0) then
	print*,'Trouble Inverting for investment ', xc,hc,it
	pause
else
	func_invfk=(marg/gamma)**(-1.0/alpha)
endif
   
end function func_invfk
! --------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------
function func_surv(hdash)

! Inverse of the Utility Function

implicit none
real(prec)::func_surv,hdash

if (hdash<=0.0) then
	print*,'Survival rate of negative health capital', xc,hc,it
	pause
else
	func_surv=1-phi/(1+hdash)
endif
   
end function func_surv
! --------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------
function func_margsurv(hdash)

! Inverse of the Utility Function

implicit none
real(prec)::func_margsurv,hdash

if (hdash<=0.0) then
	print*,'Marginal survival rate of negative health capital', xc,hc,it
	pause
else
	func_margsurv=phi/(1+hdash)**2
endif
   
end function func_margsurv
! --------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------
function func_eulererr_c(c,cprime,hprime)

! Euler Equation Error function

implicit none
real(prec)::func_eulererr_c,c,cprime,hprime

func_eulererr_c = 1.0 - func_invuc( (1+ret)*betta*func_surv(hprime)*func_MUc(cprime) ) / c
   
end function func_eulererr_c
! --------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------
function func_eulererr_k(c,k,hprime,cprime,kprime,vprime)

! Euler Equation Error function

implicit none
real(prec)::func_eulererr_k,c,k,cprime,kprime,hprime,vprime, sss
sss = func_invfk( ((1+ret)/(1-delta))*(func_margsurv(hprime)*vprime/(func_surv(hprime)*func_MUc(cprime))+netw+1/func_MFk(kprime))**(-1) )
func_eulererr_k = 1.0 - sss / k
   
end function func_eulererr_k
! --------------------------------------------------------------------------------------

end module comfuncform