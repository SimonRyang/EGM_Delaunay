module sing_dim

use comtypes
use comvar

implicit none

contains
!-----------------------------------------------------------
SUBROUTINE zbrac(func,x1,x2,succes)
use comtypes
IMPLICIT NONE
REAL(prec), INTENT(INOUT) :: x1,x2
LOGICAL, INTENT(OUT) :: succes
INTERFACE
	FUNCTION func(x)
	use comtypes
	IMPLICIT NONE
	REAL(prec), INTENT(IN) :: x
	REAL(prec) :: func
	END FUNCTION func
END INTERFACE
INTEGER, PARAMETER :: NTRY=100
REAL(prec), PARAMETER :: FACTOR=1.0_prec
INTEGER :: j
REAL(prec) :: f1,f2

if (x1 == x2) call nrerror('zbrac: you have to guess an initial range')
f1=func(x1)

f2=func(x2)


succes=.true.
do j=1,NTRY
	if ((f1 > 0.0 .and. f2 < 0.0) .or. &
		(f1 < 0.0 .and. f2 > 0.0)) RETURN
	if (abs(f1) < abs(f2)) then
		x1=x1+FACTOR*(x1-x2)
		f1=func(x1)
	else
		x2=x2+FACTOR*(x2-x1)
		f2=func(x2)
	end if
end do

succes=.false.
END SUBROUTINE zbrac

!-----------------------------------------------------------
FUNCTION zbrent(func,x1,x2,tol)

use comtypes
IMPLICIT NONE
REAL(prec), INTENT(IN) :: x1,x2,tol
REAL(prec) :: zbrent
INTERFACE
	FUNCTION func(x)
    use comtypes
	IMPLICIT NONE
	REAL(prec), INTENT(IN) :: x
	REAL(prec) :: func
	END FUNCTION func
END INTERFACE
INTEGER, PARAMETER :: ITMAX=1000
REAL(prec), PARAMETER :: EPS=epsilon(x1)
INTEGER :: iter
REAL(prec) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm

a=x1
b=x2
fa=func(a)
fb=func(b)
if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
    print*, 'zbrent: root must be bracketed for zbrent'
!	call nrerror('root must be bracketed for zbrent')
c=b
fc=fb
do iter=1,ITMAX
	if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
		c=a
		fc=fa
		d=b-a
		e=d
	end if
	if (abs(fc) < abs(fb)) then
		a=b
		b=c
		c=a
		fa=fb
		fb=fc
		fc=fa
	end if
	tol1=1.0_prec*EPS*abs(b)+0.5_prec*tol
	xm=0.5_prec*(c-b)
	if (abs(xm) <= tol1 .or. fb == 0.0) then ! this shows convergence
		zbrent=b
		RETURN
	end if
	if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
		s=fb/fa
		if (a == c) then
			p=2.0_prec*xm*s
			q=1.0_prec-s
		else
			q=fa/fc
			r=fb/fc
			p=s*(2.0_prec*xm*q*(q-r)-(b-a)*(r-1.0_prec))
			q=(q-1.0_prec)*(r-1.0_prec)*(s-1.0_prec)
		end if
		if (p > 0.0) q=-q
		p=abs(p)
		if (2.0_prec*p  <  min(3.0_prec*xm*q-abs(tol1*q),abs(e*q))) then
			e=d
			d=p/q
		else
			d=xm
			e=d
		end if
	else
		d=xm
		e=d
	end if
	a=b
	fa=fb
	b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
	fb=func(b)
end do
 print*, 'zbrent: exceeded maximum iterations'
!	call nrerror('zbrent: exceeded maximum iterations')
zbrent=b
END FUNCTION zbrent

!-------------------------------------------
SUBROUTINE nrerror(string)
CHARACTER(LEN=*), INTENT(IN) :: string
write (*,*) 'nrerror: ',string
STOP 'program terminated by nrerror'
END SUBROUTINE nrerror
!------------------------------------------


end module sing_dim