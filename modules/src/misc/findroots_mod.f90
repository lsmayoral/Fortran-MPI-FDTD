	module findRoots_mod
	use constants_mod
	implicit none
	private
	public mullermethod,lookMullerRoot

	type mullermethod
	    real*8, public :: tol,h
	    integer, public :: steps
	    complex*16 :: xr,mullerRoot
	    real*8 :: ampMullerRoot
	    contains 
		procedure :: lookMullerRoot
	end type mullermethod

	contains

	subroutine lookMullerRoot(this,evalfunc,param1) 
	    class(mullermethod) :: this
	    real*8, intent(in),optional :: param1
	    complex*16 :: evalfunc
	    complex*16 :: x0,x1,x2,x3,prevx3
	    complex*16 :: fx0,fx1,fx2,fx3
	    complex*16 :: h0,h1,d0,d1
	    complex*16 :: a,b,c,rad,radp,radm,den
	    real*8 :: ampradp,ampradm,ampfx3,ampx3,amprevx3,p_error
	    integer :: step

!		print *,"h=",this%h
!		print *,"xr=",this%xr

!		x2=this%xr
!		x1=x2-this%h*x2
!		x0=x2-2.d0*this%h*x2

		x2=this%xr
		x1=x2+this%h*x2
		x0=x2-this%h*x2

!		print *,"x0=",x0
!		print *,"x1=",x1
!		print *,"x2=",x2

	    muloop: do step=1,this%steps,1

		h0=x1-x0
		h1=x2-x1

		fx0=evalfunc(x0,param1)
		fx1=evalfunc(x1,param1)
		fx2=evalfunc(x2,param1)

!		print *,"fx0=",fx0,x0,param1
!		print *,"******************"
!		print *,"fx1=",fx1,x1,param1
!		print *,"******************"
!		print *,"fx2",fx2,x2,param1

		d0=(fx1-fx0)/h0
		d1=(fx2-fx1)/h1

		a=(d1-d0)/(h1+h0)
		b=a*h1+d1
		c=fx2

		rad=cdsqrt((b**2)-4.d0*a*c)
		radp=b+rad
		radm=b-rad
		ampradp=dsqrt((real(radp)**2)+(imag(radp)**2))
		ampradm=dsqrt((real(rad)**2)+(imag(radm)**2))
		if(ampradp.gt.ampradm)then
		    den=radp
		else
		    den=radm
		endif
		x3=x2-(2.d0*c/den)
		ampx3=dsqrt((real(x3)**2)+(imag(x3)**2))
		fx3=evalfunc(x3,param1)
		ampfx3=dsqrt((real(fx3)**2)+(imag(fx3)**2))
		if(step.ge.1)then
		    p_error=100.d0*abs(amprevx3-ampx3)/ampx3
		endif

!		print *,"i=",step,"fx3=",ampfx3,"error=",p_error
		if(ampfx3.le.this%tol.or.p_error.eq.0)then
		    exit muloop
		endif
		prevx3=x3
		amprevx3=dsqrt((real(x3)**2)+(imag(x3)**2))
		x0=x1
		x1=x2
		x2=x3
	    enddo muloop
	    this%mullerRoot=x3
	    this%ampMullerRoot=ampfx3
	end subroutine lookMullerRoot



	end module findRoots_mod
