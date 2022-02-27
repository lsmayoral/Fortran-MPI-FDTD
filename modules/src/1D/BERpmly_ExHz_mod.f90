    module BERpmly_ExHz_mod
    use constants_mod
    implicit none
    private
    public pmly_ExHz,initialize,updateE_pmy,updateH_pmy
	type :: pmly_ExHz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	PML public parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    integer :: npmly, ycoord,N_y,Ymax,yBDlow,yBDhigh
!	    real*8,  public :: order
!	    real*8,  public :: R0
!	    real*8,  public :: dy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	PML private variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    real*8, allocatable :: sigmayE(:)
	    real*8, allocatable :: alphayE(:)
	    real*8, allocatable :: f1(:)
	    real*8, allocatable :: f2(:)

	    real*8, allocatable :: sigmayH(:)
	    real*8, allocatable :: alphayH(:)
	    real*8, allocatable :: g1(:)
	    real*8, allocatable :: g2(:)

	    real*8, allocatable, public :: Ex(:),Hz(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	PML procedures (subroutines or functions)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    contains
	    procedure :: updateE => updateE_pmy
	    procedure :: setboundaryEa
	    procedure :: setboundaryEb
	    procedure :: updateE_pmy_low
	    procedure :: updateE_pmy_high
	    procedure :: updateH => updateH_pmy
	    procedure :: setboundaryHa
	    procedure :: setboundaryHb
	    procedure :: updateH_pmy_low
	    procedure :: updateH_pmy_high
	    procedure :: initialize
	    procedure :: initialize_vars
	    procedure :: allocate_vars
	end type pmly_ExHz

    contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	update fields calls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine updateE_pmy(this,Ex_)
	    class(pmly_ExHz) :: this
	    real*8, pointer, intent(inout) :: Ex_(:)
!	    real*8, pointer, intent(inout) :: Hz_(:)
	    if(this%ycoord.eq.0)then
		call this%updateE_pmy_low()
		call this%setboundaryEa(Ex_)
	    endif
	    if(this%ycoord.eq.(this%Ymax-1))then
		call this%updateE_pmy_high()
		call this%setboundaryEb(Ex_)
	    endif
	end subroutine updateE_pmy
	subroutine updateH_pmy(this,Hz_)
	    class(pmly_ExHz) :: this
!	    real*8, pointer, intent(inout) :: Ex_(:)
	    real*8, pointer, intent(inout) :: Hz_(:)
	    if(this%ycoord.eq.0)then
		call this%updateH_pmy_low()
		call this%setboundaryHa(Hz_)
	    endif
	    if(this%ycoord.eq.(this%Ymax-1))then
		call this%updateH_pmy_high()
		call this%setboundaryHb(Hz_)
	    endif
	end subroutine updateH_pmy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	update equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ex field
	subroutine updateE_pmy_low(this)
        class(pmly_ExHz) :: this
!	real*8, pointer, intent(inout) :: Ex_(:)
!        real*8, pointer, intent(inout) :: Hz_(:)
	real*8 :: curlH
	integer :: j
	do j=1,this%yBDlow,1
	    curlH=(this%Hz(j)-this%Hz(j-1))
	    this%Ex(j)=this%f1(j)*this%Ex(j)+this%f2(j)*curlH
	enddo
	end subroutine updateE_pmy_low

	subroutine updateE_pmy_high(this)
        class(pmly_ExHz) :: this
!	real*8, pointer, intent(inout) :: Ex_(:)
!        real*8, pointer, intent(inout) :: Hz_(:)
	real*8 :: curlH
	integer :: j
	do j=this%yBDhigh,this%N_y,1
	    curlH=(this%Hz(j)-this%Hz(j-1))
	    this%Ex(j)=this%f1(j)*this%Ex(j)+this%f2(j)*curlH
	enddo
	end subroutine updateE_pmy_high

! Hz field
	subroutine updateH_pmy_low(this)
        class(pmly_ExHz) :: this
!	real*8, pointer, intent(inout) :: Ex_(:)
!        real*8, pointer, intent(inout) :: Hz_(:)
	real*8 :: curlE
	integer :: j
	do j=1,this%yBDlow,1
	    curlE=(this%Ex(j+1)-this%Ex(j))
	    this%Hz(j)=this%g1(j)*this%Hz(j)+this%g2(j)*curlE
	enddo
	end subroutine updateH_pmy_low

	subroutine updateH_pmy_high(this)
        class(pmly_ExHz) :: this
!	real*8, pointer, intent(inout) :: Ex_(:)
!        real*8, pointer, intent(inout) :: Hz_(:)
	real*8 :: curlE
	integer :: j
	do j=this%yBDhigh,this%N_y,1
	    curlE=(this%Ex(j+1)-this%Ex(j))
	    this%Hz(j)=this%g1(j)*this%Hz(j)+this%g2(j)*curlE
	enddo
	end subroutine updateH_pmy_high

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine setboundaryEa(this,Ex_)
	class(pmly_ExHz) :: this
	real*8, pointer, intent(inout) :: Ex_(:)
	this%Ex(this%yBDlow+1)=Ex_(this%yBDlow+1)
	Ex_(this%yBDlow)=this%Ex(this%yBDlow)
	end subroutine setboundaryEa

	subroutine setboundaryHa(this,Hz_)
	class(pmly_ExHz) :: this
	real*8, pointer, intent(inout) :: Hz_(:)
	this%Hz(this%yBDlow+1)=Hz_(this%yBDlow+1)
	Hz_(this%yBDlow)=this%Hz(this%yBDlow)
	end subroutine setboundaryHa

	subroutine setboundaryEb(this,Ex_)
	class(pmly_ExHz) :: this
	real*8, pointer, intent(inout) :: Ex_(:)
	this%Ex(this%yBDhigh-1)=Ex_(this%yBDhigh-1)
	Ex_(this%yBDhigh)=this%Ex(this%yBDhigh)
	end subroutine setboundaryEb

	subroutine setboundaryHb(this,Hz_)
	class(pmly_ExHz) :: this
	real*8, pointer, intent(inout) :: Hz_(:)
	this%Hz(this%yBDhigh-1)=Hz_(this%yBDhigh-1)
	Hz_(this%yBDhigh)=this%Hz(this%yBDhigh)
	end subroutine setboundaryHb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine initialize(this,eplow,ephigh,sigmam,npmly,order,dt,ds,ycoord,range,Ymax)
	    class(pmly_ExHz) :: this
	    integer, intent(in) :: npmly,ycoord,range,Ymax
	    real*8, intent(in) :: eplow,ephigh,order,dt,ds,sigmam
	    this%npmly=npmly
	    this%ycoord=ycoord
	    this%N_y=range
	    this%Ymax=Ymax
	    this%yBDlow=this%npmly
	    this%yBDhigh=this%N_y-this%npmly+1
	    call this%allocate_vars(eplow/ep0,ephigh/ep0,dt*c0/ds)
	    call this%initialize_vars(eplow,ephigh,sigmam,order,dt,ds)
	end subroutine initialize

	subroutine allocate_vars(this,eprlow,eprhigh,Ku)
	    class(pmly_ExHz) :: this
	    real*8, intent(in) :: eprlow,eprhigh,Ku
	    real*8 :: epr
	    integer :: N_y
	    N_y=this%N_y

	    if(this%ycoord.eq.0.or.this%ycoord.eq.(this%Ymax-1))then
		allocate(this%Ex(0:N_y+1))
		allocate(this%sigmayE(N_y))
		allocate(this%alphayE(N_y))
		allocate(this%f1(N_y))
		allocate(this%f2(N_y))

		allocate(this%Hz(0:N_y+1))
		allocate(this%sigmayH(N_y))
		allocate(this%alphayH(N_y))
		allocate(this%g1(N_y))
		allocate(this%g2(N_y))

		if(this%ycoord.eq.0)then
		    this%f2=Ku/eprlow
		endif

		if(this%ycoord.eq.(this%Ymax-1))then
		    this%f2=Ku/eprhigh
		endif

	        this%sigmayE=0.d0; this%alphayE=0.d0; this%f1=1.d0
		this%sigmayH=0.d0; this%alphayH=0.d0; this%g1=1.d0
		this%g2=Ku
		this%Ex=0.d0
		this%Hz=0.d0
	    endif
	end subroutine allocate_vars

	subroutine initialize_vars(this,eplow,ephigh,sigmam,order,dt,ds)
	    class(pmly_ExHz) :: this
	    integer :: npml
	    real*8, intent(in) :: eplow,ephigh,order,dt,ds,sigmam
	    integer :: j,jj, lastY,npmly,ny,yBDlow,yBDhigh
	    real*8 :: z0
	    lastY=this%Ymax-1
	    yBDlow=this%yBDlow; yBDhigh=this%yBDhigh
	    npmly=this%npmly
	    ny=this%N_y
	    z0=dsqrt(mu0/ep0)
	    if(this%ycoord.eq.0.and.yBDlow.gt.0)then
		do j=1,yBDlow,1
		    this%sigmayE(j)=eplow*sigmam*((real(npmly-j,8)+0.5d0)/real(npmly,8))**order
		    this%alphayE(j)=this%sigmayE(j)*dt/eplow
			 this%f1(j)=dexp(-this%alphayE(j))
			 this%f2(j)=(1.d0-this%f1(j))/(z0*this%sigmayE(j)*ds)
		    if(j.lt.npmly)then
			this%sigmayH(j)=mu0*sigmam*((real(npmly-j,8))/real(npmly,8))**order
			this%alphayH(j)=this%sigmayH(j)*dt/mu0
			     this%g1(j)=dexp(-this%alphayH(j))
			    this%g2(j)=z0*(1.d0-this%g1(j))/(this%sigmayH(j)*ds)
		    endif
		enddo
	    endif
	    if(this%ycoord.eq.lastY.and.yBDlow.gt.0)then
		jj=npmly
		do j=yBDhigh,ny,1
		    if(jj.lt.npmly)then
			this%sigmayE(j)=ephigh*sigmam*((real(npmly-jj,8))/real(npmly,8))**order
			this%alphayE(j)=this%sigmayE(j)*dt/ephigh
			 this%f1(j)=dexp(-this%alphayE(j))
			 this%f2(j)=(1.d0-this%f1(j))/(z0*this%sigmayE(j)*ds)
		    endif
		    this%sigmayH(j)=mu0*sigmam*((real(npmly-jj,8)+0.5d0)/real(npmly,8))**order
		    this%alphayH(j)=this%sigmayH(j)*dt/mu0
			 this%g1(j)=dexp(-this%alphayH(j))
			 this%g2(j)=z0*(1.d0-this%g1(j))/(this%sigmayH(j)*ds)
		    jj=jj-1
		enddo
	    endif
	end subroutine initialize_vars

    end module
