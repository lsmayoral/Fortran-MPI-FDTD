	module graphene633_mod
	use constants_mod
	use mydatatypes
	use material_mod
	implicit none
	private
	public :: graphene633

	type, extends(material):: graphene633
	real*8,public :: dt,wNL,dnuNL,xi3,alpha
!**************************************************************
!	Variables to update Electric field
!**************************************************************
	real*8, allocatable, dimension(:) :: ga,gb,In,raman0,raman1,raman2,Sn,Snm1,Snm2
	real*8 :: n,k,sigma,ep_inf,scaleFactor
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	    contains
	    procedure :: update_E => update_E_inSLG
	    procedure :: initialize_material => initialize_SLG
	    procedure :: initialize_SLGVars
	    procedure :: allocate_SLGVars
	end type graphene633

	contains 

subroutine update_E_inSLG(this,E_pointer,D_pointer)
	class(graphene633) :: this
	real*8, pointer, intent(inout):: D_pointer(:)
        real*8,pointer, intent(inout) :: E_pointer(:)
	integer :: pos
	    if(this%posrange%init.gt.0.and.this%posrange%ends.gt.0)then
		do pos=this%posrange%init,this%posrange%ends,1
	E_pointer(pos)=(D_pointer(pos)-this%In(pos)+&
	    2.d0*this%xi3*this%alpha*(E_pointer(pos)**3))/&
	    (this%ga(pos)+this%Sn(pos)+3.d0*this%xi3*this%alpha*(E_pointer(pos)**2))
	this%In(pos)=this%In(pos)+this%gb(pos)*E_pointer(pos)

	this%Sn(pos)=this%raman1(pos)*this%Snm1(pos)+&
	    this%raman2(pos)*this%Snm2(pos)+this%raman0(pos)*(E_pointer(pos)**2)
	this%Snm2(pos)=this%Snm1(pos)
	this%Snm1(pos)=this%Sn(pos)
		enddo
	    endif
	end subroutine update_E_inSLG

	subroutine initialize_SLG(this,absrange,localdim)
	    class(graphene633),intent(inout) :: this
	    type(Irange), intent(in) :: absrange
	    integer, intent(in) :: localdim
	    this%n=2.731108d0; this%k=1.357d0
	    this%sigma=dsqrt(this%scaleFactor)*4.d0*pi*this%n*this%k*c0*ep0/(633.d0*nano)
	    this%ep_inf=this%scaleFactor*((this%n**2)-(this%k**2))
!	    this%ep_inf=((this%n**2)-(this%k**2))
	    call allocate_SLGVars(this,localdim)
	    call this%material%initialize_material_indexes(absrange,localdim)
	    call initialize_SLGVars(this)
        end subroutine initialize_SLG

	subroutine allocate_SLGVars(this,localdim)
	    class(graphene633) :: this
	    integer, intent(in) :: localdim
	    integer :: alloc_stat

	    allocate(this%ga(localdim),stat=alloc_stat)
	    allocate(this%gb(localdim),stat=alloc_stat)
	    allocate(this%In(localdim),stat=alloc_stat)
!***********************************************************
!		Nonlinear variables for one Raman peak
!***********************************************************
	    allocate(this%raman0(localdim),stat=alloc_stat)
	    allocate(this%raman1(localdim),stat=alloc_stat)
	    allocate(this%raman2(localdim),stat=alloc_stat)
	    allocate(this%Sn(localdim),stat=alloc_stat)
	    allocate(this%Snm1(localdim),stat=alloc_stat)
	    allocate(this%Snm2(localdim),stat=alloc_stat)

	    this%raman0=0.d0
	    this%raman1=0.d0
	    this%raman2=0.d0
	    this%Sn=0.d0
	    this%Snm1=0.d0
	    this%Snm2=0.d0

	    this%ga=1.d0
	    this%gb=0.d0
	    this%In=0.d0
	end subroutine allocate_SLGVars

	subroutine initialize_SLGVars(this)
	class(graphene633) :: this
	integer :: pos
	real*8 :: dt,omegaR,betaR
	dt=this%dt
	betaR=dsqrt((this%wNL**2)-(this%dnuNL**2))
	omegaR=(this%wNL**2)/betaR

	do pos=this%posrange%init,this%posrange%ends,1
	    this%ga(pos)=this%ep_inf+this%sigma*dt/ep0
	    this%gb(pos)=this%sigma*dt/ep0

	    this%raman0(pos)=(1.d0-this%alpha)*this%xi3*omegaR*dt*dexp(-this%dnuNL*dt)*dsin(betaR*dt)
	    this%raman1(pos)=2.d0*dexp(-this%dnuNL*dt)*dcos(betaR*dt)
	    this%raman2(pos)=-dexp(-2.d0*this%dnuNL*dt)
	enddo
	end subroutine initialize_SLGVars

	end module graphene633_mod


