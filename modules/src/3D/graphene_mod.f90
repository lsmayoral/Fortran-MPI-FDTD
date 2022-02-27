	module graphene_mod
	use constants_mod
	use mydatatypes
	use material_mod
	implicit none
	private
	public :: Graphene,clearGrapheneVars

	type, extends(material):: Graphene
	real*8,public :: wNL
	real*8,public :: dnuNL
	real*8,public :: xi3
	real*8,public :: dt
	real*8,public :: ds
	real*8 :: scaleFact=1.d0
!**************************************************************
!	Variables to update Electric field
!**************************************************************
	complex*16, allocatable,dimension(:,:,:) :: gax
	real*8, allocatable,dimension(:,:,:) ::ha,hb,Snx,Inx,Snz,Inz
	complex*16,allocatable,dimension(:,:,:) ::gb1,gb2,gb3,gb4
	complex*16,allocatable,dimension(:,:,:) ::gc1,gc2,gc3,gc4
	complex*16,allocatable,dimension(:,:,:) :: g1x,g2x,g3x,g4x
	complex*16,allocatable,dimension(:,:,:) :: g1z,g2z,g3z,g4z
	real*8, allocatable,dimension(:,:,:) :: SLG_raman0,SLG_raman1,SLG_raman2
	real*8, allocatable,dimension(:,:,:) :: SLG_Snx,SLG_Snxm1,SLG_Snxm2
	real*8, allocatable,dimension(:,:,:) :: SLG_Snz,SLG_Snzm1,SLG_Snzm2
!**************************************************************
	real*8 :: scaled_graphene_wD
!**************************************************************
!		Parameter of 4 Critical points model for Graphene
!**************************************************************
	real*8 :: Graphenewidth=3.4d-10
	real*8 :: graphene_epinf=1.d0
	real*8 :: graphene_dnur=9d15
	real*8 :: graphene_wD=5.2796d15
	real*8,dimension(4)::scaled_ap
	real*8,dimension(4) :: ap= (/ 13.1068d0,1.8417d0,0.4191d0,0.7737d0 /)
	real*8,dimension(4) :: phip=(/ -1.8391d0,0.3505d0,3.0565d0,-0.1199d0 /)
	real*8,dimension(4) :: omp=(/ 6.1264d14,6.6843d15,5.6289d15,1.118d16 /)
	real*8,dimension(4) :: rp=(/ 1d14,7.905d14,6.7526d14,1.0005d14 /)
	complex*16 :: expi1,expi2,expi3,expi4
	complex*16 :: exp1,exp2,exp3,exp4

	    contains
	    procedure :: update_E => update_Exz_inGraphene
	    procedure :: initialize_material => initialize_Graphene
	    procedure :: initialize_GrapheneVars
	    procedure :: allocate_GrapheneVars
	    procedure ::scaleParams
	    procedure ::setexp
	    procedure ::clearGrapheneVars
	    procedure ::getRgax
	    procedure ::getIgax
	end type Graphene

	contains 

subroutine update_Exz_inGraphene(this,Ex_pointer,Ey_pointer,Ez_pointer,Dx_pointer,Dy_pointer,Dz_pointer)
	    class(Graphene) :: this
        real*8,pointer, intent(inout) :: Ex_pointer(:,:,:),Ey_pointer(:,:,:),Ez_pointer(:,:,:)
        real*8,pointer, intent(inout) :: Dx_pointer(:,:,:),Dy_pointer(:,:,:),Dz_pointer(:,:,:)
	    integer :: i,j,k
	    if(this%xrange%init.gt.0.and.this%xrange%ends.gt.0.and.&
	       this%yrange%init.gt.0.and.this%yrange%ends.gt.0.and.&
	       this%zrange%init.gt.0.and.this%zrange%ends.gt.0)then
		do k=this%zrange%init,this%zrange%ends,1
		    do j=this%yrange%init,this%yrange%ends,1
			do i=this%xrange%init,this%xrange%ends,1
Ex_pointer(i,j,k)=(Dx_pointer(i,j,k)-this%Inx(i,j,k)+this%ha(i,j,k)*this%Snx(i,j,k)+this%gb1(i,j,k)*this%g1x(i,j,k)+&
	this%gb2(i,j,k)*this%g2x(i,j,k)+this%gb3(i,j,k)*this%g3x(i,j,k)+this%gb4(i,j,k)*this%g4x(i,j,k))/&
	(this%gax(i,j,k)+this%SLG_Snx(i,j,k))

    this%SLG_Snx(i,j,k)=this%SLG_raman1(i,j,k)*this%SLG_Snxm1(i,j,k)+&
    this%SLG_raman2(i,j,k)*this%SLG_Snxm2(i,j,k)+this%SLG_raman0(i,j,k)*(Ex_pointer(i,j,k)**2)
    this%SLG_Snxm2(i,j,k)=this%SLG_Snxm1(i,j,k)
    this%SLG_Snxm1(i,j,k)=this%SLG_Snx(i,j,k)

	this%Inx(i,j,k)=this%Inx(i,j,k)+ this%hb(i,j,k)*Ex_pointer(i,j,k)
	this%Snx(i,j,k)= this%ha(i,j,k)*this%Snx(i,j,k)+ this%hb(i,j,k)*Ex_pointer(i,j,k)
	this%g1x(i,j,k)=this%gb1(i,j,k)*this%g1x(i,j,k)+this%gc1(i,j,k)*Ex_pointer(i,j,k)
	this%g2x(i,j,k)=this%gb2(i,j,k)*this%g2x(i,j,k)+this%gc2(i,j,k)*Ex_pointer(i,j,k)
	this%g3x(i,j,k)=this%gb3(i,j,k)*this%g3x(i,j,k)+this%gc3(i,j,k)*Ex_pointer(i,j,k)
	this%g4x(i,j,k)=this%gb4(i,j,k)*this%g4x(i,j,k)+this%gc4(i,j,k)*Ex_pointer(i,j,k)

Ez_pointer(i,j,k)=(Dz_pointer(i,j,k)-this%Inz(i,j,k)+this%ha(i,j,k)*this%Snz(i,j,k)+this%gb1(i,j,k)*this%g1z(i,j,k)+&
	this%gb2(i,j,k)*this%g2z(i,j,k)+this%gb3(i,j,k)*this%g3z(i,j,k)+this%gb4(i,j,k)*this%g4z(i,j,k))/&
	(this%gax(i,j,k)+this%SLG_Snz(i,j,k))

    this%SLG_Snz(i,j,k)=this%SLG_raman1(i,j,k)*this%SLG_Snzm1(i,j,k)+&
    this%SLG_raman2(i,j,k)*this%SLG_Snzm2(i,j,k)+this%SLG_raman0(i,j,k)*(Ez_pointer(i,j,k)**2)
    this%SLG_Snzm2(i,j,k)=this%SLG_Snzm1(i,j,k)
    this%SLG_Snzm1(i,j,k)=this%SLG_Snz(i,j,k)

	this%Inz(i,j,k)=this%Inz(i,j,k)+ this%hb(i,j,k)*Ez_pointer(i,j,k)
	this%Snz(i,j,k)= this%ha(i,j,k)*this%Snz(i,j,k)+ this%hb(i,j,k)*Ez_pointer(i,j,k)
	this%g1z(i,j,k)=this%gb1(i,j,k)*this%g1z(i,j,k)+this%gc1(i,j,k)*Ez_pointer(i,j,k)
	this%g2z(i,j,k)=this%gb2(i,j,k)*this%g2z(i,j,k)+this%gc2(i,j,k)*Ez_pointer(i,j,k)
	this%g3z(i,j,k)=this%gb3(i,j,k)*this%g3z(i,j,k)+this%gc3(i,j,k)*Ez_pointer(i,j,k)
	this%g4z(i,j,k)=this%gb4(i,j,k)*this%g4z(i,j,k)+this%gc4(i,j,k)*Ez_pointer(i,j,k)
		    enddo
		enddo
	    enddo
	endif
	end subroutine update_Exz_inGraphene

	subroutine initialize_Graphene(this,rangex,rangey,rangez,range_localdim)
	    class(Graphene),intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey,rangez
	    integer, dimension(3), intent(in) :: range_localdim
	    call allocate_GrapheneVars(this,range_localdim)
	    call this%material%initialize_material_indexes(rangex,rangey,rangez,range_localdim)

	    call initialize_GrapheneVars(this)
        end subroutine initialize_Graphene

	subroutine allocate_GrapheneVars(this,range_localdim)
	    class(Graphene) :: this
	    integer,dimension(3),intent(in) :: range_localdim
	    integer :: alloc_stat,dimx,dimy,dimz
	    dimx=range_localdim(1)
	    dimy=range_localdim(2)
	    dimz=range_localdim(3)
	    allocate(this%gax(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%ha(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%hb(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gb1(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gb2(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gb3(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gb4(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gc1(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gc2(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gc3(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gc4(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Snx(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Inx(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%g1x(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%g2x(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%g3x(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%g4x(dimx,dimy,dimz),stat=alloc_stat)

	    allocate(this%Snz(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Inz(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%g1z(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%g2z(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%g3z(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%g4z(dimx,dimy,dimz),stat=alloc_stat)
!***********************************************************
!		Nonlinear variables for one Raman peak
!***********************************************************
	    allocate(this%SLG_raman0(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%SLG_raman1(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%SLG_raman2(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%SLG_Snx(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%SLG_Snxm1(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%SLG_Snxm2(dimx,dimy,dimz),stat=alloc_stat)

	    allocate(this%SLG_Snz(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%SLG_Snzm1(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%SLG_Snzm2(dimx,dimy,dimz),stat=alloc_stat)
!***********************************************************
	    this%SLG_raman0=0.d0
	    this%SLG_raman1=0.d0
	    this%SLG_raman2=0.d0

	    this%SLG_Snx=0.d0
	    this%SLG_Snxm1=0.d0
	    this%SLG_Snxm2=0.d0

	    this%SLG_Snz=0.d0
	    this%SLG_Snzm1=0.d0
	    this%SLG_Snzm2=0.d0
!***********************************************************
	    this%ha=0.d0
	    this%hb=0.d0
	    this%gb1=zero
	    this%gb2=zero
	    this%gb3=zero
	    this%gb4=zero
	    this%gc1=zero
	    this%gc2=zero
	    this%gc3=zero
	    this%gc4=zero

	    this%Snx=0.d0
	    this%Inx=0.d0
	    this%g1x=zero
	    this%g2x=zero
	    this%g3x=zero
	    this%g4x=zero

	    this%Snz=0.d0
	    this%Inz=0.d0
	    this%g1z=zero
	    this%g2z=zero
	    this%g3z=zero
	    this%g4z=zero

	    this%gax=1.d0
	end subroutine allocate_GrapheneVars

	subroutine initialize_GrapheneVars(this)
	class(Graphene) :: this
	integer :: x,y,z
	type(Irange) :: posx,posy,posz
	real*8 :: xi3,wNL,dnuNL,omegaR,alfaR,betaR,dt,graphene_dnur,graphene_epinf
	complex*16 :: exp1,exp2,exp3,exp4,expi1,expi2,expi3,expi4

	posx=this%xrange
	posy=this%yrange
	posz=this%zrange
	xi3=this%xi3
	wNL=this%wNL
	dnuNL=this%dnuNL

	alfaR=0.5d0*dnuNL
	betaR=dsqrt((wNL**2)-(alfaR**2))
	omegaR=(wNL**2)/betaR

	dt=this%dt
	graphene_dnur=this%graphene_dnur
	exp1=this%exp1
	exp2=this%exp2
	exp3=this%exp3
	exp4=this%exp4
	expi1=this%expi1
	expi2=this%expi2
	expi3=this%expi3
	expi4=this%expi4
	graphene_epinf=this%graphene_epinf
	call setexp(this)
	call this%scaleParams()
	if(this%scaleFact.eq.1.d0)then
	    this%scaled_ap=this%ap
	    this%scaled_graphene_wD=this%graphene_wD
	endif

	    do z=posz%init,posz%ends,1
		do y=posy%init,posy%ends,1
		    do x=posx%init,posx%ends,1
		this%ha(x,y,z)=dexp(-graphene_dnur*dt)
		this%hb(x,y,z)=(this%scaled_graphene_wD**2)*dt/graphene_dnur
		this%gb1(x,y,z)=exp1
		this%gb2(x,y,z)=exp2
		this%gb3(x,y,z)=exp3
		this%gb4(x,y,z)=exp4
		this%gc1(x,y,z)=dt*2.d0*unoi*this%scaled_ap(1)*this%omp(1)*expi1
		this%gc2(x,y,z)=dt*2.d0*unoi*this%scaled_ap(2)*this%omp(2)*expi2
		this%gc3(x,y,z)=dt*2.d0*unoi*this%scaled_ap(3)*this%omp(3)*expi3
		this%gc4(x,y,z)=dt*2.d0*unoi*this%scaled_ap(4)*this%omp(4)*expi4
		this%gax(x,y,z)=graphene_epinf-this%gc1(x,y,z)-this%gc2(x,y,z)-this%gc3(x,y,z)-this%gc4(x,y,z)

!		this%SLG_raman0(x,y,z)=xi3*((wNL*dt)**2)/(1.d0+dnuNL*dt)
!		this%SLG_raman1(x,y,z)=(2.d0-(wNL*dt)**2)/(1.d0+dnuNL*dt)
!		this%SLG_raman2(x,y,z)=-(1.d0-dnuNL*dt)/(1.d0+dnuNL*dt)
		this%SLG_raman0(x,y,z)=xi3*omegaR*dt*dexp(-alfaR*dt)*dsin(betaR*dt)
		this%SLG_raman1(x,y,z)=2.d0*dexp(-alfaR*dt)*dcos(betaR*dt)
		this%SLG_raman2(x,y,z)=-dexp(-2.d0*alfaR*dt)
		enddo
	    enddo
	enddo
	end subroutine initialize_GrapheneVars

	subroutine scaleParams(this)
	    class(Graphene) :: this
	    this%scaleFact=this%Graphenewidth/this%ds
	    this%scaled_ap=this%scaleFact*this%ap
	    this%scaled_graphene_wD=dsqrt(this%scaleFact)*this%graphene_wD
	end subroutine scaleParams

	subroutine setexp(this)
	    class(Graphene) :: this

	    this%expi1=cdexp(-unoi*this%phip(1))
	    this%expi2=cdexp(-unoi*this%phip(2))
	    this%expi3=cdexp(-unoi*this%phip(3))
	    this%expi4=cdexp(-unoi*this%phip(4))

	    this%exp1=cdexp((-this%rp(1)+unoi*this%omp(1))*this%dt)
	    this%exp2=cdexp((-this%rp(2)+unoi*this%omp(2))*this%dt)
	    this%exp3=cdexp((-this%rp(3)+unoi*this%omp(3))*this%dt)
	    this%exp4=cdexp((-this%rp(4)+unoi*this%omp(4))*this%dt)
	end subroutine setexp

	subroutine clearGrapheneVars(this)
	    class(Graphene) :: this
	    this%ha=0.d0
	    this%hb=0.d0
	    this%Snx=0.d0
	    this%Inx=0.d0
	    this%gb1=zero
	    this%gb2=zero
	    this%gb3=zero
	    this%gb4=zero
	    this%gc1=zero
	    this%gc2=zero
	    this%gc3=zero
	    this%gc4=zero

	    this%g1x=zero
	    this%g2x=zero
	    this%g3x=zero
	    this%g4x=zero

	    this%g1z=zero
	    this%g2z=zero
	    this%g3z=zero
	    this%g4z=zero

	    this%SLG_raman0=0.d0
	    this%SLG_raman1=0.d0
	    this%SLG_raman2=0.d0

	    this%SLG_Snx=0.d0
	    this%SLG_Snxm1=0.d0
	    this%SLG_Snxm2=0.d0

	    this%SLG_Snz=0.d0
	    this%SLG_Snzm1=0.d0
	    this%SLG_Snzm2=0.d0
	end subroutine clearGrapheneVars


    function getRgax(this,posx,posy,posz) result(val)
        class(Graphene) :: this
        integer, intent(in):: posx,posy,posz
        real*8 :: val
        val=real(this%gax(posx,posy,posz))
    end function getRgax

    function getIgax(this,posx,posy,posz) result(val)
        class(Graphene) :: this
        integer, intent(in):: posx,posy,posz
        real*8 :: val
        val=imag(this%gax(posx,posy,posz))
    end function getIgax
	end module
