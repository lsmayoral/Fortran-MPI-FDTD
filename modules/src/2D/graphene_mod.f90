	module graphene_mod
	use constants_mod
	use mydatatypes
	use material_mod
	implicit none
	private
	public :: Graphene,clearGrapheneVars,getRgax,getIgax

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
	complex*16, allocatable,dimension(:,:) :: gax
	real*8, allocatable,dimension(:,:) ::ha,hb,Snx,Inx
	complex*16,allocatable,dimension(:,:) ::gb1,gb2,gb3,gb4
	complex*16,allocatable,dimension(:,:) ::gc1,gc2,gc3,gc4
	complex*16,allocatable,dimension(:,:) :: g1x,g2x,g3x,g4x
	real*8, allocatable,dimension(:,:) :: SLG_raman0,SLG_raman1,SLG_raman2,SLG_Snx,SLG_Snxm1,SLG_Snxm2
!**************************************************************
	real*8 :: scaled_graphene_wD
!**************************************************************
!		Parameter of 4 Critical points model for Graphene
!**************************************************************
	real*8 :: Graphenewidth=3.4d-10
	real*8 :: graphene_epinf=1.d0
	real*8 :: graphene_dnur=9d15
	real*8 :: graphene_wD=5.2796d15
	real*8,dimension(4) :: scaled_ap
	real*8,dimension(4) :: ap= (/ 13.1068d0,1.8417d0,0.4191d0,0.7737d0 /)
	real*8,dimension(4) :: phip=(/ -1.8391d0,0.3505d0,3.0565d0,-0.1199d0 /)
	real*8,dimension(4) :: omp=(/ 6.1264d14,6.6843d15,5.6289d15,1.118d16 /)
	real*8,dimension(4) :: rp=(/ 1d14,7.905d14,6.7526d14,1.0005d14 /)
	complex*16 :: expi1,expi2,expi3,expi4
	complex*16 :: exp1,exp2,exp3,exp4

	    contains
	    procedure :: update_E => update_Ex_inGraphene
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

subroutine update_Ex_inGraphene(this,Ex_pointer,Ey_pointer,Dx_pointer,Dy_pointer)
	    class(Graphene) :: this
	real*8, pointer, intent(inout) ::Ey_pointer(:,:),Dy_pointer(:,:)
        real*8,pointer, intent(inout) :: Ex_pointer(:,:)
        real*8,pointer, intent(inout) :: Dx_pointer(:,:)
	    integer :: i,j,ic,jc
	    if(this%xrange%init.gt.0.and.this%xrange%ends.gt.0.and.&
		this%yrange%init.gt.0.and.this%yrange%ends.gt.0)then
		ic=(this%xrange%init+this%xrange%ends)*0.5
		jc=(this%yrange%init+this%yrange%ends)*0.5
		do i=this%xrange%init,this%xrange%ends,1
		    do j=this%yrange%init,this%yrange%ends,1
Ey_pointer(i,j)=Dy_pointer(i,j)

Ex_pointer(i,j)=(Dx_pointer(i,j)-this%Inx(i,j)+this%ha(i,j)*this%Snx(i,j)+this%gb1(i,j)*this%g1x(i,j)+&
	this%gb2(i,j)*this%g2x(i,j)+this%gb3(i,j)*this%g3x(i,j)+this%gb4(i,j)*this%g4x(i,j))/&
	(this%gax(i,j)+this%SLG_Snx(i,j))

    this%SLG_Snx(i,j)=this%SLG_raman1(i,j)*this%SLG_Snxm1(i,j)+&
    this%SLG_raman2(i,j)*this%SLG_Snxm2(i,j)+this%SLG_raman0(i,j)*(Ex_pointer(i,j)**2)
    this%SLG_Snxm2(i,j)=this%SLG_Snxm1(i,j)
    this%SLG_Snxm1(i,j)=this%SLG_Snx(i,j)

	if(i.eq.ic.and.j.eq.jc)then
	    write(9,*) this%SLG_Snx(i,j)
	endif

	this%Inx(i,j)=this%Inx(i,j)+this%hb(i,j)*Ex_pointer(i,j)
	this%Snx(i,j)=this%ha(i,j)*this%Snx(i,j)+this%hb(i,j)*Ex_pointer(i,j)
	this%g1x(i,j)=this%gb1(i,j)*this%g1x(i,j)+this%gc1(i,j)*Ex_pointer(i,j)
	this%g2x(i,j)=this%gb2(i,j)*this%g2x(i,j)+this%gc2(i,j)*Ex_pointer(i,j)
	this%g3x(i,j)=this%gb3(i,j)*this%g3x(i,j)+this%gc3(i,j)*Ex_pointer(i,j)
	this%g4x(i,j)=this%gb4(i,j)*this%g4x(i,j)+this%gc4(i,j)*Ex_pointer(i,j)
		    enddo
		enddo
	    endif
	end subroutine update_Ex_inGraphene

	subroutine initialize_Graphene(this,rangex,rangey,range_localdim)
	    class(Graphene),intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey
	    integer, dimension(2), intent(in) :: range_localdim
	    call allocate_GrapheneVars(this,range_localdim)
	    call this%material%initialize_material_indexes(rangex,rangey,range_localdim)
	    call initialize_GrapheneVars(this)
        end subroutine initialize_Graphene

	subroutine allocate_GrapheneVars(this,range_localdim)
	    class(Graphene) :: this
	    integer,dimension(2),intent(in) :: range_localdim
	    integer :: alloc_stat,dimx,dimy
	    dimx=range_localdim(1)
	    dimy=range_localdim(2)
	    allocate(this%gax(dimx,dimy),stat=alloc_stat)
	    allocate(this%ha(dimx,dimy),stat=alloc_stat)
	    allocate(this%hb(dimx,dimy),stat=alloc_stat)
	    allocate(this%gb1(dimx,dimy),stat=alloc_stat)
	    allocate(this%gb2(dimx,dimy),stat=alloc_stat)
	    allocate(this%gb3(dimx,dimy),stat=alloc_stat)
	    allocate(this%gb4(dimx,dimy),stat=alloc_stat)
	    allocate(this%gc1(dimx,dimy),stat=alloc_stat)
	    allocate(this%gc2(dimx,dimy),stat=alloc_stat)
	    allocate(this%gc3(dimx,dimy),stat=alloc_stat)
	    allocate(this%gc4(dimx,dimy),stat=alloc_stat)
	    allocate(this%Snx(dimx,dimy),stat=alloc_stat)
	    allocate(this%Inx(dimx,dimy),stat=alloc_stat)
	    allocate(this%g1x(dimx,dimy),stat=alloc_stat)
	    allocate(this%g2x(dimx,dimy),stat=alloc_stat)
	    allocate(this%g3x(dimx,dimy),stat=alloc_stat)
	    allocate(this%g4x(dimx,dimy),stat=alloc_stat)
!***********************************************************
!		Nonlinear variables for one Raman peak
!***********************************************************
	    allocate(this%SLG_raman0(dimx,dimy),stat=alloc_stat)
	    allocate(this%SLG_raman1(dimx,dimy),stat=alloc_stat)
	    allocate(this%SLG_raman2(dimx,dimy),stat=alloc_stat)
	    allocate(this%SLG_Snx(dimx,dimy),stat=alloc_stat)
	    allocate(this%SLG_Snxm1(dimx,dimy),stat=alloc_stat)
	    allocate(this%SLG_Snxm2(dimx,dimy),stat=alloc_stat)
!***********************************************************
	    this%SLG_raman0=0.d0
	    this%SLG_raman1=0.d0
	    this%SLG_raman2=0.d0
	    this%SLG_Snx=0.d0
	    this%SLG_Snxm1=0.d0
	    this%SLG_Snxm2=0.d0
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

	    this%gax=1.d0
	end subroutine allocate_GrapheneVars

	subroutine initialize_GrapheneVars(this)
	class(Graphene) :: this
		integer :: x,y
	type(Irange) :: posx,posy
	real*8 :: xi3,wNL,dnuNL,dt,graphene_dnur,graphene_epinf
	complex*16 :: exp1,exp2,exp3,exp4,expi1,expi2,expi3,expi4

	call setexp(this)
	call scaleParams(this)
	posx=this%xrange
	posy=this%yrange
	xi3=this%xi3
	wNL=this%wNL
	dnuNL=this%dnuNL
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

!	if(this%scaleFact.eq.1.d0)then
!	    this%scaled_ap=this%ap
!	    this%scaled_graphene_wD=this%graphene_wD
!	else
!	    
!	endif

	    do x=posx%init,posx%ends,1
		do y=posy%init,posy%ends,1
		this%ha(x,y)=dexp(-graphene_dnur*dt)
		this%hb(x,y)=(this%scaled_graphene_wD**2)*dt/graphene_dnur
		this%gb1(x,y)=exp1
		this%gb2(x,y)=exp2
		this%gb3(x,y)=exp3
		this%gb4(x,y)=exp4
		this%gc1(x,y)=dt*2.d0*unoi*this%scaled_ap(1)*this%omp(1)*expi1
		this%gc2(x,y)=dt*2.d0*unoi*this%scaled_ap(2)*this%omp(2)*expi2
		this%gc3(x,y)=dt*2.d0*unoi*this%scaled_ap(3)*this%omp(3)*expi3
		this%gc4(x,y)=dt*2.d0*unoi*this%scaled_ap(4)*this%omp(4)*expi4
		this%gax(x,y)=graphene_epinf-this%gc1(x,y)-this%gc2(x,y)-this%gc3(x,y)-this%gc4(x,y)

		this%SLG_raman0(x,y)=xi3*((wNL*dt)**2)/(1.d0+dnuNL*dt)
		this%SLG_raman1(x,y)=(2.d0-(wNL*dt)**2)/(1.d0+dnuNL*dt)
		this%SLG_raman2(x,y)=-(1.d0-dnuNL*dt)/(1.d0+dnuNL*dt)
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

	    this%SLG_raman0=0.d0
	    this%SLG_raman1=0.d0
	    this%SLG_raman2=0.d0
	    this%SLG_Snx=0.d0
	    this%SLG_Snxm1=0.d0
	    this%SLG_Snxm2=0.d0
	end subroutine clearGrapheneVars


    function getRgax(this,posx,posy) result(val)
        class(Graphene) :: this
        integer, intent(in):: posx,posy
        real*8 :: val
        val=real(this%gax(posx,posy))
    end function getRgax

    function getIgax(this,posx,posy) result(val)
        class(Graphene) :: this
        integer, intent(in):: posx,posy
        real*8 :: val
        val=imag(this%gax(posx,posy))
    end function getIgax
	end module
