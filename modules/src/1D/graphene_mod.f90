	module graphene_mod
	use constants_mod
	use mydatatypes
	use material_mod
	implicit none
	private
	public :: Graphene,clearGrapheneVars

	type, extends(material):: Graphene
	real*8,public :: dt
	real*8,public :: ds
	real*8 :: scaleFact=1.d0
!**************************************************************
!	Variables to update Electric field
!**************************************************************
	complex*16, allocatable,dimension(:) :: gax
	real*8, allocatable,dimension(:) ::ha,hb,Snx,Inx,kerr0
	complex*16,allocatable,dimension(:) ::gb1,gb2,gb3,gb4
	complex*16,allocatable,dimension(:) ::gc1,gc2,gc3,gc4
	complex*16,allocatable,dimension(:) :: g1x,g2x,g3x,g4x
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

subroutine update_Exz_inGraphene(this,E_pointer,D_pointer)
	    class(Graphene) :: this
        real*8,pointer, intent(inout) :: E_pointer(:)
        real*8,pointer, intent(inout) :: D_pointer(:)
	    integer :: pos
	    if(this%posrange%init.gt.0.and.this%posrange%ends.gt.0)then
		do pos=this%posrange%init,this%posrange%ends,1
E_pointer(pos)=(D_pointer(pos)-this%Inx(pos)+this%ha(pos)*this%Snx(pos)+this%gb1(pos)*this%g1x(pos)+&
	this%gb2(pos)*this%g2x(pos)+this%gb3(pos)*this%g3x(pos)+this%gb4(pos)*this%g4x(pos))/&
	this%gax(pos)

	this%Inx(pos)=this%Inx(pos)+ this%hb(pos)*E_pointer(pos)
	this%Snx(pos)= this%ha(pos)*this%Snx(pos)+ this%hb(pos)*E_pointer(pos)
	this%g1x(pos)=this%gb1(pos)*this%g1x(pos)+this%gc1(pos)*E_pointer(pos)
	this%g2x(pos)=this%gb2(pos)*this%g2x(pos)+this%gc2(pos)*E_pointer(pos)
	this%g3x(pos)=this%gb3(pos)*this%g3x(pos)+this%gc3(pos)*E_pointer(pos)
	this%g4x(pos)=this%gb4(pos)*this%g4x(pos)+this%gc4(pos)*E_pointer(pos)
	    enddo
	endif
	end subroutine update_Exz_inGraphene

	subroutine initialize_Graphene(this,absrange,localdim)
	    class(Graphene),intent(inout) :: this
	    type(Irange), intent(in) :: absrange
	    integer, intent(in) :: localdim
	    call allocate_GrapheneVars(this,localdim)
	    call this%material%initialize_material_indexes(absrange,localdim)
	    call initialize_GrapheneVars(this)
        end subroutine initialize_Graphene

	subroutine allocate_GrapheneVars(this,localdim)
	    class(Graphene) :: this
	    integer, intent(in) :: localdim
	    integer :: alloc_stat

	    allocate(this%gax(localdim),stat=alloc_stat)
	    allocate(this%ha(localdim),stat=alloc_stat)
	    allocate(this%hb(localdim),stat=alloc_stat)
	    allocate(this%gb1(localdim),stat=alloc_stat)
	    allocate(this%gb2(localdim),stat=alloc_stat)
	    allocate(this%gb3(localdim),stat=alloc_stat)
	    allocate(this%gb4(localdim),stat=alloc_stat)
	    allocate(this%gc1(localdim),stat=alloc_stat)
	    allocate(this%gc2(localdim),stat=alloc_stat)
	    allocate(this%gc3(localdim),stat=alloc_stat)
	    allocate(this%gc4(localdim),stat=alloc_stat)
	    allocate(this%Snx(localdim),stat=alloc_stat)
	    allocate(this%Inx(localdim),stat=alloc_stat)
	    allocate(this%g1x(localdim),stat=alloc_stat)
	    allocate(this%g2x(localdim),stat=alloc_stat)
	    allocate(this%g3x(localdim),stat=alloc_stat)
	    allocate(this%g4x(localdim),stat=alloc_stat)

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
	integer :: y
	type(Irange) :: pos
	real*8 :: dt,graphene_dnur,graphene_epinf
	complex*16 :: exp1,exp2,exp3,exp4,expi1,expi2,expi3,expi4
	call setexp(this)
	call this%scaleParams()
	pos=this%posrange

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

	do y=pos%init,pos%ends,1
	    this%ha(y)=dexp(-graphene_dnur*dt)
	    this%hb(y)=(this%scaled_graphene_wD**2)*dt/graphene_dnur
	    this%gb1(y)=exp1
	    this%gb2(y)=exp2
	    this%gb3(y)=exp3
	    this%gb4(y)=exp4
	    this%gc1(y)=dt*2.d0*unoi*this%scaled_ap(1)*this%omp(1)*expi1
	    this%gc2(y)=dt*2.d0*unoi*this%scaled_ap(2)*this%omp(2)*expi2
	    this%gc3(y)=dt*2.d0*unoi*this%scaled_ap(3)*this%omp(3)*expi3
	    this%gc4(y)=dt*2.d0*unoi*this%scaled_ap(4)*this%omp(4)*expi4
	    this%gax(y)=graphene_epinf-this%gc1(y)-this%gc2(y)-this%gc3(y)-this%gc4(y)
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
	end subroutine clearGrapheneVars


    function getRgax(this,pos) result(val)
        class(Graphene) :: this
        integer, intent(in):: pos
        real*8 :: val
        val=real(this%gax(pos))
    end function getRgax

    function getIgax(this,pos) result(val)
        class(Graphene) :: this
        integer, intent(in):: pos
        real*8 :: val
        val=imag(this%gax(pos))
    end function getIgax
	end module
