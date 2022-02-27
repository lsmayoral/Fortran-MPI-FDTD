	module gold_mod
	use constants_mod
	use mydatatypes
	use material_mod
	implicit none
	private
	public :: Gold,clearGoldVars,getRga,getIga

	type, extends(material):: Gold
	real*8,public :: dt
!**************************************************************
!	Variables to update Electric field
!**************************************************************
	real*8, allocatable, dimension(:,:) :: ha,hb,Snx,Inx,Sny,Iny
	complex*16,allocatable,dimension(:,:) ::ga
	complex*16,allocatable,dimension(:,:) ::gb1,gb2
	complex*16,allocatable,dimension(:,:) ::gc1,gc2
	complex*16,allocatable,dimension(:,:) ::g1x,g2x
	complex*16,allocatable,dimension(:,:) ::g1y,g2y
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c		Parameter for Gold was taken from
!c		Appl Phys A (2011) 103: 849-853
!c		DOI 10.1007/s00339-010-6224-9
!c**************************************************************
	real*8 :: epinf=1.1431d0
	real*8 :: dnur=1.0805d14
	real*8 :: wD=1.3202d16
	real*8, dimension(2) :: ap= (/ 0.26698d0,3.0834d0 /)
	real*8, dimension(2) :: phip=(/ -1.2371d0,-1.0968d0 /)
	real*8, dimension(2) :: omp=(/ 3.8711d15,4.1684d15 /)
	real*8, dimension(2) :: rp=(/ 4.4642d14,2.4555d15 /)
	complex*16 :: expi1,expi2,exp1,exp2

	    contains
	    procedure :: update_E => update_E_inGold
	    procedure :: initialize_material => initialize_Gold
	    procedure :: initialize_GoldVars
	    procedure :: allocate_GoldVars
	    procedure ::setexp
	    procedure ::clearGoldVars
	    procedure ::getRga
	    procedure ::getIga
	end type Gold

	contains 

subroutine update_E_inGold(this,Ex_pointer,Ey_pointer,Dx_pointer,Dy_pointer)
	    class(Gold) :: this
	real*8, pointer, intent(inout) ::Ey_pointer(:,:),Dy_pointer(:,:)
        real*8,pointer, intent(inout) :: Ex_pointer(:,:)
        real*8,pointer, intent(inout) :: Dx_pointer(:,:)
	    integer :: i,j
	    if(this%xrange%init.gt.0.and.this%xrange%ends.gt.0.and.&
		this%yrange%init.gt.0.and.this%yrange%ends.gt.0)then
		do i=this%xrange%init,this%xrange%ends,1
		    do j=this%yrange%init,this%yrange%ends,1
Ex_pointer(i,j)=(Dx_pointer(i,j)-this%Inx(i,j)+this%ha(i,j)*this%Snx(i,j)+this%gb1(i,j)*this%g1x(i,j)+&
	this%gb2(i,j)*this%g2x(i,j))/this%ga(i,j)

Ey_pointer(i,j)=(Dy_pointer(i,j)-this%Iny(i,j)+this%ha(i,j)*this%Sny(i,j)+this%gb1(i,j)*this%g1y(i,j)+&
	this%gb2(i,j)*this%g2y(i,j))/this%ga(i,j)

	this%Inx(i,j)=this%Inx(i,j)+this%hb(i,j)*Ex_pointer(i,j)
	this%Snx(i,j)=this%ha(i,j)*this%Snx(i,j)+this%hb(i,j)*Ex_pointer(i,j)
	this%g1x(i,j)=this%gb1(i,j)*this%g1x(i,j)+this%gc1(i,j)*Ex_pointer(i,j)
	this%g2x(i,j)=this%gb2(i,j)*this%g2x(i,j)+this%gc2(i,j)*Ex_pointer(i,j)

	this%Iny(i,j)=this%Iny(i,j)+this%hb(i,j)*Ey_pointer(i,j)
	this%Sny(i,j)=this%ha(i,j)*this%Sny(i,j)+this%hb(i,j)*Ey_pointer(i,j)
	this%g1y(i,j)=this%gb1(i,j)*this%g1y(i,j)+this%gc1(i,j)*Ey_pointer(i,j)
	this%g2y(i,j)=this%gb2(i,j)*this%g2y(i,j)+this%gc2(i,j)*Ey_pointer(i,j)
		    enddo
		enddo
	    endif
	end subroutine update_E_inGold

	subroutine initialize_Gold(this,rangex,rangey,range_localdim)
	    class(Gold),intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey
	    integer, dimension(2), intent(in) :: range_localdim
	    call allocate_GoldVars(this,range_localdim)
	    call this%material%initialize_material_indexes(rangex,rangey,range_localdim)
	    call initialize_GoldVars(this)
        end subroutine initialize_Gold

	subroutine allocate_GoldVars(this,range_localdim)
	    class(Gold) :: this
	    integer,dimension(2),intent(in) :: range_localdim
	    integer :: alloc_stat,dimx,dimy
	    dimx=range_localdim(1)
	    dimy=range_localdim(2)
	    allocate(this%ga(dimx,dimy),stat=alloc_stat)
	    allocate(this%ha(dimx,dimy),stat=alloc_stat)
	    allocate(this%hb(dimx,dimy),stat=alloc_stat)
	    allocate(this%gb1(dimx,dimy),stat=alloc_stat)
	    allocate(this%gb2(dimx,dimy),stat=alloc_stat)
	    allocate(this%gc1(dimx,dimy),stat=alloc_stat)
	    allocate(this%gc2(dimx,dimy),stat=alloc_stat)
	    allocate(this%Snx(dimx,dimy),stat=alloc_stat)
	    allocate(this%Inx(dimx,dimy),stat=alloc_stat)
	    allocate(this%g1x(dimx,dimy),stat=alloc_stat)
	    allocate(this%g2x(dimx,dimy),stat=alloc_stat)
	    allocate(this%Sny(dimx,dimy),stat=alloc_stat)
	    allocate(this%Iny(dimx,dimy),stat=alloc_stat)
	    allocate(this%g1y(dimx,dimy),stat=alloc_stat)
	    allocate(this%g2y(dimx,dimy),stat=alloc_stat)
!***********************************************************
	    this%ha=0.d0
	    this%hb=0.d0
	    this%gb1=zero
	    this%gb2=zero

	    this%gc1=zero
	    this%gc2=zero

	    this%Snx=0.d0
	    this%Inx=0.d0
	    this%g1x=zero
	    this%g2x=zero

	    this%Sny=0.d0
	    this%Iny=0.d0
	    this%g1y=zero
	    this%g2y=zero

	    this%ga=1.d0
	end subroutine allocate_GoldVars

	subroutine initialize_GoldVars(this)
	class(Gold) :: this
		integer :: x,y
	type(Irange) :: posx,posy
	real*8 :: dt,dnur,epinf
	complex*16 :: exp1,exp2,expi1,expi2

	call setexp(this)
	posx=this%xrange
	posy=this%yrange
	dt=this%dt
	dnur=this%dnur
	exp1=this%exp1
	exp2=this%exp2
	expi1=this%expi1
	expi2=this%expi2
	epinf=this%epinf

	    do x=posx%init,posx%ends,1
		do y=posy%init,posy%ends,1
		this%ha(x,y)=dexp(-dnur*dt)
		this%hb(x,y)=(this%wD**2)*dt/dnur
		this%gb1(x,y)=exp1
		this%gb2(x,y)=exp2
		this%gc1(x,y)=dt*2.d0*unoi*this%ap(1)*this%omp(1)*expi1
		this%gc2(x,y)=dt*2.d0*unoi*this%ap(2)*this%omp(2)*expi2
		this%ga(x,y)=epinf-this%gc1(x,y)-this%gc2(x,y)
		enddo
	    enddo
	end subroutine initialize_GoldVars

	subroutine setexp(this)
	    class(Gold) :: this

	    this%expi1=cdexp(-unoi*this%phip(1))
	    this%expi2=cdexp(-unoi*this%phip(2))

	    this%exp1=cdexp((-this%rp(1)+unoi*this%omp(1))*this%dt)
	    this%exp2=cdexp((-this%rp(2)+unoi*this%omp(2))*this%dt)
	end subroutine setexp

	subroutine clearGoldVars(this)
	    class(Gold) :: this
	    this%ha=0.d0
	    this%hb=0.d0
	    this%Snx=0.d0
	    this%Inx=0.d0
	    this%Sny=0.d0
	    this%Iny=0.d0
	    this%gb1=zero
	    this%gb2=zero
	    this%gc1=zero
	    this%gc2=zero
	    this%g1x=zero
	    this%g2x=zero
	end subroutine clearGoldVars


    function getRga(this,posx,posy) result(val)
        class(Gold) :: this
        integer, intent(in):: posx,posy
        real*8 :: val
        val=real(this%ga(posx,posy))
    end function getRga

    function getIga(this,posx,posy) result(val)
        class(Gold) :: this
        integer, intent(in):: posx,posy
        real*8 :: val
        val=imag(this%ga(posx,posy))
    end function getIga
	end module


