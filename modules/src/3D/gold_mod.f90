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
	real*8, allocatable, dimension(:,:,:) :: ha,hb
	real*8, allocatable, dimension(:,:,:) :: Snx,Sny,Snz
	real*8, allocatable, dimension(:,:,:) :: Iny,Inx,Inz
	complex*16,allocatable,dimension(:,:,:) ::ga
	complex*16,allocatable,dimension(:,:,:) ::gb1,gb2
	complex*16,allocatable,dimension(:,:,:) ::gc1,gc2
	complex*16,allocatable,dimension(:,:,:) ::g1x,g2x
	complex*16,allocatable,dimension(:,:,:) ::g1y,g2y
	complex*16,allocatable,dimension(:,:,:) ::g1z,g2z
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
	    procedure :: update_E => update_Ex_inGold
	    procedure :: initialize_material => initialize_Gold
	    procedure :: initialize_GoldVars
	    procedure :: allocate_GoldVars
	    procedure ::setexp
	    procedure ::clearGoldVars
	    procedure ::getRga
	    procedure ::getIga
	end type Gold

	contains 

subroutine update_Ex_inGold(this,Ex_pointer,Ey_pointer,Ez_pointer,Dx_pointer,Dy_pointer,Dz_pointer)
	class(Gold) :: this
	real*8, pointer, intent(inout):: Dx_pointer(:,:,:),Dy_pointer(:,:,:),Dz_pointer(:,:,:)
        real*8,pointer, intent(inout) :: Ex_pointer(:,:,:),Ey_pointer(:,:,:),Ez_pointer(:,:,:)
	integer :: i,j,k
	    if(this%xrange%init.gt.0.and.this%xrange%ends.gt.0.and.&
		this%yrange%init.gt.0.and.this%yrange%ends.gt.0.and.&
		this%zrange%init.gt.0.and.this%zrange%ends.gt.0)then
		do i=this%xrange%init,this%xrange%ends,1
		    do j=this%yrange%init,this%yrange%ends,1
			do k=this%zrange%init,this%zrange%ends,1
Ex_pointer(i,j,k)=(Dx_pointer(i,j,k)-this%Inx(i,j,k)+this%ha(i,j,k)*this%Snx(i,j,k)+this%gb1(i,j,k)*this%g1x(i,j,k)+&
	this%gb2(i,j,k)*this%g2x(i,j,k))/this%ga(i,j,k)

Ey_pointer(i,j,k)=(Dy_pointer(i,j,k)-this%Iny(i,j,k)+this%ha(i,j,k)*this%Sny(i,j,k)+this%gb1(i,j,k)*this%g1y(i,j,k)+&
	this%gb2(i,j,k)*this%g2y(i,j,k))/this%ga(i,j,k)

Ez_pointer(i,j,k)=(Dz_pointer(i,j,k)-this%Inz(i,j,k)+this%ha(i,j,k)*this%Snz(i,j,k)+this%gb1(i,j,k)*this%g1z(i,j,k)+&
	this%gb2(i,j,k)*this%g2z(i,j,k))/this%ga(i,j,k)

	this%Inx(i,j,k)=this%Inx(i,j,k)+this%hb(i,j,k)*Ex_pointer(i,j,k)
	this%Snx(i,j,k)=this%ha(i,j,k)*this%Snx(i,j,k)+this%hb(i,j,k)*Ex_pointer(i,j,k)
	this%g1x(i,j,k)=this%gb1(i,j,k)*this%g1x(i,j,k)+this%gc1(i,j,k)*Ex_pointer(i,j,k)
	this%g2x(i,j,k)=this%gb2(i,j,k)*this%g2x(i,j,k)+this%gc2(i,j,k)*Ex_pointer(i,j,k)

	this%Iny(i,j,k)=this%Iny(i,j,k)+this%hb(i,j,k)*Ey_pointer(i,j,k)
	this%Sny(i,j,k)=this%ha(i,j,k)*this%Sny(i,j,k)+this%hb(i,j,k)*Ey_pointer(i,j,k)
	this%g1y(i,j,k)=this%gb1(i,j,k)*this%g1y(i,j,k)+this%gc1(i,j,k)*Ey_pointer(i,j,k)
	this%g2y(i,j,k)=this%gb2(i,j,k)*this%g2y(i,j,k)+this%gc2(i,j,k)*Ey_pointer(i,j,k)

	this%Inz(i,j,k)=this%Inz(i,j,k)+this%hb(i,j,k)*Ez_pointer(i,j,k)
	this%Snz(i,j,k)=this%ha(i,j,k)*this%Snz(i,j,k)+this%hb(i,j,k)*Ez_pointer(i,j,k)
	this%g1z(i,j,k)=this%gb1(i,j,k)*this%g1z(i,j,k)+this%gc1(i,j,k)*Ez_pointer(i,j,k)
	this%g2z(i,j,k)=this%gb2(i,j,k)*this%g2z(i,j,k)+this%gc2(i,j,k)*Ez_pointer(i,j,k)
			enddo
		    enddo
		enddo
	    endif
	end subroutine update_Ex_inGold

	subroutine initialize_Gold(this,rangex,rangey,rangez,range_localdim)
	    class(Gold),intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey,rangez
	    integer, dimension(3), intent(in) :: range_localdim
	    call allocate_GoldVars(this,range_localdim)
	    call this%material%initialize_material_indexes(rangex,rangey,rangez,range_localdim)
	    call initialize_GoldVars(this)
        end subroutine initialize_Gold

	subroutine allocate_GoldVars(this,range_localdim)
	    class(Gold) :: this
	    integer,dimension(3),intent(in) :: range_localdim
	    integer :: alloc_stat,dimx,dimy,dimz
	    dimx=range_localdim(1)
	    dimy=range_localdim(2)
	    dimz=range_localdim(3)
	    allocate(this%ga(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%ha(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%hb(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gb1(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gb2(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gc1(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gc2(dimx,dimy,dimz),stat=alloc_stat)

	    allocate(this%Snx(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sny(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Snz(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Inx(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Iny(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Inz(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%g1x(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%g1y(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%g1z(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%g2x(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%g2y(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%g2z(dimx,dimy,dimz),stat=alloc_stat)
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

	    this%Snz=0.d0
	    this%Inz=0.d0
	    this%g1z=zero
	    this%g2z=zero

	    this%ga=1.d0
	end subroutine allocate_GoldVars

	subroutine initialize_GoldVars(this)
	class(Gold) :: this
		integer :: x,y,z
	type(Irange) :: posx,posy,posz
	real*8 :: dt,dnur,epinf
	complex*16 :: exp1,exp2,expi1,expi2

	call setexp(this)
	posx=this%xrange
	posy=this%yrange
	posz=this%zrange
	dt=this%dt
	dnur=this%dnur
	exp1=this%exp1
	exp2=this%exp2
	expi1=this%expi1
	expi2=this%expi2
	epinf=this%epinf

	    do x=posx%init,posx%ends,1
		do y=posy%init,posy%ends,1
		    do z=posz%init,posz%ends,1
		this%ha(x,y,z)=dexp(-dnur*dt)
		this%hb(x,y,z)=(this%wD**2)*dt/dnur
		this%gb1(x,y,z)=exp1
		this%gb2(x,y,z)=exp2
		this%gc1(x,y,z)=dt*2.d0*unoi*this%ap(1)*this%omp(1)*expi1
		this%gc2(x,y,z)=dt*2.d0*unoi*this%ap(2)*this%omp(2)*expi2
		this%ga(x,y,z)=epinf-this%gc1(x,y,z)-this%gc2(x,y,z)
		    enddo
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


    function getRga(this,posx,posy,posz) result(val)
        class(Gold) :: this
        integer, intent(in):: posx,posy,posz
        real*8 :: val
        val=real(this%ga(posx,posy,posz))
    end function getRga

    function getIga(this,posx,posy,posz) result(val)
        class(Gold) :: this
        integer, intent(in):: posx,posy,posz
        real*8 :: val
        val=imag(this%ga(posx,posy,posz))
    end function getIga
	end module


