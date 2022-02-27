	module filmgoldLorentzJp_mod
	use constants_mod
	use mydatatypes
	use material_mod
	implicit none
	private
	public :: filmgold_LJp

	type, extends(material):: filmgold_LJp
	real*8, public :: dt,ds
!**************************************************************
!	Variables to update Electric field
!**************************************************************
	real*8, allocatable,dimension(:,:) :: ga,JpxT,JpyT
	real*8, allocatable,dimension(:,:,:) :: Px0,Px,Jpx
	real*8, allocatable,dimension(:,:,:) :: Py0,Py,Jpy
!**************************************************************
!		Parameter of 3 Lorenztian model for Gold to fit from 500nm to 800nm
!**************************************************************
	integer :: npols=2
	real*8 :: epinf=3.497071d0
	real*8, dimension(2) :: xi1=(/2.301113d0,400d0/)
	real*8, dimension(2) :: w0l=(2.d0*pi*c0/micro)*(/-2.393546d0,0.3441309d0/)
	real*8, dimension(2) :: dnu0l=(2.d0*pi*c0/micro)*(/0.2460773d0,0.02533959d0/)

	    contains
	    procedure :: update_EJp => update_E_inGold
	    procedure :: initialize_material => initialize_Gold
	    procedure :: initialize_Gold
	    procedure :: allocate_GoldVars
	end type filmgold_LJp

	contains 

subroutine update_E_inGold(this,Ex_pointer,Ey_pointer,Hz_pointer)
	    class(filmgold_LJp) :: this
        real*8,pointer, intent(inout) :: Ex_pointer(:,:),Ey_pointer(:,:)
        real*8,pointer, intent(inout) :: Hz_pointer(:,:)
	    integer :: i,j,l
	    real*8 :: dt,Ku,curlhzy,curlhzx,z0
	    if(this%xrange%init.gt.0.and.this%xrange%ends.gt.0.and.&
	       this%yrange%init.gt.0.and.this%yrange%ends.gt.0)then
	       dt=this%dt
	       z0=dsqrt(mu0/ep0)
		Ku=this%dt*c0/this%ds
		    do j=this%yrange%init,this%yrange%ends,1
			do i=this%xrange%init,this%xrange%ends,1

	curlhzy=Hz_pointer(i,j)-Hz_pointer(i,j-1)
	curlhzx=-Hz_pointer(i,j)+Hz_pointer(i-1,j)

        Ex_pointer(i,j)=Ex_pointer(i,j)+&
	    Ku*this%ga(i,j)*(curlhzy-this%ds*this%JpxT(i,j))
        Ey_pointer(i,j)=Ey_pointer(i,j)+&
	    Ku*this%ga(i,j)*(curlhzx-this%ds*this%JpyT(i,j))

	this%JpxT(i,j)=0.d0;this%JpyT(i,j)=0.d0
	do l=1,this%npols,1
    this%Px0(i,j,l)=this%xi1(l)*Ex_pointer(i,j)/c0
    this%Py0(i,j,l)=this%xi1(l)*Ey_pointer(i,j)/c0

    this%Jpx(i,j,l)=((1.d0-this%dnu0l(l)*dt)/(1.d0+this%dnu0l(l)*dt))*this%Jpx(i,j,l)+&
	((this%w0l(l)**2)*dt/(1.d0+this%dnu0l(l)*dt))*(this%Px0(i,j,l)-this%Px(i,j,l))
    this%Px(i,j,l)=this%Px(i,j,l)+dt*this%Jpx(i,j,l)
    this%JpxT(i,j)=this%JpxT(i,j)+this%Jpx(i,j,l)

    this%Jpy(i,j,l)=((1.d0-this%dnu0l(l)*dt)/(1.d0+this%dnu0l(l)*dt))*this%Jpy(i,j,l)+&
	((this%w0l(l)**2)*dt/(1.d0+this%dnu0l(l)*dt))*(this%Py0(i,j,l)-this%Py(i,j,l))
    this%Py(i,j,l)=this%Py(i,j,l)+dt*this%Jpy(i,j,l)
    this%JpyT(i,j)=this%JpyT(i,j)+this%Jpy(i,j,l)
		    enddo
		enddo
	enddo
	endif
	end subroutine update_E_inGold

	subroutine initialize_Gold(this,rangex,rangey,range_localdim)
	    class(filmgold_LJp),intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey
	    integer, dimension(2), intent(in) :: range_localdim
	    call allocate_GoldVars(this,range_localdim)
	    call this%material%initialize_material_indexes(rangex,rangey,range_localdim)
	    call initialize_GoldVars(this)
        end subroutine initialize_Gold

	subroutine allocate_GoldVars(this,range_localdim)
	    class(filmgold_LJp) :: this
	    integer,dimension(2),intent(in) :: range_localdim
	    integer :: alloc_stat,dimx,dimy,npols
	    dimx=range_localdim(1)
	    dimy=range_localdim(2)
	    npols=this%npols

	    allocate(this%ga(dimx,dimy),stat=alloc_stat)

	    allocate(this%Px0(dimx,dimy,npols),stat=alloc_stat)
	    allocate(this%Px(dimx,dimy,npols),stat=alloc_stat)
	    allocate(this%Jpx(dimx,dimy,npols),stat=alloc_stat)
	    allocate(this%JpxT(dimx,dimy),stat=alloc_stat)

	    allocate(this%Py0(dimx,dimy,npols),stat=alloc_stat)
	    allocate(this%Py(dimx,dimy,npols),stat=alloc_stat)
	    allocate(this%Jpy(dimx,dimy,npols),stat=alloc_stat)
	    allocate(this%JpyT(dimx,dimy),stat=alloc_stat)

	    this%ga=1.d0
	    this%Px0=0.d0; this%Px=0.d0; this%Jpx=0.d0;this%JpxT=0.d0
	    this%Py0=0.d0; this%Py=0.d0; this%Jpy=0.d0;this%JpyT=0.d0
	end subroutine allocate_GoldVars

	subroutine initialize_GoldVars(this)
	class(filmgold_LJp) :: this
	integer :: i,j
	type(Irange) :: posx,posy
	posx=this%xrange
	posy=this%yrange
	    do j=posy%init,posy%ends,1
	        do i=posx%init,posx%ends,1
		    this%ga(i,j)=1.d0/this%epinf
		enddo
	    enddo
	end subroutine initialize_GoldVars

	end module filmgoldLorentzJp_mod
