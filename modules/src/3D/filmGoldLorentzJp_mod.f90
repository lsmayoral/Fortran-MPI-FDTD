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
	real*8, allocatable,dimension(:,:,:) :: ga,JpxT,JpyT,JpzT
	real*8, allocatable,dimension(:,:,:,:) :: Px0,Px,Jpx
	real*8, allocatable,dimension(:,:,:,:) :: Py0,Py,Jpy
	real*8, allocatable,dimension(:,:,:,:) :: Pz0,Pz,Jpz
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

subroutine update_E_inGold(this,Ex_pointer,Ey_pointer,Ez_pointer,Hx_pointer,Hy_pointer,Hz_pointer)
	    class(filmgold_LJp) :: this
        real*8,pointer, intent(inout) :: Ex_pointer(:,:,:),Ey_pointer(:,:,:),Ez_pointer(:,:,:)
        real*8,pointer, intent(inout) :: Hx_pointer(:,:,:),Hy_pointer(:,:,:),Hz_pointer(:,:,:)
	    integer :: i,j,k,l
	    real*8 :: dt,Ku,curlhzy,curlhyz,curlhxz,curlhzx,curlhyx,curlhxy,z0
	    if(this%xrange%init.gt.0.and.this%xrange%ends.gt.0.and.&
	       this%yrange%init.gt.0.and.this%yrange%ends.gt.0.and.&
	       this%zrange%init.gt.0.and.this%zrange%ends.gt.0)then
	       dt=this%dt
	       z0=dsqrt(mu0/ep0)
		Ku=this%dt*c0/this%ds
		do k=this%zrange%init,this%zrange%ends,1
		    do j=this%yrange%init,this%yrange%ends,1
			do i=this%xrange%init,this%xrange%ends,1

	curlhzy=Hz_pointer(i,j,k)-Hz_pointer(i,j-1,k)
	curlhyz=-Hy_pointer(i,j,k)+Hy_pointer(i,j,k-1)
	curlhxz=Hx_pointer(i,j,k)-Hx_pointer(i,j,k-1)
	curlhzx=-Hz_pointer(i,j,k)+Hz_pointer(i-1,j,k)
	curlhyx=Hy_pointer(i,j,k)-Hy_pointer(i-1,j,k)
        curlhxy=-Hx_pointer(i,j,k)+Hx_pointer(i,j-1,k)

        Ex_pointer(i,j,k)=Ex_pointer(i,j,k)+&
	    Ku*this%ga(i,j,k)*(curlhzy+curlhyz-this%ds*this%JpxT(i,j,k))
        Ey_pointer(i,j,k)=Ey_pointer(i,j,k)+&
	    Ku*this%ga(i,j,k)*(curlhxz+curlhzx-this%ds*this%JpyT(i,j,k))
        Ez_pointer(i,j,k)=Ez_pointer(i,j,k)+&
	    Ku*this%ga(i,j,k)*(curlhyx+curlhxy-this%ds*this%JpzT(i,j,k))

	this%JpxT(i,j,k)=0.d0;this%JpyT(i,j,k)=0.d0;this%JpzT(i,j,k)=0.d0
	do l=1,this%npols,1
    this%Px0(i,j,k,l)=this%xi1(l)*Ex_pointer(i,j,k)/c0
    this%Py0(i,j,k,l)=this%xi1(l)*Ey_pointer(i,j,k)/c0
    this%Pz0(i,j,k,l)=this%xi1(l)*Ez_pointer(i,j,k)/c0

    this%Jpx(i,j,k,l)=((1.d0-this%dnu0l(l)*dt)/(1.d0+this%dnu0l(l)*dt))*this%Jpx(i,j,k,l)+&
	((this%w0l(l)**2)*dt/(1.d0+this%dnu0l(l)*dt))*(this%Px0(i,j,k,l)-this%Px(i,j,k,l))
    this%Px(i,j,k,l)=this%Px(i,j,k,l)+dt*this%Jpx(i,j,k,l)
    this%JpxT(i,j,k)=this%JpxT(i,j,k)+this%Jpx(i,j,k,l)

    this%Jpy(i,j,k,l)=((1.d0-this%dnu0l(l)*dt)/(1.d0+this%dnu0l(l)*dt))*this%Jpy(i,j,k,l)+&
	((this%w0l(l)**2)*dt/(1.d0+this%dnu0l(l)*dt))*(this%Py0(i,j,k,l)-this%Py(i,j,k,l))
    this%Py(i,j,k,l)=this%Py(i,j,k,l)+dt*this%Jpy(i,j,k,l)
    this%JpyT(i,j,k)=this%JpyT(i,j,k)+this%Jpy(i,j,k,l)

    this%Jpz(i,j,k,l)=((1.d0-this%dnu0l(l)*dt)/(1.d0+this%dnu0l(l)*dt))*this%Jpz(i,j,k,l)+&
	((this%w0l(l)**2)*dt/(1.d0+this%dnu0l(l)*dt))*(this%Pz0(i,j,k,l)-this%Pz(i,j,k,l))
    this%Pz(i,j,k,l)=this%Pz(i,j,k,l)+dt*this%Jpz(i,j,k,l)
    this%JpzT(i,j,k)=this%JpzT(i,j,k)+this%Jpz(i,j,k,l)
	enddo
		    enddo
		enddo
	    enddo
	endif
	end subroutine update_E_inGold

	subroutine initialize_Gold(this,rangex,rangey,rangez,range_localdim)
	    class(filmgold_LJp),intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey,rangez
	    integer, dimension(3), intent(in) :: range_localdim
	    call allocate_GoldVars(this,range_localdim)
	    call this%material%initialize_material_indexes(rangex,rangey,rangez,range_localdim)
	    call initialize_GoldVars(this)
        end subroutine initialize_Gold

	subroutine allocate_GoldVars(this,range_localdim)
	    class(filmgold_LJp) :: this
	    integer,dimension(3),intent(in) :: range_localdim
	    integer :: alloc_stat,dimx,dimy,dimz,npols
	    dimx=range_localdim(1)
	    dimy=range_localdim(2)
	    dimz=range_localdim(3)
	    npols=this%npols

	    allocate(this%ga(dimx,dimy,dimz),stat=alloc_stat)

	    allocate(this%Px0(dimx,dimy,dimz,npols),stat=alloc_stat)
	    allocate(this%Px(dimx,dimy,dimz,npols),stat=alloc_stat)
	    allocate(this%Jpx(dimx,dimy,dimz,npols),stat=alloc_stat)
	    allocate(this%JpxT(dimx,dimy,dimz),stat=alloc_stat)

	    allocate(this%Py0(dimx,dimy,dimz,npols),stat=alloc_stat)
	    allocate(this%Py(dimx,dimy,dimz,npols),stat=alloc_stat)
	    allocate(this%Jpy(dimx,dimy,dimz,npols),stat=alloc_stat)
	    allocate(this%JpyT(dimx,dimy,dimz),stat=alloc_stat)

	    allocate(this%Pz0(dimx,dimy,dimz,npols),stat=alloc_stat)
	    allocate(this%Pz(dimx,dimy,dimz,npols),stat=alloc_stat)
	    allocate(this%Jpz(dimx,dimy,dimz,npols),stat=alloc_stat)
	    allocate(this%JpzT(dimx,dimy,dimz),stat=alloc_stat)

	    this%ga=1.d0
	    this%Px0=0.d0; this%Px=0.d0; this%Jpx=0.d0;this%JpxT=0.d0
	    this%Py0=0.d0; this%Py=0.d0; this%Jpy=0.d0;this%JpyT=0.d0
	    this%Pz0=0.d0; this%Pz=0.d0; this%Jpz=0.d0;this%JpzT=0.d0
	end subroutine allocate_GoldVars

	subroutine initialize_GoldVars(this)
	class(filmgold_LJp) :: this
	integer :: i,j,k
	type(Irange) :: posx,posy,posz
	posx=this%xrange
	posy=this%yrange
	posz=this%zrange
	    do k=posz%init,posz%ends,1
		do j=posy%init,posy%ends,1
		    do i=posx%init,posx%ends,1
		    this%ga(i,j,k)=1.d0/this%epinf
		enddo
	    enddo
	enddo
	end subroutine initialize_GoldVars

	end module filmgoldLorentzJp_mod
