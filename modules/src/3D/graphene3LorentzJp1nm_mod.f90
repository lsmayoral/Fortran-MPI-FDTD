	module graphene3LorentzJp1nm_mod
	use constants_mod
	use mydatatypes
	use material_mod
	implicit none
	private
	public :: SLG_3LJp1nm

	type, extends(material):: SLG_3LJp1nm
!**************************************************************
!	Nonlinear variables
!**************************************************************
	real*8,public :: wNL
	real*8,public :: dnuNL
	real*8,public :: xi3
	real*8,public :: dt
	real*8,public :: ds
!**************************************************************
!	Variables to update Electric field
!**************************************************************
	real*8, allocatable,dimension(:,:,:) :: ga
	real*8, allocatable,dimension(:,:,:,:) :: Px0,Px,Jpx,Pz0,Pz,Jpz
	real*8, allocatable,dimension(:,:,:) :: JpxT,JpzT,Gx,Gz,Qx,Qz
!**************************************************************
!		Parameter of 4 Lorenztian model for Graphene
!**************************************************************
	integer :: npols=3
	real*8 :: epinf=1.023242d0
	real*8, dimension(3) :: xi1=(/4.957366d0,4.455092d0,1.905429d0/)
	real*8, dimension(3) :: w0l=(2.d0*pi*c0/micro)*(/0.573075d0,1.713164d0,4.598546d0/)
	real*8, dimension(3) :: dnu0l=(2.d0*pi*c0/micro)*(/0.368926d0,2.42099d0,1.35581d0/)

	    contains
	    procedure :: update_EJp => update_Exz_inGraphene
	    procedure :: initialize_material => initialize_Graphene
	    procedure :: initialize_GrapheneVars
	    procedure :: allocate_GrapheneVars
	end type SLG_3LJp1nm

	contains 

subroutine update_Exz_inGraphene(this,Ex_pointer,Ey_pointer,Ez_pointer,Hx_pointer,Hy_pointer,Hz_pointer)
	    class(SLG_3LJp1nm) :: this
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
        Ey_pointer(i,j,k)=Ey_pointer(i,j,k)+Ku*(curlhxz+curlhzx)
        Ez_pointer(i,j,k)=Ez_pointer(i,j,k)+&
	    Ku*this%ga(i,j,k)*(curlhyx+curlhxy-this%ds*this%JpzT(i,j,k))

	this%JpxT(i,j,k)=0.d0;this%JpzT(i,j,k)=0.d0
	do l=1,this%npols,1
	    this%Px0(i,j,k,l)=this%xi1(l)*Ex_pointer(i,j,k)/c0
	    this%Pz0(i,j,k,l)=this%xi1(l)*Ez_pointer(i,j,k)/c0
	    if(l.eq.1)then
		this%Px0(i,j,k,l)=this%Px0(i,j,k,l)+this%xi3*this%Qx(i,j,k)*Ex_pointer(i,j,k)/c0
		this%Pz0(i,j,k,l)=this%Pz0(i,j,k,l)+this%xi3*this%Qz(i,j,k)*Ez_pointer(i,j,k)/c0
	    endif
    this%Jpx(i,j,k,l)=((1.d0-this%dnu0l(l)*dt)/(1.d0+this%dnu0l(l)*dt))*this%Jpx(i,j,k,l)+&
	((this%w0l(l)**2)*dt/(1.d0+this%dnu0l(l)*dt))*(this%Px0(i,j,k,l)-this%Px(i,j,k,l))
    this%Px(i,j,k,l)=this%Px(i,j,k,l)+dt*this%Jpx(i,j,k,l)
    this%JpxT(i,j,k)=this%JpxT(i,j,k)+this%Jpx(i,j,k,l)

    this%Jpz(i,j,k,l)=((1.d0-this%dnu0l(l)*dt)/(1.d0+this%dnu0l(l)*dt))*this%Jpz(i,j,k,l)+&
	((this%w0l(l)**2)*dt/(1.d0+this%dnu0l(l)*dt))*(this%Pz0(i,j,k,l)-this%Pz(i,j,k,l))
    this%Pz(i,j,k,l)=this%Pz(i,j,k,l)+dt*this%Jpz(i,j,k,l)
    this%JpzT(i,j,k)=this%JpzT(i,j,k)+this%Jpz(i,j,k,l)
	enddo

    this%Gx(i,j,k)=((1.d0-this%dnuNL*dt)/(1.d0+this%dnuNL*dt))*this%Gx(i,j,k)+&
	((this%wNL**2)*dt/(1.d0+this%dnuNL*dt))*(((z0*Ex_pointer(i,j,k))**2)-this%Qx(i,j,k))
    this%Qx(i,j,k)=this%Qx(i,j,k)+dt*this%Gx(i,j,k)

    this%Gz(i,j,k)=((1.d0-this%dnuNL*dt)/(1.d0+this%dnuNL*dt))*this%Gz(i,j,k)+&
        ((this%wNL**2)*dt/(1.d0+this%dnuNL*dt))*(((z0*Ez_pointer(i,j,k))**2)-this%Qz(i,j,k))
    this%Qz(i,j,k)=this%Qz(i,j,k)+dt*this%Gz(i,j,k)
		    enddo
		enddo
	    enddo
	endif
	end subroutine update_Exz_inGraphene

	subroutine initialize_Graphene(this,rangex,rangey,rangez,range_localdim)
	    class(SLG_3LJp1nm),intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey,rangez
	    integer, dimension(3), intent(in) :: range_localdim
	    call allocate_GrapheneVars(this,range_localdim)
	    call this%material%initialize_material_indexes(rangex,rangey,rangez,range_localdim)

	    call initialize_GrapheneVars(this)
        end subroutine initialize_Graphene

	subroutine allocate_GrapheneVars(this,range_localdim)
	    class(SLG_3LJp1nm) :: this
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
	    allocate(this%JPxT(dimx,dimy,dimz),stat=alloc_stat)

	    allocate(this%Pz0(dimx,dimy,dimz,npols),stat=alloc_stat)
	    allocate(this%Pz(dimx,dimy,dimz,npols),stat=alloc_stat)
	    allocate(this%Jpz(dimx,dimy,dimz,npols),stat=alloc_stat)
	    allocate(this%JPzT(dimx,dimy,dimz),stat=alloc_stat)

	    allocate(this%Gx(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Gz(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Qx(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Qz(dimx,dimy,dimz),stat=alloc_stat)
	    this%Gx=0.d0;this%Gz=0.d0;this%Qx=0.d0;this%Qz=0.d0
	    this%ga=1.d0
	    this%Px0=0.d0; this%Px=0.d0; this%Jpx=0.d0; this%JPxT=0.d0
	    this%Pz0=0.d0; this%Pz=0.d0; this%Jpz=0.d0; this%JPzT=0.d0
	end subroutine allocate_GrapheneVars

	subroutine initialize_GrapheneVars(this)
	class(SLG_3LJp1nm) :: this
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
	end subroutine initialize_GrapheneVars

	end module graphene3LorentzJp1nm_mod
