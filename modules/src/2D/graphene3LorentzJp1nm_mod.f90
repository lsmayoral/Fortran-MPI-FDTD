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
	real*8, allocatable,dimension(:,:) :: ga
	real*8, allocatable,dimension(:,:,:) :: Px0,Px,Jpx
	real*8, allocatable,dimension(:,:) :: JpxT,Gx,Qx
!**************************************************************
!		Parameter of 4 Lorenztian model for Graphene
!**************************************************************
	integer :: npols=3
	real*8 :: epinf=1.023242d0
	real*8, dimension(3) :: xi1=(/4.957366d0,4.455092d0,1.905429d0/)
	real*8, dimension(3) :: w0l=(2.d0*pi*c0/micro)*(/0.573075d0,1.713164d0,4.598546d0/)
	real*8, dimension(3) :: dnu0l=(2.d0*pi*c0/micro)*(/0.368926d0,2.42099d0,1.35581d0/)

	    contains
	    procedure :: update_EJp => update_Ex_inGraphene
	    procedure :: initialize_material => initialize_Graphene
	    procedure :: initialize_GrapheneVars
	    procedure :: allocate_GrapheneVars
	end type SLG_3LJp1nm

	contains 

subroutine update_Ex_inGraphene(this,Ex_pointer,Ey_pointer,Hz_pointer)
	    class(SLG_3LJp1nm) :: this
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
        Ey_pointer(i,j)=Ey_pointer(i,j)+Ku*curlhzx

	this%JpxT(i,j)=0.d0
	do l=1,this%npols,1
	    this%Px0(i,j,l)=this%xi1(l)*Ex_pointer(i,j)/c0
	    if(l.eq.1)then
		this%Px0(i,j,l)=this%Px0(i,j,l)+this%xi3*this%Qx(i,j)*Ex_pointer(i,j)/c0
	    endif
    this%Jpx(i,j,l)=((1.d0-this%dnu0l(l)*dt)/(1.d0+this%dnu0l(l)*dt))*this%Jpx(i,j,l)+&
	((this%w0l(l)**2)*dt/(1.d0+this%dnu0l(l)*dt))*(this%Px0(i,j,l)-this%Px(i,j,l))
    this%Px(i,j,l)=this%Px(i,j,l)+dt*this%Jpx(i,j,l)
    this%JpxT(i,j)=this%JpxT(i,j)+this%Jpx(i,j,l)
	enddo

    this%Gx(i,j)=((1.d0-this%dnuNL*dt)/(1.d0+this%dnuNL*dt))*this%Gx(i,j)+&
	((this%wNL**2)*dt/(1.d0+this%dnuNL*dt))*(((z0*Ex_pointer(i,j))**2)-this%Qx(i,j))
    this%Qx(i,j)=this%Qx(i,j)+dt*this%Gx(i,j)
		    enddo
		enddo
	endif
	end subroutine update_Ex_inGraphene

	subroutine initialize_Graphene(this,rangex,rangey,range_localdim)
	    class(SLG_3LJp1nm),intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey
	    integer, dimension(2), intent(in) :: range_localdim
	    call allocate_GrapheneVars(this,range_localdim)
	    call this%material%initialize_material_indexes(rangex,rangey,range_localdim)

	    call initialize_GrapheneVars(this)
        end subroutine initialize_Graphene

	subroutine allocate_GrapheneVars(this,range_localdim)
	    class(SLG_3LJp1nm) :: this
	    integer,dimension(2),intent(in) :: range_localdim
	    integer :: alloc_stat,dimx,dimy,npols
	    dimx=range_localdim(1)
	    dimy=range_localdim(2)
	    npols=this%npols

	    allocate(this%ga(dimx,dimy),stat=alloc_stat)

	    allocate(this%Px0(dimx,dimy,npols),stat=alloc_stat)
	    allocate(this%Px(dimx,dimy,npols),stat=alloc_stat)
	    allocate(this%Jpx(dimx,dimy,npols),stat=alloc_stat)
	    allocate(this%JPxT(dimx,dimy),stat=alloc_stat)

	    allocate(this%Gx(dimx,dimy),stat=alloc_stat)
	    allocate(this%Qx(dimx,dimy),stat=alloc_stat)
	    this%Gx=0.d0;this%Qx=0.d0
	    this%ga=1.d0
	    this%Px0=0.d0; this%Px=0.d0; this%Jpx=0.d0; this%JPxT=0.d0
	end subroutine allocate_GrapheneVars

	subroutine initialize_GrapheneVars(this)
	class(SLG_3LJp1nm) :: this
	integer :: i,j
	type(Irange) :: posx,posy
	posx=this%xrange
	posy=this%yrange
		do j=posy%init,posy%ends,1
		    do i=posx%init,posx%ends,1
		    this%ga(i,j)=1.d0/this%epinf
		enddo
	    enddo
	end subroutine initialize_GrapheneVars

	end module graphene3LorentzJp1nm_mod
