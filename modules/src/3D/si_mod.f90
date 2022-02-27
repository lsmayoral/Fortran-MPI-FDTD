	module si_mod
	use constants_mod
	use mydatatypes
	use material_mod
	implicit none
	private
	public :: Si

	type, extends(material):: Si
	real*8,public :: dt
!**************************************************************
!	Variables to update Electric field
!**************************************************************
	real*8, allocatable, dimension(:,:,:) :: ga1,ga2,ga3
	real*8, allocatable, dimension(:,:,:) :: gb1,gb2,gb3
	real*8, allocatable, dimension(:,:,:) :: gc1,gc2,gc3
	real*8, allocatable, dimension(:,:,:) :: Sn1x,Sn1y,Sn1z
	real*8, allocatable, dimension(:,:,:) :: Sn2x,Sn2y,Sn2z
	real*8, allocatable, dimension(:,:,:) :: Sn3x,Sn3y,Sn3z

	real*8, allocatable, dimension(:,:,:) :: Sn1m1x,Sn1m1y,Sn1m1z
	real*8, allocatable, dimension(:,:,:) :: Sn2m1x,Sn2m1y,Sn2m1z
	real*8, allocatable, dimension(:,:,:) :: Sn3m1x,Sn3m1y,Sn3m1z

	real*8, allocatable, dimension(:,:,:) :: Sn1m2x,Sn1m2y,Sn1m2z
	real*8, allocatable, dimension(:,:,:) :: Sn2m2x,Sn2m2y,Sn2m2z
	real*8, allocatable, dimension(:,:,:) :: Sn3m2x,Sn3m2y,Sn3m2z
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!
!c**************************************************************
!	real*8 :: fact=(2.d0*pi*c0/micro)
	real*8, dimension(3) :: dep=(/ 8.d0,2.85d0,-0.107d0 /)
	real*8, dimension(3) :: wp=(2.d0*pi*c0/micro)*(/ 3.64d0,2.76d0,1.73d0 /)
	real*8, dimension(3) :: dnup=(2.d0*pi*c0/micro)*(/ 0d0,0.063d0,2.5d0/)
	    contains
	    procedure :: update_E => update_E_inSi
	    procedure :: initialize_material => initialize_Si
	    procedure :: initialize_SiVars
	    procedure :: allocate_SiVars
	end type Si

	contains 

subroutine update_E_inSi(this,Ex_pointer,Ey_pointer,Ez_pointer,Dx_pointer,Dy_pointer,Dz_pointer)
	class(Si) :: this
	real*8, pointer, intent(inout):: Dx_pointer(:,:,:),Dy_pointer(:,:,:),Dz_pointer(:,:,:)
        real*8,pointer, intent(inout) :: Ex_pointer(:,:,:),Ey_pointer(:,:,:),Ez_pointer(:,:,:)
	integer :: i,j,k
	    if(this%xrange%init.gt.0.and.this%xrange%ends.gt.0.and.&
		this%yrange%init.gt.0.and.this%yrange%ends.gt.0.and.&
		this%zrange%init.gt.0.and.this%zrange%ends.gt.0)then
		do i=this%xrange%init,this%xrange%ends,1
		    do j=this%yrange%init,this%yrange%ends,1
			do k=this%zrange%init,this%zrange%ends,1
	Ex_pointer(i,j,k)=Dx_pointer(i,j,k)-this%Sn1x(i,j,k)-this%sn2x(i,j,k)-this%sn3x(i,j,k)
	this%Sn1x(i,j,k)=this%ga1(i,j,k)*this%Sn1m1x(i,j,k)+this%gb1(i,j,k)*this%Sn1m2x(i,j,k)+this%gc1(i,j,k)*Ex_pointer(i,j,k)
	this%Sn2x(i,j,k)=this%ga2(i,j,k)*this%Sn2m1x(i,j,k)+this%gb2(i,j,k)*this%Sn2m2x(i,j,k)+this%gc2(i,j,k)*Ex_pointer(i,j,k)
	this%Sn3x(i,j,k)=this%ga3(i,j,k)*this%Sn3m1x(i,j,k)+this%gb3(i,j,k)*this%Sn3m2x(i,j,k)+this%gc3(i,j,k)*Ex_pointer(i,j,k)
	    this%Sn1m2x(i,j,k)=this%Sn1m1x(i,j,k); this%Sn2m2x(i,j,k)=this%Sn2m1x(i,j,k); this%Sn3m2x(i,j,k)=this%Sn3m1x(i,j,k)
	    this%Sn1m1x(i,j,k)=this%Sn1x(i,j,k); this%Sn2m1x(i,j,k)=this%Sn2x(i,j,k); this%Sn3m1x(i,j,k)=this%Sn3x(i,j,k)

	Ey_pointer(i,j,k)=Dy_pointer(i,j,k)-this%Sn1y(i,j,k)-this%sn2y(i,j,k)-this%sn3y(i,j,k)
	this%Sn1y(i,j,k)=this%ga1(i,j,k)*this%Sn1m1y(i,j,k)+this%gb1(i,j,k)*this%Sn1m2y(i,j,k)+this%gc1(i,j,k)*Ey_pointer(i,j,k)
	this%Sn2y(i,j,k)=this%ga2(i,j,k)*this%Sn2m1y(i,j,k)+this%gb2(i,j,k)*this%Sn2m2y(i,j,k)+this%gc2(i,j,k)*Ey_pointer(i,j,k)
	this%Sn3y(i,j,k)=this%ga3(i,j,k)*this%Sn3m1y(i,j,k)+this%gb3(i,j,k)*this%Sn3m2y(i,j,k)+this%gc3(i,j,k)*Ey_pointer(i,j,k)
	    this%Sn1m2y(i,j,k)=this%Sn1m1y(i,j,k); this%Sn2m2y(i,j,k)=this%Sn2m1y(i,j,k); this%Sn3m2y(i,j,k)=this%Sn3m1y(i,j,k)
	    this%Sn1m1y(i,j,k)=this%Sn1y(i,j,k); this%Sn2m1y(i,j,k)=this%Sn2y(i,j,k); this%Sn3m1y(i,j,k)=this%Sn3y(i,j,k)

	Ez_pointer(i,j,k)=Dz_pointer(i,j,k)-this%Sn1z(i,j,k)-this%sn2z(i,j,k)-this%sn3z(i,j,k)
	this%Sn1z(i,j,k)=this%ga1(i,j,k)*this%Sn1m1z(i,j,k)+this%gb1(i,j,k)*this%Sn1m2z(i,j,k)+this%gc1(i,j,k)*Ez_pointer(i,j,k)
	this%Sn2z(i,j,k)=this%ga2(i,j,k)*this%Sn2m1z(i,j,k)+this%gb2(i,j,k)*this%Sn2m2z(i,j,k)+this%gc2(i,j,k)*Ez_pointer(i,j,k)
	this%Sn3z(i,j,k)=this%ga3(i,j,k)*this%Sn3m1z(i,j,k)+this%gb3(i,j,k)*this%Sn3m2z(i,j,k)+this%gc3(i,j,k)*Ez_pointer(i,j,k)
	    this%Sn1m2z(i,j,k)=this%Sn1m1z(i,j,k); this%Sn2m2z(i,j,k)=this%Sn2m1z(i,j,k); this%Sn3m2z(i,j,k)=this%Sn3m1z(i,j,k)
	    this%Sn1m1z(i,j,k)=this%Sn1z(i,j,k); this%Sn2m1z(i,j,k)=this%Sn2z(i,j,k); this%Sn3m1z(i,j,k)=this%Sn3z(i,j,k)

			enddo
		    enddo
		enddo
	    endif
	end subroutine update_E_inSi

	subroutine initialize_Si(this,rangex,rangey,rangez,range_localdim)
	    class(Si),intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey,rangez
	    integer, dimension(3), intent(in) :: range_localdim
	    call allocate_SiVars(this,range_localdim)
	    call this%material%initialize_material_indexes(rangex,rangey,rangez,range_localdim)
	    call initialize_SiVars(this)
        end subroutine initialize_Si

	subroutine allocate_SiVars(this,range_localdim)
	    class(Si) :: this
	    integer,dimension(3),intent(in) :: range_localdim
	    integer :: alloc_stat,dimx,dimy,dimz
	    dimx=range_localdim(1)
	    dimy=range_localdim(2)
	    dimz=range_localdim(3)

	    allocate(this%ga1(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%ga2(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%ga3(dimx,dimy,dimz),stat=alloc_stat)

	    allocate(this%gb1(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gb2(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gb3(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gc1(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gc2(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%gc3(dimx,dimy,dimz),stat=alloc_stat)

	    allocate(this%Sn1x(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn1y(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn1z(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn2x(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn2y(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn2z(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn3x(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn3y(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn3z(dimx,dimy,dimz),stat=alloc_stat)

	    allocate(this%Sn1m1x(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn1m1y(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn1m1z(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn2m1x(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn2m1y(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn2m1z(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn3m1x(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn3m1y(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn3m1z(dimx,dimy,dimz),stat=alloc_stat)

	    allocate(this%Sn1m2x(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn1m2y(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn1m2z(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn2m2x(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn2m2y(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn2m2z(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn3m2x(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn3m2y(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Sn3m2z(dimx,dimy,dimz),stat=alloc_stat)

	    this%ga1=0.d0
	    this%ga2=0.d0
	    this%ga3=0.d0
	    this%gb1=0.d0
	    this%gb2=0.d0
	    this%gb3=0.d0
	    this%gc1=0.d0
	    this%gc2=0.d0
	    this%gc3=0.d0
	    this%Sn1x=0.d0
	    this%Sn1y=0.d0
	    this%Sn1z=0.d0
	    this%Sn2x=0.d0
	    this%Sn2y=0.d0
	    this%Sn2z=0.d0
	    this%Sn3x=0.d0
	    this%Sn3y=0.d0
	    this%Sn3z=0.d0
	    this%Sn1m1x=0.d0
	    this%Sn1m1y=0.d0
	    this%Sn1m1z=0.d0
	    this%Sn2m1x=0.d0
	    this%Sn2m1y=0.d0
	    this%Sn2m1z=0.d0
	    this%Sn3m1x=0.d0
	    this%Sn3m1y=0.d0
	    this%Sn3m1z=0.d0
	    this%Sn1m2x=0.d0
	    this%Sn1m2y=0.d0
	    this%Sn1m2z=0.d0
	    this%Sn2m2x=0.d0
	    this%Sn2m2y=0.d0
	    this%Sn2m2z=0.d0
	    this%Sn3m2x=0.d0
	    this%Sn3m2y=0.d0
	    this%Sn3m2z=0.d0
	end subroutine allocate_SiVars

	subroutine initialize_SiVars(this)
	class(Si) :: this
	integer :: x,y,z
	type(Irange) :: posx,posy,posz
	real*8 :: dt

	posx=this%xrange
	posy=this%yrange
	posz=this%zrange
	dt=this%dt

	    do x=posx%init,posx%ends,1
		do y=posy%init,posy%ends,1
		    do z=posz%init,posz%ends,1
			this%ga1(x,y,z)=(2.d0-(dt*this%wp(1))**2)/(1.d0+dt*this%dnup(1))
			this%ga2(x,y,z)=(2.d0-(dt*this%wp(2))**2)/(1.d0+dt*this%dnup(2))
			this%ga3(x,y,z)=(2.d0-(dt*this%wp(3))**2)/(1.d0+dt*this%dnup(3))
			this%gb1(x,y,z)=(dt*this%dnup(1)-1.d0)/(1.d0+dt*this%dnup(1))
			this%gb2(x,y,z)=(dt*this%dnup(2)-1.d0)/(1.d0+dt*this%dnup(2))
			this%gb3(x,y,z)=(dt*this%dnup(3)-1.d0)/(1.d0+dt*this%dnup(3))
			this%gc1(x,y,z)=(this%dep(1)*(this%wp(1)*dt)**2)/(1.d0+dt*this%dnup(1))
			this%gc2(x,y,z)=(this%dep(2)*(this%wp(2)*dt)**2)/(1.d0+dt*this%dnup(2))
			this%gc3(x,y,z)=(this%dep(3)*(this%wp(3)*dt)**2)/(1.d0+dt*this%dnup(3))
		    enddo
		enddo
	    enddo
	end subroutine initialize_SiVars


	end module


