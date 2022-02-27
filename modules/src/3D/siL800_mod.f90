	module siL800_mod
	use constants_mod
	use mydatatypes
	use material_mod
	implicit none
	private
	public :: SiL800

	type, extends(material):: SiL800
	real*8,public :: dt
!**************************************************************
!	Variables to update Electric field
!**************************************************************
	real*8, allocatable, dimension(:,:,:) :: ga,gb
	real*8, allocatable, dimension(:,:,:) :: Inx,Iny,Inz
	real*8 :: sigma=112.8134634d0
	real*8 :: ep_inf=(3.6750d0)**2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	    contains
	    procedure :: update_E => update_E_inSi
	    procedure :: initialize_material => initialize_Si
	    procedure :: initialize_SiVars
	    procedure :: allocate_SiVars
	end type SiL800

	contains 

subroutine update_E_inSi(this,Ex_pointer,Ey_pointer,Ez_pointer,Dx_pointer,Dy_pointer,Dz_pointer)
	class(SiL800) :: this
	real*8, pointer, intent(inout):: Dx_pointer(:,:,:),Dy_pointer(:,:,:),Dz_pointer(:,:,:)
        real*8,pointer, intent(inout) :: Ex_pointer(:,:,:),Ey_pointer(:,:,:),Ez_pointer(:,:,:)
	integer :: i,j,k
	    if(this%xrange%init.gt.0.and.this%xrange%ends.gt.0.and.&
		this%yrange%init.gt.0.and.this%yrange%ends.gt.0.and.&
		this%zrange%init.gt.0.and.this%zrange%ends.gt.0)then
		do k=this%zrange%init,this%zrange%ends,1
		    do j=this%yrange%init,this%yrange%ends,1
			do i=this%xrange%init,this%xrange%ends,1
	Ex_pointer(i,j,k)=(Dx_pointer(i,j,k)-this%Inx(i,j,k))/this%ga(i,j,k)
	this%Inx(i,j,k)=this%Inx(i,j,k)+this%gb(i,j,k)*Ex_pointer(i,j,k)

	Ey_pointer(i,j,k)=(Dy_pointer(i,j,k)-this%Iny(i,j,k))/this%ga(i,j,k)
	this%Iny(i,j,k)=this%Iny(i,j,k)+this%gb(i,j,k)*Ey_pointer(i,j,k)

	Ey_pointer(i,j,k)=(Dy_pointer(i,j,k)-this%Iny(i,j,k))/this%ga(i,j,k)
	this%Iny(i,j,k)=this%Iny(i,j,k)+this%gb(i,j,k)*Ey_pointer(i,j,k)
			enddo
		    enddo
		enddo
	    endif
	end subroutine update_E_inSi

	subroutine initialize_Si(this,rangex,rangey,rangez,range_localdim)
	    class(SiL800),intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey,rangez
	    integer, dimension(3), intent(in) :: range_localdim
	    call allocate_SiVars(this,range_localdim)
	    call this%material%initialize_material_indexes(rangex,rangey,rangez,range_localdim)
	    call initialize_SiVars(this)
        end subroutine initialize_Si

	subroutine allocate_SiVars(this,range_localdim)
	    class(SiL800) :: this
	    integer,dimension(3),intent(in) :: range_localdim
	    integer :: alloc_stat,dimx,dimy,dimz
	    dimx=range_localdim(1)
	    dimy=range_localdim(2)
	    dimz=range_localdim(3)

	    allocate(this%ga(dimx,dimy,dimz),stat=alloc_stat)

	    allocate(this%gb(dimx,dimy,dimz),stat=alloc_stat)

	    allocate(this%Inx(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Iny(dimx,dimy,dimz),stat=alloc_stat)
	    allocate(this%Inz(dimx,dimy,dimz),stat=alloc_stat)

	    this%ga=1.d0
	    this%gb=0.d0
	    this%Inx=0.d0
	    this%Iny=0.d0
	    this%Inz=0.d0
	end subroutine allocate_SiVars

	subroutine initialize_SiVars(this)
	class(SiL800) :: this
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
			this%ga(x,y,z)=this%ep_inf+this%sigma*dt/ep0
			this%gb(x,y,z)=this%sigma*dt/ep0
		    enddo
		enddo
	    enddo
	end subroutine initialize_SiVars

	end module siL800_mod


