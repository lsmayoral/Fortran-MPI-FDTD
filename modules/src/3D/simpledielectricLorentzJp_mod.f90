	module simpledielectricLorentzJp_mod
	use mydatatypes
	use constants_mod
	use material_mod
	implicit none
	private
	public :: ConstantDielectricLJp !,setDielectric_vars
	type, extends(material):: ConstantDielectricLJp
	    real*8,public :: constant_ep,dt,ds
	    real*8, allocatable,dimension(:,:,:), public :: ga
	    contains
	    procedure :: update_EJp => update_E_inConstantDielectric
	    procedure :: initialize_material => initialize_ConstantDielectric
	    procedure :: initialize_ConstantDielectricVars
	    procedure :: allocate_ConstantDielectricVars
	end type ConstantDielectricLJp

	contains 

subroutine update_E_inConstantDielectric(this,Ex_pointer,Ey_pointer,Ez_pointer,Hx_pointer,Hy_pointer,Hz_pointer)
	    class(ConstantDielectricLJp) :: this
        real*8,pointer, intent(inout) :: Ex_pointer(:,:,:),Ey_pointer(:,:,:),Ez_pointer(:,:,:)
        real*8,pointer, intent(inout) :: Hx_pointer(:,:,:),Hy_pointer(:,:,:),Hz_pointer(:,:,:)
	real*8 :: Ku,curlhzy,curlhyz,curlhxz,curlhzx,curlhyx,curlhxy
	integer :: i,j,k

	    if(this%xrange%init.gt.0.and.this%xrange%ends.gt.0.and.&
		this%yrange%init.gt.0.and.this%yrange%ends.gt.0.and.&
		this%zrange%init.gt.0.and.this%zrange%ends.gt.0)then
		Ku=this%dt*c0/this%ds
		do i=this%xrange%init,this%xrange%ends,1
		    do j=this%yrange%init,this%yrange%ends,1
			do k=this%zrange%init,this%zrange%ends,1
		curlhzy=Hz_pointer(i,j,k)-Hz_pointer(i,j-1,k)
		curlhyz=-Hy_pointer(i,j,k)+Hy_pointer(i,j,k-1)
		curlhxz=Hx_pointer(i,j,k)-Hx_pointer(i,j,k-1)
		curlhzx=-Hz_pointer(i,j,k)+Hz_pointer(i-1,j,k)
		curlhyx=Hy_pointer(i,j,k)-Hy_pointer(i-1,j,k)
		curlhxy=-Hx_pointer(i,j,k)+Hx_pointer(i,j-1,k)

		Ex_pointer(i,j,k)=Ex_pointer(i,j,k)+Ku*this%ga(i,j,k)*(curlhzy+curlhyz)
		Ey_pointer(i,j,k)=Ey_pointer(i,j,k)+Ku*this%ga(i,j,k)*(curlhxz+curlhzx)
		Ez_pointer(i,j,k)=Ez_pointer(i,j,k)+Ku*this%ga(i,j,k)*(curlhyx+curlhxy)
			enddo
		    enddo
		enddo
	    endif
	end subroutine update_E_inConstantDielectric

	subroutine initialize_ConstantDielectric(this,rangex,rangey,rangez,range_localdim)
	    class(ConstantDielectricLJp), intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey,rangez
	    integer, dimension(3), intent(in) :: range_localdim

	    call allocate_ConstantdielectricVars(this,range_localdim)
	    call this%material%initialize_material_indexes(rangex,rangey,rangez,range_localdim)
	    call initialize_ConstantDielectricVars(this)
        end subroutine initialize_ConstantDielectric

	subroutine allocate_ConstantDielectricVars(this,range_localdim)
	    class(ConstantDielectricLJp) :: this
	    integer,dimension(3),intent(in) :: range_localdim
	    integer :: alloc_stat,nx,ny,nz
	    nx=range_localdim(1)
	    ny=range_localdim(2)
	    nz=range_localdim(3)

	    allocate(this%ga(nx,ny,nz),stat=alloc_stat)
	    this%ga=1.d0
	end subroutine allocate_ConstantDielectricVars

	subroutine initialize_ConstantDielectricVars(this)
	class(ConstantDielectricLJp) :: this
	type(Irange) :: posx,posy,posz
	integer :: x,y,z
	posx=this%xrange
	posy=this%yrange
	posz=this%zrange

	    do x=posx%init,posx%ends,1
		do y=posy%init,posy%ends,1
		    do z=posz%init,posz%ends,1
			if(this%constant_ep.eq.0.d0)then
			    this%ga(x,y,z)=0.d0
			else
			    this%ga(x,y,z)=1.d0/this%constant_ep
			endif
		    enddo
		enddo
	    enddo
	end subroutine initialize_ConstantDielectricVars

	end module simpledielectricLorentzJp_mod
