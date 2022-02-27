	module simpledielectricLorentzJp_mod
	use mydatatypes
	use constants_mod
	use material_mod
	implicit none
	private
	public :: ConstantDielectricLJp !,setDielectric_vars
	type, extends(material):: ConstantDielectricLJp
	    real*8,public :: constant_ep,dt,ds
	    real*8, allocatable,dimension(:,:), public :: ga
	    contains
	    procedure :: update_EJp => update_E_inConstantDielectric
	    procedure :: initialize_material => initialize_ConstantDielectric
	    procedure :: initialize_ConstantDielectricVars
	    procedure :: allocate_ConstantDielectricVars
	end type ConstantDielectricLJp

	contains 

subroutine update_E_inConstantDielectric(this,Ex_pointer,Ey_pointer,Hz_pointer)
	    class(ConstantDielectricLJp) :: this
        real*8,pointer, intent(inout) :: Ex_pointer(:,:),Ey_pointer(:,:)
        real*8,pointer, intent(inout) :: Hz_pointer(:,:)
	real*8 :: Ku,curlhzy,curlhyz,curlhxz,curlhzx,curlhyx,curlhxy
	integer :: i,j,k

	    if(this%xrange%init.gt.0.and.this%xrange%ends.gt.0.and.&
		this%yrange%init.gt.0.and.this%yrange%ends.gt.0)then
		Ku=this%dt*c0/this%ds
		do i=this%xrange%init,this%xrange%ends,1
		    do j=this%yrange%init,this%yrange%ends,1
		curlhzy=Hz_pointer(i,j)-Hz_pointer(i,j-1)
		curlhzx=-Hz_pointer(i,j)+Hz_pointer(i-1,j)

		Ex_pointer(i,j)=Ex_pointer(i,j)+Ku*this%ga(i,j)*curlhzy
		Ey_pointer(i,j)=Ey_pointer(i,j)+Ku*this%ga(i,j)*curlhzx
			enddo
		    enddo
	    endif
	end subroutine update_E_inConstantDielectric

	subroutine initialize_ConstantDielectric(this,rangex,rangey,range_localdim)
	    class(ConstantDielectricLJp), intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey
	    integer, dimension(2), intent(in) :: range_localdim

	    call allocate_ConstantdielectricVars(this,range_localdim)
	    call this%material%initialize_material_indexes(rangex,rangey,range_localdim)
	    call initialize_ConstantDielectricVars(this)
        end subroutine initialize_ConstantDielectric

	subroutine allocate_ConstantDielectricVars(this,range_localdim)
	    class(ConstantDielectricLJp) :: this
	    integer,dimension(2),intent(in) :: range_localdim
	    integer :: alloc_stat,nx,ny
	    nx=range_localdim(1)
	    ny=range_localdim(2)

	    allocate(this%ga(nx,ny),stat=alloc_stat)
	    this%ga=1.d0
	end subroutine allocate_ConstantDielectricVars

	subroutine initialize_ConstantDielectricVars(this)
	class(ConstantDielectricLJp) :: this
	type(Irange) :: posx,posy
	integer :: x,y
	posx=this%xrange
	posy=this%yrange

	    do x=posx%init,posx%ends,1
		do y=posy%init,posy%ends,1
			if(this%constant_ep.eq.0.d0)then
			    this%ga(x,y)=0.d0
			else
			    this%ga(x,y)=1.d0/this%constant_ep
			endif
		    enddo
		enddo
	end subroutine initialize_ConstantDielectricVars

	end module simpledielectricLorentzJp_mod
