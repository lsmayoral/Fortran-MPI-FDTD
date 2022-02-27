	module simpledielectric_mod
	use mydatatypes
	use material_mod
	implicit none
	private
	public :: ConstantDielectric !,setDielectric_vars
	type, extends(material):: ConstantDielectric
	    real*8 :: constant_ep
	    real*8, allocatable,dimension(:,:),private :: gax,gay
	    contains
	    procedure :: update_E => update_E_inConstantDielectric
	    procedure :: initialize_material => initialize_ConstantDielectric
	    procedure :: initialize_ConstantDielectricVars
	    procedure :: allocate_ConstantDielectricVars
	    procedure :: getgax
	    procedure :: getgay
	end type ConstantDielectric

	contains 

subroutine update_E_inConstantDielectric(this,Ex_pointer,Ey_pointer,Dx_pointer,Dy_pointer)
	    class(ConstantDielectric) :: this
        real*8,pointer, intent(inout) :: Ex_pointer(:,:),Ey_pointer(:,:)
        real*8,pointer, intent(inout) :: Dx_pointer(:,:),Dy_pointer(:,:)
	    integer :: i,j

	    if(this%xrange%init.gt.0.and.this%xrange%ends.gt.0.and.&
		this%yrange%init.gt.0.and.this%yrange%ends.gt.0)then

!		call this%material%material_show_ranges()
		do i=this%xrange%init,this%xrange%ends,1
		    do j=this%yrange%init,this%yrange%ends,1
			Ex_pointer(i,j)=Dx_pointer(i,j)*this%gax(i,j)
			Ey_pointer(i,j)=Dy_pointer(i,j)*this%gay(i,j)
		    enddo
		enddo
	    endif
	end subroutine update_E_inConstantDielectric

	subroutine initialize_ConstantDielectric(this,rangex,rangey,range_localdim)
	    class(ConstantDielectric), intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey
	    integer, dimension(2), intent(in) :: range_localdim

	    call allocate_ConstantdielectricVars(this,range_localdim)
	    call this%material%initialize_material_indexes(rangex,rangey,range_localdim)
	    call initialize_ConstantDielectricVars(this)
        end subroutine initialize_ConstantDielectric

	subroutine allocate_ConstantDielectricVars(this,range_localdim)
	    class(ConstantDielectric) :: this
	    integer,dimension(2),intent(in) :: range_localdim
	    integer :: alloc_stat
	    allocate(this%gax(range_localdim(1),range_localdim(2)),stat=alloc_stat)
	    allocate(this%gay(range_localdim(1),range_localdim(2)),stat=alloc_stat)
	    this%gax=1.d0
	    this%gay=1.d0
	end subroutine allocate_ConstantDielectricVars

	subroutine initialize_ConstantDielectricVars(this)
	class(ConstantDielectric) :: this
	type(Irange) :: posx,posy
	integer :: x,y
	posx=this%xrange
	posy=this%yrange
	    do x=posx%init,posx%ends,1
		do y=posy%init,posy%ends,1
		    if(this%constant_ep.eq.0.d0)then
			this%gax(x,y)=0.d0
			this%gay(x,y)=0.d0
		    else
			this%gax(x,y)=1.d0/this%constant_ep
			this%gay(x,y)=1.d0/this%constant_ep
		    endif
		enddo
	    enddo
	end subroutine initialize_ConstantDielectricVars

	function getgax(this,posx,posy) result(val)
	    class(ConstantDielectric) :: this
	    integer, intent(in):: posx,posy
	    real*8 :: val
	    val=this%gax(posx,posy)
	end function getgax

	function getgay(this,posx,posy) result(val)
	    class(ConstantDielectric) :: this
	    integer, intent(in):: posx,posy
	    real*8 :: val
	    val=this%gay(posx,posy)
	end function getgay

	end module
