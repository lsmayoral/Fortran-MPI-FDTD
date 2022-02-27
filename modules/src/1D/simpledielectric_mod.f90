	module simpledielectric_mod
	use mydatatypes
	use material_mod
	implicit none
	private
	public :: ConstantDielectric
	type, extends(material):: ConstantDielectric
	    real*8 :: constant_ep
	    real*8, allocatable,dimension(:), public :: ga
	    contains
	    procedure :: update_E => update_E_inConstantDielectric
	    procedure :: initialize_material => initialize_ConstantDielectric
	    procedure :: initialize_ConstantDielectricVars
	    procedure :: allocate_ConstantDielectricVars
	    procedure :: getga
	end type ConstantDielectric

	contains 

subroutine update_E_inConstantDielectric(this,E_pointer,D_pointer)
	    class(ConstantDielectric) :: this
        real*8,pointer, intent(inout) :: E_pointer(:)
        real*8,pointer, intent(inout) :: D_pointer(:)
	integer :: pos

	    if(this%posrange%init.gt.0.and.this%posrange%ends.gt.0)then
		do pos=this%posrange%init,this%posrange%ends,1
		    E_pointer(pos)=D_pointer(pos)*this%ga(pos)
		enddo
	    endif
	end subroutine update_E_inConstantDielectric

	subroutine initialize_ConstantDielectric(this,absrange,localdim)
	    class(ConstantDielectric), intent(inout) :: this
	    type(Irange), intent(in) :: absrange
	    integer, intent(in) :: localdim

	    call allocate_ConstantdielectricVars(this,localdim)
	    call this%material%initialize_material_indexes(absrange,localdim)
	    call initialize_ConstantDielectricVars(this)
        end subroutine initialize_ConstantDielectric

	subroutine allocate_ConstantDielectricVars(this,localdim)
	    class(ConstantDielectric) :: this
	    integer,intent(in) :: localdim
	    integer :: alloc_stat
	    allocate(this%ga(localdim),stat=alloc_stat)
	    this%ga=1.d0
	end subroutine allocate_ConstantDielectricVars

	subroutine initialize_ConstantDielectricVars(this)
	class(ConstantDielectric) :: this
	integer :: pos

	    do pos=this%posrange%init,this%posrange%ends,1
		if(this%constant_ep.eq.0.d0)then
		    this%ga(pos)=0.d0
		else
		    this%ga(pos)=1.d0/this%constant_ep
		endif
	    enddo
	end subroutine initialize_ConstantDielectricVars

	function getga(this,pos) result(val)
	    class(ConstantDielectric) :: this
	    integer, intent(in):: pos
	    real*8 :: val
	    val=this%ga(pos)
	end function getga

	end module
