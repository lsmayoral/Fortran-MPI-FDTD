	module material_mod
	use mydatatypes
	implicit none
	private
	public material,material_show_ranges,update_E_inVacuum
	type material
	    integer :: posinit,posend
	    type(Irange),public :: posrange
	    integer :: npos
	    contains 
		procedure :: initialize_material
		procedure :: initialize_material_indexes
		procedure :: material_show_ranges
		procedure :: update_E => update_E_inVacuum
	end type material

	contains 
	subroutine material_show_ranges(this,msg)
	    class(material),intent(in) :: this
	    character(*), intent(in) :: msg
	    write(*,*)"*************************"
	    write(*,*) msg
	    write(*,*) this%posrange
	    write(*,*)"*************************"
	end subroutine material_show_ranges


	subroutine update_E_inVacuum(this,E_pointer,D_pointer)
	    class(material) :: this
	    real*8,pointer, intent(inout) :: E_pointer(:)
	    real*8,pointer, intent(inout) :: D_pointer(:)
	    integer :: pos

	    if(this%posrange%init.gt.0.and.this%posrange%ends.gt.0)then
		do pos=this%posrange%init,this%posrange%ends,1
		    E_pointer(pos)=D_pointer(pos)
		enddo
	    endif
	end subroutine update_E_inVacuum

	subroutine initialize_material(this,absrange,localdim)
	    class(material), intent(inout) :: this
	    type(Irange), intent(in) :: absrange
	    integer, intent(in) :: localdim
	
	    call initialize_material_indexes(this,absrange,localdim)
	end subroutine initialize_material

	subroutine initialize_material_indexes(this,absrange,localdim)
	    class(material), intent(inout) :: this
	    type(Irange), intent(in) :: absrange
	    integer, intent(in) :: localdim
	    integer :: mat_init,mat_end
	    integer :: pos

		this%npos=localdim

		this%posrange%init=0
		this%posrange%ends=-1
	    if(absrange%init.le.this%posend.and.absrange%ends.ge.this%posinit)then
	        mat_init=1
	        mat_end=localdim

		if(absrange%init.le.this%posinit)then
		    mat_init=this%posinit-absrange%init+1
		endif
		if(absrange%ends.ge.this%posend)then
		    mat_end=this%posend-absrange%init+1
		endif

		this%posrange%init=mat_init
		this%posrange%ends=mat_end
	    endif
	    end subroutine initialize_material_indexes

!	    subroutine initialize_none_vars(this,posx,posy)
!	    subroutine initialize_vars(this,posx,posy)
!		class(material) :: this
!		type(Irange), intent(in) :: posx,posy
!	    end subroutine initialize_none_vars
!
!	    subroutine allocate_None(this,range_localdim)
!	    subroutine allocate_materialVars(this,range_localdim)
!		class(material) :: this
!		integer,dimension(2),intent(in) :: range_localdim
!		print *,"not allocating..."
!	    end subroutine allocate_None

	end module
