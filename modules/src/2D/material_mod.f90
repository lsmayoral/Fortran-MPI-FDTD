	module material_mod
	use mydatatypes
	implicit none
	private
!	public::material,initialize_material,update_E
	public material,material_show_ranges
	type material
	    integer :: xinit,yinit,xend,yend
	    type(Irange) :: xrange,yrange
	    type(Irange) :: rankxrange,rankyrange
	    contains 
		procedure :: initialize_material
		procedure :: initialize_material_indexes
		procedure :: material_show_ranges
		procedure :: update_E
	end type material

	contains 
	subroutine material_show_ranges(this)
	    class(material),intent(in) :: this
	    write(*,*)"*************************"
	    write(*,*) this%xrange
	    write(*,*) this%yrange
	    write(*,*)"*************************"
	end subroutine material_show_ranges


	subroutine update_E(this,Ex_pointer,Ey_pointer,Dx_pointer,Dy_pointer)
	    class(material) :: this
	    real*8,pointer, intent(inout) :: Ex_pointer(:,:),Ey_pointer(:,:)
	    real*8,pointer, intent(inout) :: Dx_pointer(:,:),Dy_pointer(:,:)
	    print *,"this never happend"
!	    integer :: i,j
!	    if(this%xrange%init.gt.0.and.this%xrange%ends.gt.0.and.&
!		this%yrange%init.gt.0.and.this%yrange%ends.gt.0)then
!		do i=this%xrange%init,this%xrange%ends,1
!		    do j=this%yrange%init,this%yrange%ends,1
!		    Ex_pointer(i,j)=Dx_pointer(i,j)
!		    Ey_pointer(i,j)=Dy_pointer(i,j)
!		    enddo
!		enddo
!	    endif
	end subroutine update_E

	subroutine initialize_material(this,rangex,rangey,range_localdim)
	    class(material), intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey
	    integer, dimension(2), intent(in) :: range_localdim
	
	    call initialize_material_indexes(this,rangex,rangey,range_localdim)
	end subroutine initialize_material

	subroutine initialize_material_indexes(this,rangex,rangey,range_localdim)
	    class(material), intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey
	    integer, dimension(2), intent(in) :: range_localdim
	    integer :: mat_xinit,mat_xend,mat_yinit,mat_yend
	    integer :: posx,posy

		this%xrange%init=0
		this%xrange%ends=-1
		this%yrange%init=0
		this%yrange%ends=-1

		this%rankxrange=rangex
		this%rankyrange=rangey

	    if(rangex%init.le.this%xend.and.rangex%ends.ge.this%xinit.and. &
		rangey%init.le.this%yend.and.rangey%ends.ge.this%yinit)then
	        mat_xinit=1
	        mat_xend=range_localdim(1)
	        mat_yinit=1
	        mat_yend=range_localdim(2)

		if(rangex%init.le.this%xinit)then
		    mat_xinit=this%xinit-rangex%init+1
		endif
		if(rangex%ends.ge.this%xend)then
		    mat_xend=this%xend-rangex%init+1
		endif
		if(rangey%init.le.this%yinit)then
		    mat_yinit=this%yinit-rangey%init+1
		endif
		if(rangey%ends.ge.this%yend)then
		    mat_yend=this%yend-rangey%init+1
		endif
		this%xrange%init=mat_xinit
		this%xrange%ends=mat_xend
		this%yrange%init=mat_yinit
		this%yrange%ends=mat_yend
	    endif
!	    call material_show_ranges(this)
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
