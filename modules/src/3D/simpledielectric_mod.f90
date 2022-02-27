	module simpledielectric_mod
	use mydatatypes
	use material_mod
	implicit none
	private
	public :: ConstantDielectric !,setDielectric_vars
	type, extends(material):: ConstantDielectric
	    real*8 :: constant_ep
	    real*8, allocatable,dimension(:,:,:), public :: gax,gay,gaz
	    contains
	    procedure :: update_E => update_E_inConstantDielectric
	    procedure :: initialize_material => initialize_ConstantDielectric
	    procedure :: initialize_ConstantDielectricVars
	    procedure :: allocate_ConstantDielectricVars
	    procedure :: getgax
	    procedure :: getgay
	    procedure :: getgaz
	end type ConstantDielectric

	contains 

subroutine update_E_inConstantDielectric(this,Ex_pointer,Ey_pointer,Ez_pointer,Dx_pointer,Dy_pointer,Dz_pointer)
	    class(ConstantDielectric) :: this
        real*8,pointer, intent(inout) :: Ex_pointer(:,:,:),Ey_pointer(:,:,:),Ez_pointer(:,:,:)
        real*8,pointer, intent(inout) :: Dx_pointer(:,:,:),Dy_pointer(:,:,:),Dz_pointer(:,:,:)
	integer :: i,j,k

	    if(this%xrange%init.gt.0.and.this%xrange%ends.gt.0.and.&
		this%yrange%init.gt.0.and.this%yrange%ends.gt.0.and.&
		this%zrange%init.gt.0.and.this%zrange%ends.gt.0)then
!		write(*,*)"**********************************"
!		write(*,*)this%myrank,"Exassociated=",associated(Ex_pointer)
!		write(*,*)this%myrank,"Eyassociated=",associated(Ey_pointer)
!		write(*,*)this%myrank,"Ezassociated=",associated(Ez_pointer)
!		write(*,*)this%myrank,"Dxassociated=",associated(Dx_pointer)
!		write(*,*)this%myrank,"Dyassociated=",associated(Dy_pointer)
!		write(*,*)this%myrank,"Dzassociated=",associated(Dz_pointer)
!		write(*,*)this%myrank,"xrange:", this%xrange
!		write(*,*)this%myrank,"yrange:", this%yrange
!		write(*,*)this%myrank,"zrange:",this%zrange
!		write(*,*)"**********************************"
		do i=this%xrange%init,this%xrange%ends,1
		    do j=this%yrange%init,this%yrange%ends,1
			do k=this%zrange%init,this%zrange%ends,1
			Ex_pointer(i,j,k)=Dx_pointer(i,j,k)*this%gax(i,j,k)
			Ey_pointer(i,j,k)=Dy_pointer(i,j,k)*this%gay(i,j,k)
			Ez_pointer(i,j,k)=Dz_pointer(i,j,k)*this%gaz(i,j,k)
			enddo
		    enddo
		enddo
	    endif
	end subroutine update_E_inConstantDielectric

	subroutine initialize_ConstantDielectric(this,rangex,rangey,rangez,range_localdim)
	    class(ConstantDielectric), intent(inout) :: this
	    type(Irange), intent(in) :: rangex,rangey,rangez
	    integer, dimension(3), intent(in) :: range_localdim

	    call allocate_ConstantdielectricVars(this,range_localdim)
	    call this%material%initialize_material_indexes(rangex,rangey,rangez,range_localdim)
	    call initialize_ConstantDielectricVars(this)
        end subroutine initialize_ConstantDielectric

	subroutine allocate_ConstantDielectricVars(this,range_localdim)
	    class(ConstantDielectric) :: this
	    integer,dimension(3),intent(in) :: range_localdim
	    integer :: alloc_stat,nx,ny,nz
	    nx=range_localdim(1)
	    ny=range_localdim(2)
	    nz=range_localdim(3)

	    allocate(this%gax(nx,ny,nz),stat=alloc_stat)
	    allocate(this%gay(nx,ny,nz),stat=alloc_stat)
	    allocate(this%gaz(nx,ny,nz),stat=alloc_stat)
	    this%gax=1.d0
	    this%gay=1.d0
	    this%gaz=1.d0
	end subroutine allocate_ConstantDielectricVars

	subroutine initialize_ConstantDielectricVars(this)
	class(ConstantDielectric) :: this
	type(Irange) :: posx,posy,posz
	integer :: x,y,z
	posx=this%xrange
	posy=this%yrange
	posz=this%zrange

	    do x=posx%init,posx%ends,1
		do y=posy%init,posy%ends,1
		    do z=posz%init,posz%ends,1
			if(this%constant_ep.eq.0.d0)then
			    this%gax(x,y,z)=0.d0
			    this%gay(x,y,z)=0.d0
			    this%gaz(x,y,z)=0.d0
			else
			    this%gax(x,y,z)=1.d0/this%constant_ep
			    this%gay(x,y,z)=1.d0/this%constant_ep
			    this%gaz(x,y,z)=1.d0/this%constant_ep
			endif
		    enddo
		enddo
	    enddo
	end subroutine initialize_ConstantDielectricVars

	function getgax(this,posx,posy,posz) result(val)
	    class(ConstantDielectric) :: this
	    integer, intent(in):: posx,posy,posz
	    real*8 :: val
	    val=this%gax(posx,posy,posz)
	end function getgax

	function getgay(this,posx,posy,posz) result(val)
	    class(ConstantDielectric) :: this
	    integer, intent(in):: posx,posy,posz
	    real*8 :: val
	    val=this%gay(posx,posy,posz)
	end function getgay

	function getgaz(this,posx,posy,posz) result(val)
	    class(ConstantDielectric) :: this
	    integer, intent(in):: posx,posy,posz
	    real*8 :: val
	    val=this%gaz(posx,posy,posz)
	end function getgaz
	end module
