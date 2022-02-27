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
	real*8, allocatable, dimension(:) :: ga1,ga2,ga3
	real*8, allocatable, dimension(:) :: gb1,gb2,gb3
	real*8, allocatable, dimension(:) :: gc1,gc2,gc3
	real*8, allocatable, dimension(:) :: Sn1,Sn2,Sn3
	real*8, allocatable, dimension(:) :: Sn1m1,Sn2m1,Sn3m1
	real*8, allocatable, dimension(:) :: Sn1m2,Sn2m2,Sn3m2
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

subroutine update_E_inSi(this,E_pointer,D_pointer)
	class(Si) :: this
	real*8, pointer, intent(inout):: D_pointer(:)
        real*8,pointer, intent(inout) :: E_pointer(:)
	integer :: pos
	    if(this%posrange%init.gt.0.and.this%posrange%ends.gt.0)then
		do pos=this%posrange%init,this%posrange%ends,1
	E_pointer(pos)=D_pointer(pos)-this%Sn1(pos)-this%sn2(pos)-this%sn3(pos)
	this%Sn1(pos)=this%ga1(pos)*this%Sn1m1(pos)+this%gb1(pos)*this%Sn1m2(pos)+this%gc1(pos)*E_pointer(pos)
	this%Sn2(pos)=this%ga2(pos)*this%Sn2m1(pos)+this%gb2(pos)*this%Sn2m2(pos)+this%gc2(pos)*E_pointer(pos)
	this%Sn3(pos)=this%ga3(pos)*this%Sn3m1(pos)+this%gb3(pos)*this%Sn3m2(pos)+this%gc3(pos)*E_pointer(pos)
	this%Sn1m2(pos)=this%Sn1m1(pos); this%Sn2m2(pos)=this%Sn2m1(pos); this%Sn3m2(pos)=this%Sn3m1(pos)
	this%Sn1m1(pos)=this%Sn1(pos); this%Sn2m1(pos)=this%Sn2(pos); this%Sn3m1(pos)=this%Sn3(pos)
		enddo
	    endif
	end subroutine update_E_inSi

	subroutine initialize_Si(this,absrange,localdim)
	    class(Si),intent(inout) :: this
	    type(Irange), intent(in) :: absrange
	    integer, intent(in) :: localdim
	    call allocate_SiVars(this,localdim)
	    call this%material%initialize_material_indexes(absrange,localdim)
	    call initialize_SiVars(this)
        end subroutine initialize_Si

	subroutine allocate_SiVars(this,localdim)
	    class(Si) :: this
	    integer, intent(in) :: localdim
	    integer :: alloc_stat

	    allocate(this%ga1(localdim),stat=alloc_stat)
	    allocate(this%ga2(localdim),stat=alloc_stat)
	    allocate(this%ga3(localdim),stat=alloc_stat)

	    allocate(this%gb1(localdim),stat=alloc_stat)
	    allocate(this%gb2(localdim),stat=alloc_stat)
	    allocate(this%gb3(localdim),stat=alloc_stat)
	    allocate(this%gc1(localdim),stat=alloc_stat)
	    allocate(this%gc2(localdim),stat=alloc_stat)
	    allocate(this%gc3(localdim),stat=alloc_stat)

	    allocate(this%Sn1(localdim),stat=alloc_stat)
	    allocate(this%Sn1m1(localdim),stat=alloc_stat)
	    allocate(this%Sn1m2(localdim),stat=alloc_stat)
	    allocate(this%Sn2(localdim),stat=alloc_stat)
	    allocate(this%Sn2m1(localdim),stat=alloc_stat)
	    allocate(this%Sn2m2(localdim),stat=alloc_stat)
	    allocate(this%Sn3(localdim),stat=alloc_stat)
	    allocate(this%Sn3m1(localdim),stat=alloc_stat)
	    allocate(this%Sn3m2(localdim),stat=alloc_stat)

	    this%ga1=0.d0
	    this%ga2=0.d0
	    this%ga3=0.d0
	    this%gb1=0.d0
	    this%gb2=0.d0
	    this%gb3=0.d0
	    this%gc1=0.d0
	    this%gc2=0.d0
	    this%gc3=0.d0
	    this%Sn1=0.d0
	    this%Sn2=0.d0
	    this%Sn3=0.d0
	    this%Sn1m1=0.d0
	    this%Sn2m1=0.d0
	    this%Sn3m1=0.d0
	    this%Sn1m2=0.d0
	    this%Sn2m2=0.d0
	    this%Sn3m2=0.d0
	end subroutine allocate_SiVars

	subroutine initialize_SiVars(this)
	class(Si) :: this
	integer :: pos
	real*8 :: dt
	dt=this%dt
        do pos=this%posrange%init,this%posrange%ends,1
	    this%ga1(pos)=(2.d0-(dt*this%wp(1))**2)/(1.d0+dt*this%dnup(1))
	    this%ga2(pos)=(2.d0-(dt*this%wp(2))**2)/(1.d0+dt*this%dnup(2))
	    this%ga3(pos)=(2.d0-(dt*this%wp(3))**2)/(1.d0+dt*this%dnup(3))
	    this%gb1(pos)=(dt*this%dnup(1)-1.d0)/(1.d0+dt*this%dnup(1))
	    this%gb2(pos)=(dt*this%dnup(2)-1.d0)/(1.d0+dt*this%dnup(2))
	    this%gb3(pos)=(dt*this%dnup(3)-1.d0)/(1.d0+dt*this%dnup(3))
	    this%gc1(pos)=(this%dep(1)*(this%wp(1)*dt)**2)/(1.d0+dt*this%dnup(1))
	    this%gc2(pos)=(this%dep(2)*(this%wp(2)*dt)**2)/(1.d0+dt*this%dnup(2))
	    this%gc3(pos)=(this%dep(3)*(this%wp(3)*dt)**2)/(1.d0+dt*this%dnup(3))
	enddo
	end subroutine initialize_SiVars
	end module


