    module BERpmly_mod
    use constants_mod
    implicit none
    private
    public pmly,initialize,updateE_pml,updateH_pml
	type :: pmly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	PML public parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    integer :: npmlx,npmly,npmlz
	    integer :: xcord,ycord,zcord
	    integer :: N_x,N_y,N_z
	    integer :: Xmax,Ymax,Zmax
	    integer :: xBDlow,xBDhigh
	    integer :: yBDlow,yBDhigh
	    integer :: zBDlow,zBDhigh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	PML private variables 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		PML parameters
	    real*8, allocatable :: sigmaEx(:,:,:),sigmaEy(:,:,:),sigmaEz(:,:,:)
	    real*8, allocatable :: alphaEx(:,:,:),alphaEy(:,:,:),alphaEz(:,:,:)
	    real*8, allocatable :: f1x(:,:,:),f2x(:,:,:)
	    real*8, allocatable :: f1y(:,:,:),f2y(:,:,:)
	    real*8, allocatable :: f1z(:,:,:),f2z(:,:,:)

	    real*8, allocatable :: sigmaHx(:,:,:),sigmaHy(:,:,:),sigmaHz(:,:,:)
	    real*8, allocatable :: alphaHx(:,:,:),alphaHy(:,:,:),alphaHz(:,:,:)
	    real*8, allocatable :: g1x(:,:,:),g2x(:,:,:)
	    real*8, allocatable :: g1y(:,:,:),g2y(:,:,:)
	    real*8, allocatable :: g1z(:,:,:),g2z(:,:,:)
!		PML splitted fields
	    real*8, allocatable, public :: Exy(:,:,:),Exz(:,:,:)
	    real*8, allocatable, public  :: Eyz(:,:,:),Eyx(:,:,:)
	    real*8, allocatable, public  :: Ezx(:,:,:),Ezy(:,:,:)
	    real*8, allocatable, public  :: Hxy(:,:,:),Hxz(:,:,:)
	    real*8, allocatable, public  :: Hyz(:,:,:),Hyx(:,:,:)
	    real*8, allocatable, public  :: Hzx(:,:,:),Hzy(:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	PML procedures (subroutines or functions)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    contains
	    procedure :: updateE => updateE_pml
	    procedure :: setboundayEya
	    procedure :: setboundayEyb
	    procedure :: setPecEa
	    procedure :: setPecEb
	    procedure :: setPeriodicEa
	    procedure :: setPeriodicEb
	    procedure :: updateEs
	    procedure :: updateH => updateH_pml
	    procedure :: setboundayHya
	    procedure :: setboundayHyb
	    procedure :: setPecHa
	    procedure :: setPecHb
	    procedure :: setPeriodicHa
	    procedure :: setPeriodicHb
	    procedure :: updateHs
	    procedure :: initialize
	    procedure :: initialize_vars
	    procedure :: allocate_vars
	    procedure :: init_varx1
	    procedure :: init_varx2
	    procedure :: init_varz1
	    procedure :: init_varz2
	end type pmly
    contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	update fields calls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine updateE_pml(this,Ex_,Ey_,Ez_,Hx_,Hy_,Hz_)
	    class(pmly) :: this
	    real*8, pointer, intent(inout) :: Ex_(:,:,:),Ey_(:,:,:),Ez_(:,:,:)
	    real*8, pointer, intent(inout) :: Hx_(:,:,:),Hy_(:,:,:),Hz_(:,:,:)
	    integer :: xi,xf,yi,yf,zi,zf
	    if(this%ycord.eq.0.and.this%npmly.gt.0)then
		xi=1; xf=this%N_x; yi=1; yf=this%yBDlow; zi=1; zf=this%N_z
		call this%updateEs(xi,xf,yi,yf,zi,zf,Hx_,Hz_)
		call this%setboundayEya(Ex_,Ey_,Ez_)
!		call this%setPeriodicEa()
!		call this%setPecEa()
	    endif
	    if(this%ycord.eq.(this%Ymax-1).and.this%npmly.gt.0)then
		xi=1; xf=this%N_x; yi=this%yBDhigh; yf=this%N_y; zi=1; zf=this%N_z
		call this%updateEs(xi,xf,yi,yf,zi,zf,Hx_,Hz_)
		call this%setboundayEyb(Ex_,Ey_,Ez_)
!		call this%setPeriodicEb()
!		call this%setPecEb()
	    endif
	end subroutine updateE_pml

	subroutine updateH_pml(this,Ex_,Ey_,Ez_,Hx_,Hy_,Hz_)
	    class(pmly) :: this
	    real*8, pointer, intent(inout) :: Ex_(:,:,:),Ey_(:,:,:),Ez_(:,:,:)
	    real*8, pointer, intent(inout) :: Hx_(:,:,:),Hy_(:,:,:),Hz_(:,:,:)
	    integer :: xi,xf,yi,yf,zi,zf
	    if(this%ycord.eq.0.and.this%npmly.gt.0)then
		xi=1; xf=this%N_x; yi=1; yf=this%yBDlow; zi=1; zf=this%N_z
		call this%updateHs(xi,xf,yi,yf,zi,zf,Ex_,Ez_)
		call this%setboundayHya(Hx_,Hy_,Hz_)
!		call this%setPeriodicHa()
!		call this%setPecHa()
	    endif
	    if(this%ycord.eq.(this%Ymax-1).and.this%npmly.gt.0)then
		xi=1; xf=this%N_x; yi=this%yBDhigh; yf=this%N_y; zi=1; zf=this%N_z
		call this%updateHs(xi,xf,yi,yf,zi,zf,Ex_,Ez_)
		call this%setboundayHyb(Hx_,Hy_,Hz_)
!		call this%setPeriodicHb()
!		call this%setPecHb()
	    endif
	end subroutine updateH_pml

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Update equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine updateEs(this,xi,xf,yi,yf,zi,zf,Hx_,Hz_)
	    class(pmly) :: this
	    integer, intent(in) :: xi,xf,yi,yf,zi,zf
	    real*8, pointer, intent(inout) :: Hx_(:,:,:),Hz_(:,:,:)
	    real*8 :: curlHzy_p,curlHyz_p,curlHxz_p,curlHzx_p,curlHyx_p,curlHxy_p
	    real*8 :: curlHzy_m,curlHyz_m,curlHxz_m,curlHzx_m,curlHyx_m,curlHxy_m
	    integer :: i,j,ja,jb,k
	    ja=this%yBDlow
	    jb=this%yBDhigh
	    do k=zi,zf,1
		do j=yi,yf,1
		    do i=xi,xf,1
			curlHzy_p= this%Hzx(i,  j,  k)  +this%Hzy(i,  j,  k)
			curlHzy_m=-this%Hzx(i,  j-1,k)  -this%Hzy(i,  j-1,k)
			curlHyz_p=-this%Hyz(i,  j,  k)  -this%Hyx(i,  j,  k)
			curlHyz_m= this%Hyz(i,  j,  k-1)+this%Hyx(i,  j,  k-1)
			curlHxz_p= this%Hxy(i,  j,  k)  +this%Hxz(i,  j,  k)
			curlHxz_m=-this%Hxy(i,  j,  k-1)-this%Hxz(i,  j,  k-1)
			curlHzx_p=-this%Hzx(i,  j,  k)  -this%Hzy(i,  j,  k)
			curlHzx_m= this%Hzx(i-1,j,  k)  +this%Hzy(i-1,j,  k)
			curlHyx_p= this%Hyz(i,  j,  k)  +this%Hyx(i,  j,  k)
			curlHyx_m=-this%Hyz(i-1,j,  k)  -this%Hyx(i-1,j,  k)
			curlHxy_p=-this%Hxy(i,  j,  k)  -this%Hxz(i,  j,  k)
			curlHxy_m= this%Hxy(i,  j-1,k)  +this%Hxz(i,  j-1,k)

!			if(j.eq.ja.and.this%ycord.eq.0)then
!			    curlHzy_p=Hz_(i,j,k)
!			    curlHxy_p=-Hx_(i,j,k)
!			endif
!			if(j.eq.jb.and.this%ycord.eq.(this%Ymax-1))then
!			    curlHzy_m=-Hz_(i,j-1,k)
!			    curlHxy_m=Hx_(i,j-1,k)
!			endif
		this%Eyx(i,j,k)=this%f1x(i,j,k)*this%Eyx(i,j,k)+this%f2x(i,j,k)*(curlHzx_p+curlHzx_m)
		this%Ezx(i,j,k)=this%f1x(i,j,k)*this%Ezx(i,j,k)+this%f2x(i,j,k)*(curlHyx_p+curlHyx_m)
		this%Exy(i,j,k)=this%f1y(i,j,k)*this%Exy(i,j,k)+this%f2y(i,j,k)*(curlHzy_p+curlHzy_m)
		this%Ezy(i,j,k)=this%f1y(i,j,k)*this%Ezy(i,j,k)+this%f2y(i,j,k)*(curlHxy_p+curlHxy_m)
		this%Exz(i,j,k)=this%f1z(i,j,k)*this%Exz(i,j,k)+this%f2z(i,j,k)*(curlHyz_p+curlHyz_m)
		this%Eyz(i,j,k)=this%f1z(i,j,k)*this%Eyz(i,j,k)+this%f2z(i,j,k)*(curlHxz_p+curlHxz_m)
		    enddo
		enddo
	    enddo
	end subroutine updateEs

	subroutine updateHs(this,xi,xf,yi,yf,zi,zf,Ex_,Ez_)
	    class(pmly) :: this
	    integer, intent(in) :: xi,xf,yi,yf,zi,zf
	    real*8, pointer, intent(inout) :: Ex_(:,:,:),Ez_(:,:,:)
	    real*8 :: curlEzy_p,curlEyz_p,curlExz_p,curlEzx_p,curlEyx_p,curlExy_p
	    real*8 :: curlEzy_m,curlEyz_m,curlExz_m,curlEzx_m,curlEyx_m,curlExy_m
	    integer :: i,j,ja,jb,k
	    ja=this%yBDlow
	    jb=this%yBDhigh
	    do k=zi,zf,1
		do j=yi,yf,1
		    do i=xi,xf,1
			curlEzy_p=-this%Ezx(i,  j+1,k)  -this%Ezy(i,  j+1,k)
			curlEzy_m= this%Ezx(i,  j,  k)  +this%Ezy(i,  j,  k)
			curlEyz_p= this%Eyz(i,  j,  k+1)+this%Eyx(i,  j,  k+1)
			curlEyz_m=-this%Eyz(i,  j,  k)  -this%Eyx(i,  j,  k)
			curlExz_p=-this%Exy(i,  j,  k+1)-this%Exz(i,  j,  k+1)
			curlExz_m= this%Exy(i,  j,  k)  +this%Exz(i,  j,  k)
			curlEzx_p= this%Ezx(i+1,j,  k)  +this%Ezy(i+1,j,  k)
			curlEzx_m=-this%Ezx(i,  j,  k)  -this%Ezy(i  ,j,  k)
			curlEyx_p=-this%Eyz(i+1,j,  k)  -this%Eyx(i+1,j,  k)
			curlEyx_m= this%Eyz(i,  j,  k)  +this%Eyx(i,  j,  k)
			curlExy_p= this%Exy(i,  j+1,k)  +this%Exz(i,  j+1,k)
			curlExy_m=-this%Exy(i,  j,  k)  -this%Exz(i,  j,  k)

		this%Hyx(i,j,k)=this%g1x(i,j,k)*this%Hyx(i,j,k)+this%g2x(i,j,k)*(curlEzx_m+curlEzx_p)
		this%Hzx(i,j,k)=this%g1x(i,j,k)*this%Hzx(i,j,k)+this%g2x(i,j,k)*(curlEyx_m+curlEyx_p)
		this%Hxy(i,j,k)=this%g1y(i,j,k)*this%Hxy(i,j,k)+this%g2y(i,j,k)*(curlEzy_m+curlEzy_p)
		this%Hzy(i,j,k)=this%g1y(i,j,k)*this%Hzy(i,j,k)+this%g2y(i,j,k)*(curlExy_m+curlExy_p)
		this%Hxz(i,j,k)=this%g1z(i,j,k)*this%Hxz(i,j,k)+this%g2z(i,j,k)*(curlEyz_m+curlEyz_p)
		this%Hyz(i,j,k)=this%g1z(i,j,k)*this%Hyz(i,j,k)+this%g2z(i,j,k)*(curlExz_m+curlExz_p)
		    enddo
		enddo
	    enddo
	end subroutine updateHs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine setboundayEya(this,Ex_,Ey_,Ez_)
	    class(pmly) :: this
	    real*8, pointer, intent(inout) :: Ex_(:,:,:),Ey_(:,:,:),Ez_(:,:,:)
	    integer :: i,j,ja,k
	    ja=this%yBDlow

	    do k=1,this%N_z,1
		do i=1,this%N_x,1
		    this%Exy(i,ja+1,k)=Ex_(i,ja+1,k)
		    Ex_(i,ja,k)=this%Exy(i,ja,k)

		    this%Ezy(i,ja+1,k)=Ez_(i,ja+1,k)
		    Ez_(i,ja,k)=this%Ezy(i,ja,k)
		enddo
	    enddo

!	    do k=1,this%N_z,1
!		do j=1,ja-1,1
!		    do i=1,this%N_x,1
!			Ex_(i,j,k)=this%Exz(i,j,k)+this%Exy(i,j,k)
!			Ey_(i,j,k)=this%Eyz(i,j,k)+this%Eyx(i,j,k)
!			Ez_(i,j,k)=this%Ezy(i,j,k)+this%Ezx(i,j,k)
!		    enddo
!		enddo
!	    enddo
	    end subroutine setboundayEya

	subroutine setboundayHya(this,Hx_,Hy_,Hz_)
	    class(pmly) :: this
	    real*8, pointer, intent(inout) :: Hx_(:,:,:),Hy_(:,:,:),Hz_(:,:,:)
	    integer :: i,j,ja,k
	    ja=this%yBDlow

	    do k=1,this%N_z,1
		do i=1,this%N_x,1
		    this%Hxy(i,ja+1,k)=Hx_(i,ja+1,k)
		    Hx_(i,ja,k)=this%Hxy(i,ja,k)

		    this%Hzy(i,ja+1,k)=Hz_(i,ja+1,k)
		    Hz_(i,ja,k)=this%Hzy(i,ja,k)
		enddo
	    enddo

!	    do k=1,this%N_z,1
!		do j=1,ja-1,1
!		    do i=1,this%N_x,1
!			Hx_(i,j,k)=this%Hxy(i,j,k)+this%Hxz(i,j,k)
!			Hy_(i,j,k)=this%Hyz(i,j,k)+this%Hyx(i,j,k)
!			Hz_(i,j,k)=this%Hzy(i,j,k)+this%Hzx(i,j,k)
!		    enddo
!		enddo
!	    enddo

	end subroutine setboundayHya

	subroutine setboundayEyb(this,Ex_,Ey_,Ez_)
	    class(pmly) :: this
	    real*8, pointer, intent(inout) :: Ex_(:,:,:),Ey_(:,:,:),Ez_(:,:,:)
	    integer :: i,j,jb,k
	    jb=this%yBDhigh
	    do k=1,this%N_z,1
		do i=1,this%N_x,1
		    this%Exy(i,jb-1,k)=Ex_(i,jb-1,k)
		    Ex_(i,jb,k)=this%Exy(i,jb,k)

		    this%Ezy(i,jb-1,k)=Ez_(i,jb-1,k)
		    Ez_(i,jb,k)=this%Ezy(i,jb,k)
		enddo
	    enddo
!	    do k=1,this%N_z,1
!		do j=jb,this%N_y,1
!		    do i=1,this%N_x,1
!			Ex_(i,j,k)=this%Exy(i,j,k)+this%Exz(i,j,k)
!			Ey_(i,j,k)=this%Eyz(i,j,k)+this%Eyx(i,j,k)
!			Ez_(i,j,k)=this%Ezx(i,j,k)+this%Ezy(i,j,k)
!		    enddo
!		enddo
!	    enddo
	    end subroutine setboundayEyb

	subroutine setboundayHyb(this,Hx_,Hy_,Hz_)
	    class(pmly) :: this
	    real*8, pointer, intent(inout) :: Hx_(:,:,:),Hy_(:,:,:),Hz_(:,:,:)
	    real*8 :: curlEzy_p,curlEyz_p,curlExz_p,curlEzx_p,curlEyx_p,curlExy_p
	    real*8 :: curlEzy_m,curlEyz_m,curlExz_m,curlEzx_m,curlEyx_m,curlExy_m
	    integer :: i,j,jb,k
	    jb=this%yBDhigh

	    do k=1,this%N_z,1
		do i=1,this%N_x,1
		    this%Hxy(i,jb-1,k)=Hx_(i,jb-1,k)
		    Hx_(i,jb,k)=this%Hxy(i,jb,k)

		    this%Hzy(i,jb-1,k)=Hz_(i,jb-1,k)
		    Hz_(i,jb,k)=this%Hzy(i,jb,k)
		enddo
	    enddo
!	    do k=1,this%N_z,1
!		do j=jb,this%N_y
!		    do i=1,this%N_x,1
!			Hx_(i,j,k)=this%Hxy(i,j,k)+this%Hxz(i,j,k)
!			Hy_(i,j,k)=this%Hyz(i,j,k)+this%Hyx(i,j,k)
!			Hz_(i,j,k)=this%Hzy(i,j,k)+this%Hzx(i,j,k)
!		    enddo
!		enddo
!	    enddo
	end subroutine setboundayHyb

!set Pec boundary conditions
	subroutine setPecEa(this)
	class(pmly) :: this
	integer :: i,k
	    do k=0,this%N_z+1
		do i=0,this%N_x+1
		    this%Exy(i,0,k)=0.d0
		    this%Exz(i,0,k)=0.d0
		    this%Ezx(i,0,k)=0.d0
		    this%Ezy(i,0,k)=0.d0
		enddo
	    enddo
	end subroutine setPecEa

	subroutine setPecEb(this)
	class(pmly) :: this
	integer :: i,k
	    do k=0,this%N_z+1
		do i=0,this%N_x+1
		    this%Exy(i,this%N_y+1,k)=0.d0
		    this%Exz(i,this%N_y+1,k)=0.d0
		    this%Ezx(i,this%N_y+1,k)=0.d0
		    this%Ezy(i,this%N_y+1,k)=0.d0
		enddo
	    enddo
	end subroutine setPecEb

	subroutine setPecHa(this)
	class(pmly) :: this
	integer :: i,k
	    do k=0,this%N_z+1
		do i=0,this%N_x+1
		    this%Hyz(i,0,k)=0.d0
		    this%Hyx(i,0,k)=0.d0
		enddo
	    enddo

	end subroutine setPecHa

	subroutine setPecHb(this)
	class(pmly) :: this
	integer :: i,k
	    do k=0,this%N_z+1
		do i=0,this%N_x+1
		    this%Hyz(i,this%N_y+1,k)=0.d0
		    this%Hyx(i,this%N_y+1,k)=0.d0
		enddo
	    enddo
	end subroutine setPecHb

	subroutine setPeriodicEa(this)
	class(pmly) :: this
	integer :: i,j,k
	    do k=0,this%N_z+1
		do j=1,this%npmly
		    this%eyz(this%N_x+1,j,k)=this%eyz(1,j,k)
		    this%eyx(this%N_x+1,j,k)=this%eyx(1,j,k)
		    this%eyz(0,j,k)=this%eyz(this%N_x,j,k)
		    this%eyx(0,j,k)=this%eyx(this%N_x,j,k)

		    this%ezx(this%N_x+1,j,k)=this%ezx(1,j,k)
		    this%ezy(this%N_x+1,j,k)=this%ezy(1,j,k)
		    this%ezx(0,j,k)=this%ezx(this%N_x,j,k)
		    this%ezy(0,j,k)=this%ezy(this%N_x,j,k)
		enddo
	    enddo
	    do j=1,this%npmly
		do i=0,this%N_x+1
		    this%exy(i,j,this%N_z+1)=this%exy(i,j,1)
		    this%exz(i,j,this%N_z+1)=this%exz(i,j,1)
		    this%exy(i,j,0)=this%exy(i,j,this%N_z)
		    this%exz(i,j,0)=this%exz(i,j,this%N_z)

		    this%eyz(i,j,this%N_z+1)=this%eyz(i,j,1)
		    this%eyx(i,j,this%N_z+1)=this%eyx(i,j,1)
		    this%eyz(i,j,0)=this%eyz(i,j,this%N_z)
		    this%eyx(i,j,0)=this%eyx(i,j,this%N_z)
		enddo
	    enddo
	end subroutine setPeriodicEa

	subroutine setPeriodicEb(this)
	class(pmly) :: this
	integer :: i,j,k
	    do k=0,this%N_z+1
		do j=this%yBDhigh,this%N_y
		    this%eyz(this%N_x+1,j,k)=this%eyz(1,j,k)
		    this%eyx(this%N_x+1,j,k)=this%eyx(1,j,k)
		    this%eyz(0,j,k)=this%eyz(this%N_x,j,k)
		    this%eyx(0,j,k)=this%eyx(this%N_x,j,k)

		    this%ezx(this%N_x+1,j,k)=this%ezx(1,j,k)
		    this%ezy(this%N_x+1,j,k)=this%ezy(1,j,k)
		    this%ezx(0,j,k)=this%ezx(this%N_x,j,k)
		    this%ezy(0,j,k)=this%ezy(this%N_x,j,k)
		enddo
	    enddo
	    do j=this%yBDhigh,this%N_y
		do i=0,this%N_x+1
		    this%exy(i,j,this%N_z+1)=this%exy(i,j,1)
		    this%exz(i,j,this%N_z+1)=this%exz(i,j,1)
		    this%exy(i,j,0)=this%exy(i,j,this%N_z)
		    this%exz(i,j,0)=this%exz(i,j,this%N_z)

		    this%eyz(i,j,this%N_z+1)=this%eyz(i,j,1)
		    this%eyx(i,j,this%N_z+1)=this%eyx(i,j,1)
		    this%eyz(i,j,0)=this%eyz(i,j,this%N_z)
		    this%eyx(i,j,0)=this%eyx(i,j,this%N_z)
		enddo
	    enddo
	end subroutine setPeriodicEb

	subroutine setPeriodicHa(this)
	class(pmly) :: this
	integer :: i,j,k
	    do k=0,this%N_z+1
		do j=1,this%npmly
		    this%hyz(this%N_x+1,j,k)=this%hyz(1,j,k)
		    this%hyx(this%N_x+1,j,k)=this%hyx(1,j,k)
		    this%hyz(0,j,k)=this%hyz(this%N_x,j,k)
		    this%hyx(0,j,k)=this%hyx(this%N_x,j,k)

		    this%hzx(this%N_x+1,j,k)=this%hzx(1,j,k)
		    this%hzy(this%N_x+1,j,k)=this%hzy(1,j,k)
		    this%hzx(0,j,k)=this%hzx(this%N_x,j,k)
		    this%hzy(0,j,k)=this%hzy(this%N_x,j,k)
		enddo
	    enddo
	    do j=1,this%npmly
		do i=0,this%N_x+1
		    this%hxy(i,j,this%N_z+1)=this%hxy(i,j,1)
		    this%hxz(i,j,this%N_z+1)=this%hxz(i,j,1)
		    this%hxy(i,j,0)=this%hxy(i,j,this%N_z)
		    this%hxz(i,j,0)=this%hxz(i,j,this%N_z)

		    this%hyz(i,j,this%N_z+1)=this%hyz(i,j,1)
		    this%hyx(i,j,this%N_z+1)=this%hyx(i,j,1)
		    this%hyz(i,j,0)=this%hyz(i,j,this%N_z)
		    this%hyx(i,j,0)=this%hyx(i,j,this%N_z)
		enddo
	    enddo
	end subroutine setPeriodicHa

	subroutine setPeriodicHb(this)
	class(pmly) :: this
	integer :: i,j,k
	    do k=0,this%N_z+1
		do j=this%yBDhigh,this%N_y
		    this%hyz(this%N_x+1,j,k)=this%hyz(1,j,k)
		    this%hyx(this%N_x+1,j,k)=this%hyx(1,j,k)
		    this%hyz(0,j,k)=this%hyz(this%N_x,j,k)
		    this%hyx(0,j,k)=this%hyx(this%N_x,j,k)

		    this%hzx(this%N_x+1,j,k)=this%hzx(1,j,k)
		    this%hzy(this%N_x+1,j,k)=this%hzy(1,j,k)
		    this%hzx(0,j,k)=this%hzx(this%N_x,j,k)
		    this%hzy(0,j,k)=this%hzy(this%N_x,j,k)
		enddo
	    enddo
	    do j=this%yBDhigh,this%N_y
		do i=0,this%N_x+1
		    this%hxy(i,j,this%N_z+1)=this%hxy(i,j,1)
		    this%hxz(i,j,this%N_z+1)=this%hxz(i,j,1)
		    this%hxy(i,j,0)=this%hxy(i,j,this%N_z)
		    this%hxz(i,j,0)=this%hxz(i,j,this%N_z)

		    this%hyz(i,j,this%N_z+1)=this%hyz(i,j,1)
		    this%hyx(i,j,this%N_z+1)=this%hyx(i,j,1)
		    this%hyz(i,j,0)=this%hyz(i,j,this%N_z)
		    this%hyx(i,j,0)=this%hyx(i,j,this%N_z)
		enddo
	    enddo
	end subroutine setPeriodicHb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Initialization 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine initialize(this,eplow,ephigh,sigmam,npmls,order,dt,ds,coords,ranks,maxs)
	    class(pmly) :: this
	    integer, intent(in) :: npmls(3),coords(3),ranks(3),maxs(3)
	    real*8, intent(in) :: eplow,ephigh,order,dt,ds,sigmam
	    this%npmlx=npmls(1)
	    this%npmly=npmls(2)
	    this%npmlz=npmls(3)
	    this%xcord=coords(1)
	    this%ycord=coords(2)
	    this%zcord=coords(3)
	    this%N_x=ranks(1)
	    this%N_y=ranks(2)
	    this%N_z=ranks(3)
	    this%Xmax=maxs(1)
	    this%Ymax=maxs(2)
	    this%Zmax=maxs(3)

	    this%xBDlow=this%npmlx;this%xBDhigh=this%N_x-this%npmlx+1
	    this%yBDlow=this%npmly;this%yBDhigh=this%N_y-this%npmly+1
	    this%zBDlow=this%npmlz;this%zBDhigh=this%N_z-this%npmlz+1

	    call this%allocate_vars(eplow/ep0,ephigh/ep0,dt*c0/ds)
	    call this%initialize_vars(eplow,ephigh,sigmam,order,dt,ds)
	end subroutine initialize

	subroutine allocate_vars(this,eprlow,eprhigh,Ku)
	    class(pmly) :: this
	    real*8, intent(in) :: eprlow,eprhigh,Ku
	    real*8 :: epr
	    integer :: N_x,N_y,N_z
	    N_x=this%N_x
	    N_y=this%N_y
	    N_z=this%N_z

!    if((this%xcord.eq.0.or.this%xcord.eq.(this%Xmax-1)).and.&
!       (this%ycord.eq.0.or.this%ycord.eq.(this%Ymax-1)).and.&
!       (this%zcord.eq.0.or.this%zcord.eq.(this%Zmax-1)))then
    if(this%ycord.eq.0.or.this%ycord.eq.(this%Ymax-1))then
!	Along X
		allocate(this%sigmaEx(N_x,N_y,N_z))
		allocate(this%alphaEx(N_x,N_y,N_z))
		allocate(this%f1x(N_x,N_y,N_z))
		allocate(this%f2x(N_x,N_y,N_z))

		allocate(this%sigmaHx(N_x,N_y,N_z))
		allocate(this%alphaHx(N_x,N_y,N_z))
		allocate(this%g1x(N_x,N_y,N_z))
		allocate(this%g2x(N_x,N_y,N_z))
!	Along Y
		allocate(this%sigmaEy(N_x,N_y,N_z))
		allocate(this%alphaEy(N_x,N_y,N_z))
		allocate(this%f1y(N_x,N_y,N_z))
		allocate(this%f2y(N_x,N_y,N_z))

		allocate(this%sigmaHy(N_x,N_y,N_z))
		allocate(this%alphaHy(N_x,N_y,N_z))
		allocate(this%g1y(N_x,N_y,N_z))
		allocate(this%g2y(N_x,N_y,N_z))
!	Along Z
		allocate(this%sigmaEz(N_x,N_y,N_z))
		allocate(this%alphaEz(N_x,N_y,N_z))
		allocate(this%f1z(N_x,N_y,N_z))
		allocate(this%f2z(N_x,N_y,N_z))

		allocate(this%sigmaHz(N_x,N_y,N_z))
		allocate(this%alphaHz(N_x,N_y,N_z))
		allocate(this%g1z(N_x,N_y,N_z))
		allocate(this%g2z(N_x,N_y,N_z))

	    allocate(this%Exy(0:N_x+1,0:N_y+1,0:N_z+1))
	    allocate(this%Exz(0:N_x+1,0:N_y+1,0:N_z+1))
	    allocate(this%Eyz(0:N_x+1,0:N_y+1,0:N_z+1))
	    allocate(this%Eyx(0:N_x+1,0:N_y+1,0:N_z+1))
	    allocate(this%Ezx(0:N_x+1,0:N_y+1,0:N_z+1))
	    allocate(this%Ezy(0:N_x+1,0:N_y+1,0:N_z+1))
	    allocate(this%Hxy(0:N_x+1,0:N_y+1,0:N_z+1))
	    allocate(this%Hxz(0:N_x+1,0:N_y+1,0:N_z+1))
	    allocate(this%Hyz(0:N_x+1,0:N_y+1,0:N_z+1))
	    allocate(this%Hyx(0:N_x+1,0:N_y+1,0:N_z+1))
	    allocate(this%Hzx(0:N_x+1,0:N_y+1,0:N_z+1))
	    allocate(this%Hzy(0:N_x+1,0:N_y+1,0:N_z+1))

	    this%f2x=0.d0
	    this%f2y=0.d0!Ku/epr
	    this%f2z=0.d0

	    if(this%ycord.eq.0)then
!		this%f2x(:,1:this%yBDlow+1,:)=Ku/eprlow
!		this%f2y(:,1:this%yBDlow+1,:)=Ku/eprlow
!		this%f2z(:,1:this%yBDlow+1,:)=Ku/eprlow
		this%f2x=Ku/eprlow
		this%f2y=Ku/eprlow
		this%f2z=Ku/eprlow
	    endif
	    if(this%ycord.eq.(this%Ymax-1))then
!		this%f2x(:,this%yBDhigh-1,:)=Ku/eprhigh
!		this%f2y(:,this%yBDhigh-1,:)=Ku/eprhigh
!		this%f2z(:,this%yBDhigh-1,:)=Ku/eprhigh
		this%f2x=Ku/eprhigh
		this%f2y=Ku/eprhigh
		this%f2z=Ku/eprhigh
	    endif

	    this%sigmaEx=0.d0; this%alphaEx=0.d0
	    this%f1x=1.d0; 
	    this%sigmaHx=0.d0; this%alphaHx=0.d0
	    this%g1x=1.d0; this%g2x=Ku

	    this%sigmaEy=0.d0; this%alphaEy=0.d0
	    this%f1y=1.d0; 
	    this%sigmaHy=0.d0; this%alphaHy=0.d0
	    this%g1y=1.d0; this%g2y=Ku

	    this%sigmaEz=0.d0; this%alphaEz=0.d0
	    this%f1z=1.d0; 
	    this%sigmaHz=0.d0; this%alphaHz=0.d0
	    this%g1z=1.d0; this%g2z=Ku

	    this%Exy=0.d0; this%Exz=0.d0
	    this%Eyz=0.d0; this%Eyx=0.d0
	    this%Ezx=0.d0; this%Ezy=0.d0
	    this%Hxy=0.d0; this%Hxz=0.d0
	    this%Hyz=0.d0; this%Hyx=0.d0
	    this%Hzx=0.d0; this%Hzy=0.d0
	    endif
	end subroutine allocate_vars

	subroutine initialize_vars(this,eplow,ephigh,sigmam,order,dt,ds)
	    class(pmly) :: this
	    real*8, intent(in) :: eplow,ephigh,order,dt,ds,sigmam
	    integer :: i,j,k,jj
	    integer :: lastX,lastY,lastZ, npmlx, npmly, npmlz, nx, ny, nz
	    integer :: xBDlow,xBDhigh,yBDlow,yBDhigh,zBDlow,zBDhigh
	    real*8 :: z0
	    lastX=this%Xmax-1; lastY=this%Ymax-1; lastZ=this%Zmax-1
	    xBDlow=this%xBDlow; xBDhigh=this%xBDhigh
	    yBDlow=this%yBDlow; yBDhigh=this%yBDhigh
	    zBDlow=this%zBDlow; zBDhigh=this%zBDhigh
	    npmlx=this%npmlx;npmly=this%npmly;npmlz=this%npmlz
	    nx=this%N_x;ny=this%N_y; nz=this%N_z
	    z0=dsqrt(mu0/ep0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	pmlx-profile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	low profile
	if(this%xcord.eq.0.and.xBDlow.gt.0.and.this%ycord.eq.0)then
	    call this%init_varx1(1,npmlx,1,npmly,1,nz,eplow,sigmam,order,dt,ds)
	endif
	if(this%xcord.eq.0.and.xBDlow.gt.0.and.this%ycord.eq.lastY)then
	    call this%init_varx1(1,npmlx,yBDhigh,ny,1,nz,ephigh,sigmam,order,dt,ds)
	endif

!	if(this%xcord.eq.0.and.xBDlow.gt.0.and.this%zcord.eq.0.and.this%ycord.eq.0)then
!	    call this%init_varx1(1,npmlx,1,npmly,1,nz,eplow,sigmam,order,dt,ds)
!	endif
!	if(this%xcord.eq.0.and.xBDlow.gt.0.and.this%zcord.eq.lastZ.and.this%ycord.eq.0)then
!	    call this%init_varx1(1,npmlx,1,npmly,1,nz,eplow,sigmam,order,dt,ds)
!	endif
!	if(this%xcord.eq.0.and.xBDlow.gt.0.and.this%zcord.eq.0.and.this%ycord.eq.lastY)then
!	    call this%init_varx1(1,npmlx,yBDhigh,ny,1,nz,ephigh,sigmam,order,dt,ds)
!	endif
!	if(this%xcord.eq.0.and.xBDlow.gt.0.and.this%zcord.eq.lastZ.and.this%ycord.eq.lastY)then
!	    call this%init_varx1(1,npmlx,yBDhigh,ny,1,nz,ephigh,sigmam,order,dt,ds)
!	endif

!	high profile
	if(this%xcord.eq.lastX.and.xBDlow.gt.0.and.this%ycord.eq.0)then
	    call this%init_varx2(xBDhigh,nx,1,npmly,1,nz,eplow,sigmam,order,dt,ds)
	endif
	if(this%xcord.eq.lastX.and.xBDlow.gt.0.and.this%ycord.eq.lastY)then
	    call this%init_varx2(xBDhigh,nx,yBDhigh,ny,1,nz,ephigh,sigmam,order,dt,ds)
	endif

!	if(this%xcord.eq.lastX.and.xBDlow.gt.0.and.this%zcord.eq.0.and.this%ycord.eq.0)then
!	    call this%init_varx2(xBDhigh,nx,1,npmly,1,nz,eplow,sigmam,order,dt,ds)
!	endif
!	if(this%xcord.eq.lastX.and.xBDlow.gt.0.and.this%zcord.eq.lastZ.and.this%ycord.eq.0)then
!	    call this%init_varx2(xBDhigh,nx,1,npmly,1,nz,eplow,sigmam,order,dt,ds)
!	endif
!	if(this%xcord.eq.lastX.and.xBDlow.gt.0.and.this%zcord.eq.0.and.this%ycord.eq.lastY)then
!	    call this%init_varx2(xBDhigh,nx,yBDhigh,ny,1,nz,ephigh,sigmam,order,dt,ds)
!	endif
!	if(this%xcord.eq.lastX.and.xBDlow.gt.0.and.this%zcord.eq.lastZ.and.this%ycord.eq.lastY)then
!	    call this%init_varx2(xBDhigh,nx,yBDhigh,ny,1,nz,ephigh,sigmam,order,dt,ds)
!	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	pmly-profile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	low profile
	if(this%ycord.eq.0.and.yBDlow.gt.0)then
	    do k=1,nz,1
		do j=1,yBDlow,1
		    do i=1,nx,1
			this%sigmaEy(i,j,k)=eplow*sigmam*((real(npmly-j,8)+0.5d0)/real(npmly,8))**order
			this%alphaEy(i,j,k)=this%sigmaEy(i,j,k)*dt/eplow
			    this%f1y(i,j,k)=dexp(-this%alphaEy(i,j,k))
			    this%f2y(i,j,k)=(1.d0-this%f1y(i,j,k))/(z0*this%sigmaEy(i,j,k)*ds)

			if(j.lt.npmly)then
			    this%sigmaHy(i,j,k)=mu0*sigmam*((real(npmly-j,8))/real(npmly,8))**order
			    this%alphaHy(i,j,k)=this%sigmaHy(i,j,k)*dt/mu0
			        this%g1y(i,j,k)=dexp(-this%alphaHy(i,j,k))
			        this%g2y(i,j,k)=z0*(1.d0-this%g1y(i,j,k))/(this%sigmaHy(i,j,k)*ds)
			endif
		    enddo
		enddo
	    enddo
	endif

!	high profile
	if(this%ycord.eq.lastY.and.yBDlow.gt.0)then
	    do k=1,nz,1
		do i=1,nx,1
		    jj=npmly
		    do j=yBDhigh,ny,1
			if(jj.lt.npmly)then
			    this%sigmaEy(i,j,k)=ephigh*sigmam*((real(npmly-jj,8))/real(npmly,8))**order
			    this%alphaEy(i,j,k)=this%sigmaEy(i,j,k)*dt/ephigh
				this%f1y(i,j,k)=dexp(-this%alphaEy(i,j,k))
				this%f2y(i,j,k)=(1.d0-this%f1y(i,j,k))/(z0*this%sigmaEy(i,j,k)*ds)
			endif
			this%sigmaHy(i,j,k)=mu0*sigmam*((real(npmly-jj,8)+0.5d0)/real(npmly,8))**order
			this%alphaHy(i,j,k)=this%sigmaHy(i,j,k)*dt/mu0
			this%g1y(i,j,k)=dexp(-this%alphaHy(i,j,k))
			this%g2y(i,j,k)=z0*(1.d0-this%g1y(i,j,k))/(this%sigmaHy(i,j,k)*ds)
			jj=jj-1
		    enddo
		enddo
	    enddo
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	pmlz-profile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	low profile
	if(this%zcord.eq.0.and.zBDlow.gt.0.and.this%ycord.eq.0)then
	    call this%init_varz1(1,nx,1,npmly,1,npmlz,eplow,sigmam,order,dt,ds)
	endif
	if(this%zcord.eq.0.and.zBDlow.gt.0.and.this%ycord.eq.lastY)then
	    call this%init_varz1(1,nx,yBDhigh,ny,1,npmlz,ephigh,sigmam,order,dt,ds)
	endif

!	if(this%zcord.eq.0.and.zBDlow.gt.0.and.this%xcord.eq.0.and.this%ycord.eq.0)then
!	    call this%init_varz1(1,nx,1,npmly,1,npmlz,eplow,sigmam,order,dt,ds)
!	endif
!	if(this%zcord.eq.0.and.zBDlow.gt.0.and.this%xcord.eq.lastX.and.this%ycord.eq.0)then
!	    call this%init_varz1(1,nx,1,npmly,1,npmlz,eplow,sigmam,order,dt,ds)
!	endif
!	if(this%zcord.eq.0.and.zBDlow.gt.0.and.this%xcord.eq.0.and.this%ycord.eq.lastY)then
!	    call this%init_varz1(1,nx,yBDhigh,ny,1,npmlz,ephigh,sigmam,order,dt,ds)
!	endif
!	if(this%zcord.eq.0.and.zBDlow.gt.0.and.this%xcord.eq.lastX.and.this%ycord.eq.lastY)then
!	    call this%init_varz1(1,nx,yBDhigh,ny,1,npmlz,ephigh,sigmam,order,dt,ds)
!	endif

!	high profile
	if(this%zcord.eq.lastZ.and.zBDlow.gt.0.and.this%ycord.eq.0)then
	    call this%init_varz2(1,nx,1,npmly,zBDhigh,nz,eplow,sigmam,order,dt,ds)
	endif
	if(this%zcord.eq.lastZ.and.zBDlow.gt.0.and.this%ycord.eq.lastY)then
	    call this%init_varz2(1,nx,yBDhigh,ny,zBDhigh,nz,ephigh,sigmam,order,dt,ds)
	endif

!	if(this%zcord.eq.lastZ.and.zBDlow.gt.0.and.this%xcord.eq.0.and.this%ycord.eq.0)then
!	    call this%init_varz2(1,nx,1,npmly,zBDhigh,nz,eplow,sigmam,order,dt,ds)
!	endif
!	if(this%zcord.eq.lastZ.and.zBDlow.gt.0.and.this%xcord.eq.lastX.and.this%ycord.eq.0)then
!	    call this%init_varz2(1,nx,1,npmly,zBDhigh,nz,eplow,sigmam,order,dt,ds)
!	endif
!	if(this%zcord.eq.lastZ.and.zBDlow.gt.0.and.this%xcord.eq.0.and.this%ycord.eq.lastY)then
!	    call this%init_varz2(1,nx,yBDhigh,ny,zBDhigh,nz,ephigh,sigmam,order,dt,ds)
!	endif
!	if(this%zcord.eq.lastZ.and.zBDlow.gt.0.and.this%xcord.eq.lastX.and.this%ycord.eq.lastY)then
!	    call this%init_varz2(1,nx,yBDhigh,ny,zBDhigh,nz,ephigh,sigmam,order,dt,ds)
!	endif
    end subroutine initialize_vars

    subroutine init_varx1(this,xi,xf,yi,yf,zi,zf,ep,sigmam,order,dt,ds)
    class(pmly) :: this
    integer, intent(in) :: xi,xf,yi,yf,zi,zf
    real*8, intent(in) :: ep,sigmam,order,dt,ds
    integer :: i,j,k,npmlx
    real*8 :: z0
    npmlx=this%npmlx
    z0=dsqrt(mu0/ep0)
	do k=zi,zf,1
	    do j=yi,yf,1
		do i=xi,xf,1
		    this%sigmaEx(i,j,k)=ep*sigmam*((real(npmlx-i,8)+0.5d0)/real(npmlx,8))**order
		    this%alphaEx(i,j,k)=this%sigmaEx(i,j,k)*dt/ep
			this%f1x(i,j,k)=dexp(-this%alphaEx(i,j,k))
			this%f2x(i,j,k)=(1.d0-this%f1x(i,j,k))/(z0*this%sigmaEx(i,j,k)*ds)

		    if(i.lt.npmlx)then
			this%sigmaHx(i,j,k)=mu0*sigmam*((real(npmlx-i,8))/real(npmlx,8))**order
			this%alphaHx(i,j,k)=this%sigmaHx(i,j,k)*dt/mu0
			    this%g1x(i,j,k)=dexp(-this%alphaHx(i,j,k))
			    this%g2x(i,j,k)=z0*(1.d0-this%g1x(i,j,k))/(this%sigmaHx(i,j,k)*ds)
		    endif
		enddo
	    enddo
	enddo
    end subroutine init_varx1

    subroutine init_varx2(this,xi,xf,yi,yf,zi,zf,ep,sigmam,order,dt,ds)
    class(pmly) :: this
    integer, intent(in) :: xi,xf,yi,yf,zi,zf
    real*8, intent(in) :: ep,sigmam,order,dt,ds
    integer :: i,ii,j,k,npmlx
    real*8 :: z0
    npmlx=this%npmlx
    z0=dsqrt(mu0/ep0)
	do k=zi,zf,1
	    do j=yi,yf,1
		ii=npmlx
		do i=xi,xf,1
		    if(ii.lt.npmlx)then
			this%sigmaEx(i,j,k)=ep*sigmam*((real(npmlx-ii,8))/real(npmlx,8))**order
			this%alphaEx(i,j,k)=this%sigmaEx(i,j,k)*dt/ep
			    this%f1x(i,j,k)=dexp(-this%alphaEx(i,j,k))
			    this%f2x(i,j,k)=(1.d0-this%f1x(i,j,k))/(z0*this%sigmaEx(i,j,k)*ds)
		    endif

		    this%sigmaHx(i,j,k)=mu0*sigmam*((real(npmlx-ii,8)+0.5d0)/real(npmlx,8))**order
		    this%alphaHx(i,j,k)=this%sigmaHx(i,j,k)*dt/mu0
		        this%g1x(i,j,k)=dexp(-this%alphaHx(i,j,k))
		        this%g2x(i,j,k)=z0*(1.d0-this%g1x(i,j,k))/(this%sigmaHx(i,j,k)*ds)
		    ii=ii-1
		enddo
	    enddo
	enddo
    end subroutine init_varx2

    subroutine init_varz1(this,xi,xf,yi,yf,zi,zf,ep,sigmam,order,dt,ds)
    class(pmly) :: this
    integer, intent(in) :: xi,xf,yi,yf,zi,zf
    real*8, intent(in) :: ep,sigmam,order,dt,ds
    integer :: i,j,k,npmlz
    real*8 :: z0
    npmlz=this%npmlz
    z0=dsqrt(mu0/ep0)
	do k=zi,zf,1
	    do j=yi,yf,1
		do i=xi,xf,1
		    this%sigmaEz(i,j,k)=ep*sigmam*((real(npmlz-k,8)+0.5d0)/real(npmlz,8))**order
		    this%alphaEz(i,j,k)=this%sigmaEz(i,j,k)*dt/ep
			this%f1z(i,j,k)=dexp(-this%alphaEz(i,j,k))
			this%f2z(i,j,k)=(1.d0-this%f1z(i,j,k))/(z0*this%sigmaEz(i,j,k)*ds)

		    if(k.lt.npmlz)then
			this%sigmaHz(i,j,k)=mu0*sigmam*((real(npmlz-k,8))/real(npmlz,8))**order
			this%alphaHz(i,j,k)=this%sigmaHz(i,j,k)*dt/mu0
			    this%g1z(i,j,k)=dexp(-this%alphaHz(i,j,k))
			    this%g2z(i,j,k)=z0*(1.d0-this%g1z(i,j,k))/(this%sigmaHz(i,j,k)*ds)
		    endif
		enddo
	    enddo
	enddo

    end subroutine init_varz1

    subroutine init_varz2(this,xi,xf,yi,yf,zi,zf,ep,sigmam,order,dt,ds)
    class(pmly) :: this
    integer, intent(in) :: xi,xf,yi,yf,zi,zf
    real*8, intent(in) :: ep,sigmam,order,dt,ds
    integer :: i,j,k,kk,npmlz
    real*8 :: z0
    npmlz=this%npmlz
    z0=dsqrt(mu0/ep0)

	do i=xi,xf,1
	    do j=yi,yf,1
		kk=npmlz
		do k=zi,zf,1
		    if(kk.lt.npmlz)then
			this%sigmaEz(i,j,k)=ep*sigmam*((real(npmlz-kk,8))/real(npmlz,8))**order
			this%alphaEz(i,j,k)=this%sigmaEz(i,j,k)*dt/ep
			    this%f1z(i,j,k)=dexp(-this%alphaEz(i,j,k))
			    this%f2z(i,j,k)=(1.d0-this%f1z(i,j,k))/(z0*this%sigmaEz(i,j,k)*ds)
		    endif

		    this%sigmaHz(i,j,k)=mu0*sigmam*((real(npmlz-kk,8)+0.5d0)/real(npmlz,8))**order
		    this%alphaHz(i,j,k)=this%sigmaHz(i,j,k)*dt/mu0
		        this%g1z(i,j,k)=dexp(-this%alphaHz(i,j,k))
		        this%g2z(i,j,k)=z0*(1.d0-this%g1z(i,j,k))/(this%sigmaHz(i,j,k)*ds)
		    kk=kk-1
		enddo
	    enddo
	enddo

    end subroutine init_varz2

    end module
