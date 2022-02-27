    module BERpmly_HzExEy_mod
    use constants_mod
    implicit none
    private
    public pmly,initialize,updateE_pml,updateH_pml
	type :: pmly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	PML public parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    integer :: npmlx,npmly
	    integer :: xcord,ycord
	    integer :: N_x,N_y
	    integer :: Xmax,Ymax
	    integer :: xBDlow,xBDhigh
	    integer :: yBDlow,yBDhigh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	PML private variables 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		PML parameters
	    real*8, allocatable :: sigmaEx(:,:),sigmaEy(:,:)
	    real*8, allocatable :: alphaEx(:,:),alphaEy(:,:)
	    real*8, allocatable :: f1x(:,:),f2x(:,:)
	    real*8, allocatable :: f1y(:,:),f2y(:,:)

	    real*8, allocatable :: sigmaHx(:,:),sigmaHy(:,:)
	    real*8, allocatable :: alphaHx(:,:),alphaHy(:,:)
	    real*8, allocatable :: g1x(:,:),g2x(:,:)
	    real*8, allocatable :: g1y(:,:),g2y(:,:)
!		PML splitted fields
	    real*8, allocatable :: Ex(:,:),Ey(:,:)
	    real*8, allocatable :: Hzx(:,:),Hzy(:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	PML procedures (subroutines or functions)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    contains
	    procedure :: updateE => updateE_pml
	    procedure :: setboundayEya
	    procedure :: setboundayEyb
	    procedure :: setPecEa
	    procedure :: setPecEb
	    procedure :: updateEs
	    procedure :: updateH => updateH_pml
	    procedure :: setboundayHya
	    procedure :: setboundayHyb
	    procedure :: setPecHa
	    procedure :: setPecHb
	    procedure :: updateHs
	    procedure :: initialize
	    procedure :: initialize_vars
	    procedure :: init_varx1
	    procedure :: init_varx2
	    procedure :: allocate_vars
	end type pmly
    contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	update fild calls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine updateE_pml(this,Ex_,Hz_)
	    class(pmly) :: this
	    real*8, pointer, intent(inout) :: Ex_(:,:),Hz_(:,:)
	    integer :: xi,xf,yi,yf
	    if(this%ycord.eq.0.and.this%npmly.gt.0)then
		xi=1; xf=this%N_x; yi=1; yf=this%yBDlow
		call this%updateEs(xi,xf,yi,yf,Hz_)
		call this%setboundayEya(Ex_)
		call this%setPecEa()
	    endif
	    if(this%ycord.eq.(this%Ymax-1).and.this%npmly.gt.0)then
		xi=1; xf=this%N_x; yi=this%yBDhigh; yf=this%N_y
		call this%updateEs(xi,xf,yi,yf,Hz_)
		call this%setboundayEyb(Ex_)
		call this%setPecEb()
	    endif
	end subroutine updateE_pml

	subroutine updateH_pml(this,Hz_)
	    class(pmly) :: this
	    real*8, pointer, intent(inout) :: Hz_(:,:)
	    integer :: xi,xf,yi,yf
	    if(this%ycord.eq.0.and.this%npmly.gt.0)then
		xi=1; xf=this%N_x; yi=1; yf=this%yBDlow
		call this%updateHs(xi,xf,yi,yf)
		call this%setboundayHya(Hz_)
		call this%setPecHa()
	    endif
	    if(this%ycord.eq.(this%Ymax-1).and.this%npmly.gt.0)then
		xi=1; xf=this%N_x; yi=this%yBDhigh; yf=this%N_y
		call this%updateHs(xi,xf,yi,yf)
		call this%setboundayHyb(Hz_)
		call this%setPecHb()
	    endif
	end subroutine updateH_pml
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	update equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine updateEs(this,xi,xf,yi,yf,Hz_)
	    class(pmly) :: this
	    integer, intent(in) :: xi,xf,yi,yf
	    real*8, pointer, intent(inout) :: Hz_(:,:)
	    real*8 :: curlHzx_p,curlHzy_p
	    real*8 :: curlHzx_m,curlHzy_m
	    integer :: i,j,ja,jb
	    ja=this%yBDlow
	    jb=this%yBDhigh
	    do j=yi,yf,1
	        do i=xi,xf,1
		    curlHzx_p=-this%Hzx(i,  j)  -this%Hzy(i,  j)
		    curlHzx_m= this%Hzx(i-1,j)  +this%Hzy(i-1,j)
		    curlHzy_p= this%Hzx(i,  j)  +this%Hzy(i,  j)
		    curlHzy_m=-this%Hzx(i,  j-1)-this%Hzy(i,  j-1)

		    if(j.eq.ja)then
			curlHzy_p=Hz_(i,j)
		    endif
		    if(j.eq.jb)then
			curlHzy_m=-Hz_(i,j-1)
		    endif

		    this%Ex(i,j)=this%f1y(i,j)*this%Ex(i,j)+this%f2y(i,j)*(curlHzy_p+curlHzy_m)
		    this%Ey(i,j)=this%f1x(i,j)*this%Ey(i,j)+this%f2x(i,j)*(curlHzx_p+curlHzx_m)
		enddo
	    enddo
	end subroutine updateEs

	subroutine updateHs(this,xi,xf,yi,yf)
	    class(pmly) :: this
	    integer, intent(in) :: xi,xf,yi,yf
	    real*8 :: curlEyx,curlExy
	    integer :: i,j,ja,ia
	    ja=this%yBDlow
	    do j=yi,yf,1
	        do i=xi,xf,1
		    curlExy= this%Ex(i  ,j+1)-this%Ex(i,j)
		    curlEyx=-this%Ey(i+1,j)  +this%Ey(i,j)

		    this%Hzy(i,j)=this%g1y(i,j)*this%Hzy(i,j)+this%g2y(i,j)*curlExy
		    this%Hzx(i,j)=this%g1x(i,j)*this%Hzx(i,j)+this%g2x(i,j)*curlEyx
		enddo
	    enddo
	end subroutine updateHs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine setboundayEya(this,Ex_)
	    class(pmly) :: this
	    real*8, pointer, intent(inout) :: Ex_(:,:)
	    integer :: i,j,ja
	    ja=this%yBDlow
	    do j=1,ja,1
	        do i=0,this%N_x+1,1
			Ex_(i,j)=this%Ex(i,j)
		enddo
	    enddo
	    end subroutine setboundayEya

	subroutine setboundayHya(this,Hz_)
	    class(pmly) :: this
	    real*8, pointer, intent(inout) :: Hz_(:,:)
	    integer :: i,j,ja
	    ja=this%yBDlow
	    do j=1,ja-1,1
		do i=0,this%N_x+1,1
		    Hz_(i,j)=this%Hzy(i,j)+this%Hzx(i,j)
		enddo
	    enddo
	end subroutine setboundayHya

	subroutine setboundayEyb(this,Ex_)
	    class(pmly) :: this
	    real*8, pointer, intent(inout) :: Ex_(:,:)
	    integer :: i,j,jb
	    jb=this%yBDhigh
	    do j=jb,this%N_y,1
		do i=0,this%N_x+1,1
		    Ex_(i,j)=this%Ex(i,j)
		enddo
	    enddo
	    end subroutine setboundayEyb

	subroutine setboundayHyb(this,Hz_)
	    class(pmly) :: this
	    real*8, pointer, intent(inout) :: Hz_(:,:)
	    integer :: i,j,jb
	    jb=this%yBDhigh
	    do j=jb,this%N_y
		do i=0,this%N_x+1,1
		    Hz_(i,j)=this%Hzy(i,j)+this%Hzx(i,j)
		enddo
	    enddo
	end subroutine setboundayHyb

!	set Pec boundary conditions
	subroutine setPecEa(this)
	class(pmly) :: this
	integer :: i
	do i=0,this%N_x+1
	    this%Ex(i,0)=0.d0
	enddo
	end subroutine setPecEa

	subroutine setPecEb(this)
	class(pmly) :: this
	integer :: i
	do i=0,this%N_x+1
	    this%Ex(i,this%N_y+1)=0.d0
	enddo
	end subroutine setPecEb

	subroutine setPecHa(this)
	class(pmly) :: this
	integer :: i
	do i=0,this%N_x+1
	    this%Hzx(i,0)=0.d0
	    this%Hzy(i,0)=0.d0
	enddo
	end subroutine setPecHa

	subroutine setPecHb(this)
	class(pmly) :: this
	integer :: i
	do i=0,this%N_x+1
	    this%Hzx(i,this%N_y+1)=0.d0
	    this%Hzy(i,this%N_y+1)=0.d0
	enddo
	end subroutine setPecHb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Initialization 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine initialize(this,eplow,ephigh,sigmam,npmls,order,dt,ds,coords,ranks,maxs)
	    class(pmly) :: this
	    integer, intent(in) :: npmls(2),coords(2),ranks(2),maxs(2)
	    real*8, intent(in) :: eplow,ephigh,order,dt,ds,sigmam
	    this%npmlx=npmls(1)
	    this%npmly=npmls(2)
	    this%xcord=coords(1)
	    this%ycord=coords(2)
	    this%N_x=ranks(1)
	    this%N_y=ranks(2)
	    this%Xmax=maxs(1)
	    this%Ymax=maxs(2)

	    this%xBDlow=this%npmlx;this%xBDhigh=this%N_x-this%npmlx+1
	    this%yBDlow=this%npmly;this%yBDhigh=this%N_y-this%npmly+1

	    call this%allocate_vars(eplow/ep0,ephigh/ep0,dt*c0/ds)
	    call this%initialize_vars(eplow,ephigh,sigmam,order,dt,ds)
	end subroutine initialize

	subroutine allocate_vars(this,eprlow,eprhigh,Ku)
	    class(pmly) :: this
	    real*8, intent(in) :: eprlow,eprhigh,Ku
	    real*8 :: epr
	    integer :: N_x,N_y
	    N_x=this%N_x
	    N_y=this%N_y
    if((this%xcord.eq.0.or.this%xcord.eq.(this%Xmax-1)).and.&
       (this%ycord.eq.0.or.this%ycord.eq.(this%Ymax-1)))then
!	Along X
		allocate(this%sigmaEx(N_x,N_y))
		allocate(this%alphaEx(N_x,N_y))
		allocate(this%f1x(N_x,N_y))
		allocate(this%f2x(N_x,N_y))

		allocate(this%sigmaHx(N_x,N_y))
		allocate(this%alphaHx(N_x,N_y))
		allocate(this%g1x(N_x,N_y))
		allocate(this%g2x(N_x,N_y))
!	Along Y
		allocate(this%sigmaEy(N_x,N_y))
		allocate(this%alphaEy(N_x,N_y))
		allocate(this%f1y(N_x,N_y))
		allocate(this%f2y(N_x,N_y))

		allocate(this%sigmaHy(N_x,N_y))
		allocate(this%alphaHy(N_x,N_y))
		allocate(this%g1y(N_x,N_y))
		allocate(this%g2y(N_x,N_y))

	    allocate(this%Ex(0:N_x+1,0:N_y+1))
	    allocate(this%Ey(0:N_x+1,0:N_y+1))
	    allocate(this%Hzx(0:N_x+1,0:N_y+1))
	    allocate(this%Hzy(0:N_x+1,0:N_y+1))

	    this%f2x=0.d0
	    this%f2y=0.d0
	    if(this%ycord.eq.0)then
!		this%f2x(:,1:this%yBDlow+5)=Ku/eprlow
!		this%f2y(:,1:this%yBDlow+5)=Ku/eprlow
		this%f2x(:,1:this%N_y/2-1)=Ku/eprlow
		this%f2y(:,1:this%N_y/2-1)=Ku/eprlow
	    endif
	    if(this%ycord.eq.(this%Ymax-1))then
		this%f2x(:,this%N_y/2:)=Ku/eprhigh
		this%f2y(:,this%N_y/2:)=Ku/eprhigh
!		this%f2x(:,this%yBDhigh-5:)=Ku/eprhigh
!		this%f2y(:,this%yBDhigh-5:)=Ku/eprhigh
	    endif

	    this%sigmaEx=0.d0; this%alphaEx=0.d0
	    this%f1x=1.d0; 
	    this%sigmaHx=0.d0; this%alphaHx=0.d0
	    this%g1x=1.d0; this%g2x=Ku

	    this%sigmaEy=0.d0; this%alphaEy=0.d0
	    this%f1y=1.d0;
	    this%sigmaHy=0.d0; this%alphaHy=0.d0
	    this%g1y=1.d0; this%g2y=Ku

	    this%Ex=0.d0; this%Ey=0.d0
	    this%Hzx=0.d0; this%Hzy=0.d0
	    endif
	end subroutine allocate_vars

	subroutine initialize_vars(this,eplow,ephigh,sigmam,order,dt,ds)
	    class(pmly) :: this
	    integer :: npml
	    real*8, intent(in) :: eplow,ephigh,order,dt,ds,sigmam
	    integer :: i,j,ii,jj
	    integer :: lastX,lastY,npmlx,npmly,nx,ny
	    integer :: xBDlow,xBDhigh,yBDlow,yBDhigh
	    real*8 :: z0,sigma0
	    lastX=this%Xmax-1; lastY=this%Ymax-1
	    xBDlow=this%xBDlow; yBDlow=this%yBDlow
	    xBDhigh=this%xBDhigh; yBDhigh=this%yBDhigh
	    npmlx=this%npmlx; npmly=this%npmly
	    nx=this%N_x; ny=this%N_y
	    z0=dsqrt(mu0/ep0)
	if(this%xcord.eq.0.and.xBDlow.gt.0.and.this%ycord.eq.0)then
	    call this%init_varx1(1,npmlx,1,npmly,eplow,sigmam,order,dt,ds)
	endif
	if(this%xcord.eq.0.and.xBDlow.gt.0.and.this%ycord.eq.lastY)then
	    call this%init_varx1(1,npmlx,yBDhigh,ny,ephigh,sigmam,order,dt,ds)
	endif
	if(this%xcord.eq.lastX.and.xBDlow.gt.0.and.this%ycord.eq.0)then
	    call this%init_varx2(this%xBDhigh,nx,1,npmly,eplow,sigmam,order,dt,ds)
	endif
	if(this%xcord.eq.lastX.and.xBDlow.gt.0.and.this%ycord.eq.lastY)then
	    call this%init_varx2(this%xBDhigh,nx,yBDhigh,ny,ephigh,sigmam,order,dt,ds)
	endif
	if(this%ycord.eq.0.and.yBDlow.gt.0)then
	    do j=1,yBDlow,1
		do i=1,nx,1
		    this%sigmaEy(i,j)=eplow*sigmam*((real(npmly-j,8)+0.5d0)/real(npmly,8))**order
		    this%alphaEy(i,j)=this%sigmaEy(i,j)*dt/eplow
		    this%f1y(i,j)=dexp(-this%alphaEy(i,j))
		    this%f2y(i,j)=(1.d0-this%f1y(i,j))/(z0*this%sigmaEy(i,j)*ds)
		    if(j.lt.npmly)then
			this%sigmaHy(i,j)=mu0*sigmam*((real(npmly-j,8))/real(npmly,8))**order
			this%alphaHy(i,j)=this%sigmaHy(i,j)*dt/mu0
			this%g1y(i,j)=dexp(-this%alphaHy(i,j))
			this%g2y(i,j)=z0*(1.d0-this%g1y(i,j))/(this%sigmaHy(i,j)*ds)
		    endif
		enddo
	    enddo
	endif
	if(this%ycord.eq.lastY.and.yBDlow.gt.0)then
	    do i=1,nx,1
		jj=npmly
		do j=yBDhigh,ny,1
		    if(jj.lt.npmly)then
			this%sigmaEy(i,j)=ephigh*sigmam*((real(npmly-jj,8))/real(npmly,8))**order
			this%alphaEy(i,j)=this%sigmaEy(i,j)*dt/ephigh
			this%f1y(i,j)=dexp(-this%alphaEy(i,j))
			this%f2y(i,j)=(1.d0-this%f1y(i,j))/(z0*this%sigmaEy(i,j)*ds)
		    endif

			this%sigmaHy(i,j)=mu0*sigmam*((real(npmly-jj,8)+0.5d0)/real(npmly,8))**order
			this%alphaHy(i,j)=this%sigmaHy(i,j)*dt/mu0
			this%g1y(i,j)=dexp(-this%alphaHy(i,j))
			this%g2y(i,j)=z0*(1.d0-this%g1y(i,j))/(this%sigmaHy(i,j)*ds)
		    jj=jj-1
		enddo
	    enddo
	endif
	end subroutine initialize_vars

	subroutine init_varx1(this,xi,xf,yi,yf,ep,sigmam,order,dt,ds)
	class(pmly) :: this
	integer, intent(in) :: xi,xf,yi,yf
	real*8, intent(in) :: ep,sigmam,order,dt,ds
	integer :: i,j,npmlx
	real*8 :: z0
	npmlx=this%npmlx
	z0=dsqrt(mu0/ep0)
	do j=yi,yf,1
	    do i=xi,xf,1
		this%sigmaEx(i,j)=ep*sigmam*((real(npmlx-i,8)+0.5d0)/real(npmlx,8))**order
		this%alphaEx(i,j)=this%sigmaEx(i,j)*dt/ep
		this%f1x(i,j)=dexp(-this%alphaEx(i,j))
		this%f2x(i,j)=(1.d0-this%f1x(i,j))/(z0*this%sigmaEx(i,j)*ds)
		if(i.lt.npmlx)then
		    this%sigmaHx(i,j)=mu0*sigmam*((real(npmlx-i,8))/real(npmlx,8))**order
		    this%alphaHx(i,j)=this%sigmaHx(i,j)*dt/mu0
		    this%g1x(i,j)=dexp(-this%alphaHx(i,j))
		    this%g2x(i,j)=z0*(1.d0-this%g1x(i,j))/(this%sigmaHx(i,j)*ds)
		endif
	    enddo
	enddo
	end subroutine init_varx1

	subroutine init_varx2(this,xi,xf,yi,yf,ep,sigmam,order,dt,ds)
	class(pmly) :: this
	integer, intent(in) :: xi,xf,yi,yf
	real*8, intent(in) :: ep,sigmam,order,dt,ds
	integer :: i,ii,j,npmlx
	real*8 :: z0
	npmlx=this%npmlx
	z0=dsqrt(mu0/ep0)
	    do j=yi,yf,1
		ii=npmlx
		do i=xi,xf,1
		    if(ii.lt.npmlx)then
			this%sigmaEx(i,j)=ep*sigmam*((real(npmlx-ii,8))/real(npmlx,8))**order
		        this%alphaEx(i,j)=this%sigmaEx(i,j)*dt/ep
			this%f1x(i,j)=dexp(-this%alphaEx(i,j))
			this%f2x(i,j)=(1.d0-this%f1x(i,j))/(z0*this%sigmaEx(i,j)*ds)
		    endif
		    this%sigmaHx(i,j)=mu0*sigmam*((real(npmlx-ii,8)+0.5d0)/real(npmlx,8))**order
		    this%alphaHx(i,j)=this%sigmaHx(i,j)*dt/mu0
		    this%g1x(i,j)=dexp(-this%alphaHx(i,j))
		    this%g2x(i,j)=z0*(1.d0-this%g1x(i,j))/(this%sigmaHx(i,j)*ds)
		    ii=ii-1
		enddo
	    enddo
	end subroutine init_varx2
    end module
