	program my3dfdtd
	use constants_mod
	use mydatatypes
	use simpleDielectric_mod
	use gold_mod
	use BERpmly_mod
	implicit none
	include 'mpif.h'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	control variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	logical, parameter :: PBC_x=.true.,PBC_y=.false.,PBC_z=.true.
	logical, parameter :: printTimesMap=.false.
	logical, parameter :: printTimesLine=.true.
	logical, parameter :: printTimesVTK=.false.

	character(len=4) :: framestr
	character(len=6) :: rankstr,rankstr2
	character(len=600) :: tmpFilename,names,msg
	integer, parameter :: pf=1
	integer :: nout,dn,n,frame,m,i,j,k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	MPI constants and variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer :: nproces,ncolumn,ierror,myrank,status(mpi_status_size)
	integer :: COMM_3D
	integer :: tagE1_WE,tagE2_WE,tagH1_WE,tagH2_WE
	integer :: tagE1_NS,tagE2_NS,tagH1_NS,tagH2_NS
	integer :: tagE1_UD,tagE2_UD,tagH1_UD,tagH2_UD
	integer :: tagE3_WE,tagE4_WE,tagH3_WE,tagH4_WE
	integer :: tagE3_NS,tagE4_NS,tagH3_NS,tagH4_NS
	integer :: tagE3_UD,tagE4_UD,tagH3_UD,tagH4_UD
	integer, parameter :: NORTH=1
	integer, parameter :: SOUTH=2
	integer, parameter :: EAST=3
	integer, parameter :: WEST=4
	integer, parameter :: UP=5
	integer, parameter :: DOWN=6
	integer:: neigborN,neigborS,neigborE,neigborW,neigborU,neigborD
	integer :: TypeMatrXZ,TypeMatrYZ,countXY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Grid process Layout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real*8, parameter :: lambda0=633.d0*nano
	real*8, parameter :: lambda1=lambda0/3.88d0
	real*8, parameter :: f0=c0/lambda0
	real*8, parameter :: stabilityFactor=0.99d0 !0.866~Ku=0.5
	real*8, parameter :: Ku=0.5d0!stabilityFactor/dsqrt(3.d0)
	real*8, parameter :: ds=3.999d0*nano 
	real*8, parameter :: dt=ku*ds/c0
	integer, parameter :: lambda0ds=floor(lambda0/ds)
	integer, parameter :: lambda1ds=floor(lambda1/ds)
	real*8, parameter :: xreal=500.d0*nano
	real*8, parameter :: yreal=5064.d0*nano
	real*8, parameter :: zreal=25.d0*ds
	integer, parameter :: Xmax=1,Ymax=6,Zmax=1	!number of processor in each direction
	integer, parameter :: nx=floor(xreal/ds)
	integer, parameter :: ny=floor(yreal/ds)
	integer, parameter :: nz=floor(zreal/ds)
	integer, parameter :: ic=nx/2
	integer, parameter :: jc=ny/2
	integer, parameter :: kc=nz/2
	integer, parameter :: n_x=nx/Xmax
	integer, parameter :: n_y=ny/Ymax
	integer, parameter :: n_z=nz/Zmax
	integer, parameter :: i_c=n_x/2
	integer, parameter :: j_c=n_y/2
	integer, parameter :: k_c=n_z/2
	integer, parameter :: Ndim(3)=(/ Xmax,Ymax,Zmax /)
	integer, dimension(3) :: coords
	integer :: xcord,ycord,zcord
	logical, parameter :: periodic(3)=(/ .false.,.false.,.false. /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	field variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real*8 :: t
	real*8, dimension(0:N_x+1,0:N_y+1,0:N_z+1), target :: Hx,Hy,Hz,Ex,Ey,Ez,Dx,Dy,Dz
	real*8, dimension(:,:,:), pointer :: PEx,PEy,PEz,PDx,PDy,PDz,PHx,PHy,PHz
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	Material variables
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	type(ConstantDielectric) :: Si,SiO2,air
	type(Gold) :: Au
	type(Irange) :: absIxrange,absIyrange,absIzrange
	type(Irange) :: locIxrange,locIyrange,locIzrange
	integer, parameter,dimension(3) :: rank_dims=(/ n_x,n_y,n_z /)
	real*8, parameter :: n_SiO2=1.46d0
	real*8, parameter :: n_Si=3.88d0
	real*8, parameter :: n_air=1.d0
	real*8, parameter :: SiO2_ep=n_SiO2**2
	real*8, parameter :: Si_ep=n_Si**2
	real*8, parameter :: air_ep=n_air**2
	real*8 :: goldLength, goldWidth, goldThick
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		PML variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	type(pmly) :: mypmly
	integer, parameter :: npmlx=10,npmly=10,npmlz=10
	integer, parameter, dimension(3) :: npmls=(/ npmlx,npmly,npmlz /)
	real*8, parameter :: order=3.d0
	real*8, parameter :: R0=dexp(-16.d0)
	real*8, parameter :: sigmam=(order+1.d0)*c0*log(1.d0/R0)/(2.d0*npmly*ds)
!	Ndim,coords,rank_dims
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	source variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real*8, parameter :: periodoin=1.d0/f0
	real*8, parameter :: w0=2.d0*pi*f0
	real*8, parameter :: beta=1.d0
	real*8, parameter :: st=beta*periodoin
	real*8, parameter :: t0=4.d0*st
	integer, parameter :: sposx=ic
	integer, parameter :: sposy=npmly+2
	integer, parameter :: sposz=kc
	real*8 :: gauss_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		fourier Transform variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, parameter :: freqsteps=100

    real*8, parameter :: lambdainit=500d-9
    real*8, parameter :: lambdaend=800d-9
    real*8, parameter :: finit=c0/lambdaend
    real*8, parameter :: fend=c0/lambdainit
    real*8, parameter :: df=(fend-finit)/freqsteps

    real*8,dimension(n_x,n_z,freqsteps) :: real_hzin,img_hzin
    real*8,dimension(n_x,n_z,freqsteps) :: real_hzt,img_hzt
    real*8,dimension(n_x,n_z,freqsteps) :: real_hzr,img_hzr

    real*8,dimension(n_x,n_z,freqsteps) :: real_hxin,img_hxin
    real*8,dimension(n_x,n_z,freqsteps) :: real_hxt,img_hxt
    real*8,dimension(n_x,n_z,freqsteps) :: real_hxr,img_hxr

    real*8,dimension(n_x,n_z,freqsteps) :: real_exin,img_exin
    real*8,dimension(n_x,n_z,freqsteps) :: real_ext,img_ext
    real*8,dimension(n_x,n_z,freqsteps) :: real_exr,img_exr

    real*8,dimension(n_x,n_z,freqsteps) :: real_ezin,img_ezin
    real*8,dimension(n_x,n_z,freqsteps) :: real_ezt,img_ezt
    real*8,dimension(n_x,n_z,freqsteps) :: real_ezr,img_ezr

    real*8,dimension(freqsteps) :: arg
    real*8 :: Syin,Syr,Syt,freq,lambda
!	incidence monitor
    integer, parameter :: ptinx=ic
    integer, parameter :: ptiny=sposy+2
    integer, parameter :: ptinz=kc
!	reflection monitor
    integer, parameter :: ptrx=ptinx
    integer, parameter :: ptry=ptiny
    integer, parameter :: ptrz=ptinz
!	transmission monitor
    integer, parameter :: pttz=kc
    integer, parameter :: pttx=ic
    integer :: ptty
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	time duration variables
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8 :: time4incident,distance4incident,vel4incident
	integer :: tfinStop
	real*8 :: timeSim
	integer :: Nt
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	Fourier Transform initialization
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    real_exin=0.d0
    img_exin=0.d0
    real_exr=0.d0
    img_exr=0.d0
    real_ext=0.d0
    img_ext=0.d0

    real_ezin=0.d0
    img_ezin=0.d0
    real_ezr=0.d0
    img_ezr=0.d0
    real_ezt=0.d0
    img_ezt=0.d0

    real_hzin=0.d0
    img_hzin=0.d0
    real_hzr=0.d0
    img_hzr=0.d0
    real_hzt=0.d0
    img_hzt=0.d0

    real_hxin=0.d0
    img_hxin=0.d0
    real_hxr=0.d0
    img_hxr=0.d0
    real_hxt=0.d0
    img_hxt=0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Initialize MPI environment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call my_mpi_init()
	write(rankstr,'(I4)') myrank

!	Matrix blocks creations for MPI messaging
	call create_mpi_planes()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	PML initiation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call mypmly%initialize(air_ep*ep0,Si_ep*ep0,sigmam,npmls,order,dt,ds,coords,rank_dims,Ndim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Initialize variables (set everything to zero)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Ex=0.d0
	Ey=0.d0
	Ez=0.d0
	Hx=0.d0
	Hy=0.d0
	Hz=0.d0

	PEx=>Ex
	PEy=>Ey
	PEz=>Ez

	PDx=>Dx
	PDy=>Dy
	PDz=>Dz

	PHx=>Hx
	PHy=>Hy
	PHz=>Hz
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	air initialization
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	air%xinit=1
	air%xend=Nx
	air%yinit=npmly+1
	air%yend=0.8*Ny!floor(3792.d0*nano/ds)
	air%zinit=1
	air%zend=Nz
	air%constant_ep=air_ep
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	gold initialization
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	goldLength=100.d0*nano
	goldWidth=50.d0*nano
	goldThick=20.d0*nano

	Au%xinit=ic-0.5*floor(goldLength/ds)+1
	Au%xend=ic+0.5*floor(goldLength/ds)
!	Au%zinit=kc-0.5*floor(goldWidth/ds)+1
!	Au%zend=kc+0.5*floor(goldWidth/ds)
!	Au%xinit=1
!	Au%xend=Nx
	Au%zinit=1
	Au%zend=Nz
	Au%yinit=air%yend-floor(goldThick/ds)+1
	Au%yend=air%yend
	Au%dt=dt
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	SiO2 initialization
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	SiO2%xinit=1
	SiO2%xend=Nx
	SiO2%yinit=air%yend+1
	SiO2%yend=SiO2%yinit+floor(284.d0*nano/ds)-1
	SiO2%zinit=1
	SiO2%zend=Nz
	SiO2%constant_ep=SiO2_ep
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	Si initialization
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	Si%xinit=1
	Si%xend=Nx
	Si%yinit=SiO2%yend+1
	Si%yend=Ny-npmly
	Si%zinit=1
	Si%zend=Nz
	Si%constant_ep=Si_ep
	ptty=0.5*(Si%yend+Si%yinit)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	time windowing
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    distance4incident=ds*real(air%yend-sposy,8)+beta*4.5d0*lambda0/n_air
    vel4incident=c0/n_air
    time4incident=distance4incident/vel4incident
    tfinStop=floor(time4incident/dt)
    Nt=4*tfinStop
    dn=aint(Nt/20.d0)	! printing rate per frame
    nout=dn
    frame=0

	call mpi_barrier(mpi_comm_world,ierror)
	call air%initialize_material(absIxrange,absIyrange,absIzrange,rank_dims)
	call SiO2%initialize_material(absIxrange,absIyrange,absIzrange,rank_dims)
	call Si%initialize_material(absIxrange,absIyrange,absIzrange,rank_dims)
	call Au%initialize_material(absIxrange,absIyrange,absIzrange,rank_dims)

    if(inRange(absIyrange,Au%yinit).or.inRange(absIyrange,Au%yend))then
        tmpFilename="r"//trim(adjustl(rankstr))//"goldVars.dat"
        open(unit=30,file=tmpFilename,status="unknown")
        do j=1,n_y,1
	    do i=1,n_x,1
	    write(30,'(2F15.5,12E20.10)')&
	(i+absIxrange%init-1)*ds/nano,(j+absIyrange%init-1)*ds/nano,&
	Au%ha(i,j,k_c),Au%hb(i,j,k_c),Au%ga(i,j,k_c),&
	Au%gb1(i,j,k_c),Au%gb2(i,j,k_c),&
	Au%gc1(i,j,k_c),Au%gc2(i,j,k_c)
	enddo
	write(30,*)
        enddo
    endif

	if(myrank.eq.0)then
	    print *,"sigmam=",sigmam
	    print *,"npmls=",npmls
	    print *,"Nx=",Nx,"Ny=",Ny,"Nz=",Nz
	    print *,"N_x=",N_x,"N_y=",N_y,"N_z=",N_z
	    print *,"ds=",ds
	    print *,"dt=",dt
	    print *,"dn=",dn
	    print *,"Nt=",Nt
	    print *,"Nt/dn=",Nt/dn
	    print *,"lambda0=",lambda0,floor(lambda0/ds)
	    print *,"lambda1=",lambda1,floor(lambda1/ds)
	    print *,"source_sposy=",sposy
	    print *,"Au:xi",Au%xinit,"Au:xf",Au%xend
	    print *,"Au:yi",Au%yinit,"Au:yf",Au%yend
	    print *,"Au:zi",Au%zinit,"Au:zf",Au%zend
	    print *,"Si:xi",Si%xinit,"Si:xf",Si%xend
	    print *,"Si:yi",Si%yinit,"Si:yf",Si%yend
	    print *,"Si:zi",Si%zinit,"Si:zf",Si%zend
	    print *,"Si:ep",Si%constant_ep,dsqrt(Si%constant_ep)
	    print *,"SiO2:xi",SiO2%xinit,"SiO2:xf",SiO2%xend
	    print *,"SiO2:yi",SiO2%yinit,"SiO2:yf",SiO2%yend
	    print *,"SiO2:zi",SiO2%zinit,"SiO2:zf",SiO2%zend
	    print *,"SiO2:ep",SiO2%constant_ep,dsqrt(SiO2%constant_ep)
	    print *,"air:xi",air%xinit,"air:xf",air%xend
	    print *,"air:yi",air%yinit,"air:yf",air%yend
	    print *,"air:zi",air%zinit,"air:zf",air%zend
	    print *,"air:ep",air%constant_ep,dsqrt(air%constant_ep)
	    print *,"time4incident=",time4incident
	    print *,"distance4incident=",distance4incident
	    print *,"vel4incident=",vel4incident,vel4incident/c0
	    print *,"tfinStop=",tfinStop
	    print *,"xBDlow=",mypmly%xBDlow,"xBDhigh=",mypmly%xBDhigh
	    print *,"yBDlow=",mypmly%yBDlow,"yBDhigh=",mypmly%yBDhigh
	    print *,"zBDlow=",mypmly%zBDlow,"zBDhigh=",mypmly%zBDhigh
	    print *,"ic=",ic,"jc=",jc,"kc=",kc
	    print *,"i_c=",i_c,"j_c=",j_c,"k_c=",k_c
	    print *,"ptiny=",ptiny,"ptry=",ptry,"ptty=",ptty
	endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c		open file for input fields of the FT
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	print *,myrank,"open file for input fields of the FT"
    freq=finit
    do m=1,freqsteps,1
        arg(m)=2.d0*pi*freq
        freq=freq+df
    enddo

    if(inRange(absIxrange,ptinx).and.&
       inRange(absIyrange,ptiny).and.&
       inRange(absIzrange,ptinz))then
    tmpFilename="transforms/r"//trim(adjustl(rankstr))//"infld.dat"
        open(unit=2,file=tmpFilename,status="unknown")
    endif

    if(inRange(absIxrange,ptrx).and.&
       inRange(absIyrange,ptry).and.&
       inRange(absIzrange,ptrz))then
        tmpFilename="transforms/r"//trim(adjustl(rankstr))//"rfld.dat"
        open(unit=3,file=tmpFilename,status="unknown")
    endif

    if(inRange(absIxrange,pttx).and.&
       inRange(absIyrange,ptty).and.&
       inRange(absIzrange,pttz))then
        tmpFilename="transforms/r"//trim(adjustl(rankstr))//"tfld.dat"
        open(unit=4,file=tmpFilename,status="unknown")
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Start FDTD time loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do n=1,Nt,1
!	do n=1,1,1
	t=real(n-1,8)*dt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	calculate D field in simulation space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call updateD(locIxrange%init,locIxrange%ends,&
		      locIyrange%init,locIyrange%ends,&
		      locIzrange%init,locIzrange%ends)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	source definition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	plane wave
	if(inRange(absIyrange,sposy))then
	    if(n.eq.1)then
		print *,"source:",myrank
	    endif
	    gauss_t=dexp(-0.5d0*((t-t0)/st)**2)*dsin(w0*t)
	    do k=1,n_z,1
		do i=1,n_x,1
    Dx(i,sposy-absIyrange%init+1,k)=Dx(i,sposy-absIyrange%init+1,k)+gauss_t*n_air
		enddo
	    enddo
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Calculate E from D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	call updateE_vaccum(1,n_x,1,n_y,1,n_z)
        call air%update_E(pex,pey,pez,pdx,pdy,pdz)
        call Au%update_E(pex,pey,pez,pdx,pdy,pdz)
        call SiO2%update_E(pex,pey,pez,pdx,pdy,pdz)
        call Si%update_E(pex,pey,pez,pdx,pdy,pdz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Calculate E field in PML layers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call mypmly%updateE(PEx,PEy,PEz,PHx,PHy,PHz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	communication of the E-fields components at the subspace boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    call interchange_E()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Calculate H field inside simulation space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call updateH(locIxrange%init,locIxrange%ends,&
		      locIyrange%init,locIyrange%ends,&
		      locIzrange%init,locIzrange%ends)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Calculate H field in PML layers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call mypmly%updateH(PEx,PEy,PEz,PHx,PHy,PHz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	communication of the H-fields components at the subspace boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call interchange_H()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Calculate TFs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(inRange(absIyrange,ptiny).and.&
    n.lt.tfinStop)then
	if(inRange(absIxrange,ptinx).and.&
	    inRange(absIzrange,ptinz))then
	    write(2,'(7E30.19)')n*dt,&
    ex(ptinx-absIxrange%init+1,ptiny-absIyrange%init+1,ptinz-absIzrange%init+1),&
    ey(ptinx-absIxrange%init+1,ptiny-absIyrange%init+1,ptinz-absIzrange%init+1),&
    ez(ptinx-absIxrange%init+1,ptiny-absIyrange%init+1,ptinz-absIzrange%init+1),&
    hx(ptinx-absIxrange%init+1,ptiny-absIyrange%init+1,ptinz-absIzrange%init+1),&
    hy(ptinx-absIxrange%init+1,ptiny-absIyrange%init+1,ptinz-absIzrange%init+1),&
    hz(ptinx-absIxrange%init+1,ptiny-absIyrange%init+1,ptinz-absIzrange%init+1)
	endif
        if(n.eq.1)then
	print *,"TFin: myrank=",myrank,ptiny-absIyrange%init+1
        endif
    do m=1,freqsteps,1
        do k=1,n_z,1
	    do i=1,n_x,1
real_exin(i,k,m)=real_exin(i,k,m)+dcos(arg(m)*t)*ex(i,ptiny-absIyrange%init+1,k)*dt
 img_exin(i,k,m)= img_exin(i,k,m)+dsin(arg(m)*t)*ex(i,ptiny-absIyrange%init+1,k)*dt

real_ezin(i,k,m)=real_ezin(i,k,m)+dcos(arg(m)*t)*ez(i,ptiny-absIyrange%init+1,k)*dt
 img_ezin(i,k,m)= img_ezin(i,k,m)+dsin(arg(m)*t)*ez(i,ptiny-absIyrange%init+1,k)*dt

real_hxin(i,k,m)=real_hxin(i,k,m)+dcos(arg(m)*t)*hx(i,ptiny-absIyrange%init+1,k)*dt
 img_hxin(i,k,m)= img_hxin(i,k,m)+dsin(arg(m)*t)*hx(i,ptiny-absIyrange%init+1,k)*dt

real_hzin(i,k,m)=real_hzin(i,k,m)+dcos(arg(m)*t)*hz(i,ptiny-absIyrange%init+1,k)*dt
 img_hzin(i,k,m)= img_hzin(i,k,m)+dsin(arg(m)*t)*hz(i,ptiny-absIyrange%init+1,k)*dt
	    enddo
        enddo
    enddo
    endif

    if(inRange(absIyrange,ptry).and.&
    n.gt.tfinStop)then
        if(inRange(absIxrange,ptrx).and.&
           inRange(absIzrange,ptrz))then
	write(3,'(7E30.19)')n*dt,&
    ex(ptrx-absIxrange%init+1,ptry-absIyrange%init+1,ptrz-absIzrange%init+1),&
    ey(ptrx-absIxrange%init+1,ptry-absIyrange%init+1,ptrz-absIzrange%init+1),&
    ez(ptrx-absIxrange%init+1,ptry-absIyrange%init+1,ptrz-absIzrange%init+1),&
    hx(ptrx-absIxrange%init+1,ptry-absIyrange%init+1,ptrz-absIzrange%init+1),&
    hy(ptrx-absIxrange%init+1,ptry-absIyrange%init+1,ptrz-absIzrange%init+1),&
    hz(ptrx-absIxrange%init+1,ptry-absIyrange%init+1,ptrz-absIzrange%init+1)
        endif
        if(n.eq.1)then
        print *,"TFr: myrank=",myrank,ptry-absIyrange%init+1
        endif
    do m=1,freqsteps,1
	do k=1,n_z,1
	    do i=1,n_x,1
real_exr(i,k,m)=real_exr(i,k,m)+dcos(arg(m)*t)*ex(i,ptry-absIyrange%init+1,k)*dt
 img_exr(i,k,m)= img_exr(i,k,m)+dsin(arg(m)*t)*ex(i,ptry-absIyrange%init+1,k)*dt

real_ezr(i,k,m)=real_ezr(i,k,m)+dcos(arg(m)*t)*ez(i,ptry-absIyrange%init+1,k)*dt
 img_ezr(i,k,m)= img_ezr(i,k,m)+dsin(arg(m)*t)*ez(i,ptry-absIyrange%init+1,k)*dt

real_hxr(i,k,m)=real_hxr(i,k,m)+dcos(arg(m)*t)*hx(i,ptry-absIyrange%init+1,k)*dt
 img_hxr(i,k,m)= img_hxr(i,k,m)+dsin(arg(m)*t)*hx(i,ptry-absIyrange%init+1,k)*dt

real_hzr(i,k,m)=real_hzr(i,k,m)+dcos(arg(m)*t)*hz(i,ptry-absIyrange%init+1,k)*dt
 img_hzr(i,k,m)= img_hzr(i,k,m)+dsin(arg(m)*t)*hz(i,ptry-absIyrange%init+1,k)*dt
	    enddo
	enddo
    enddo
    endif

    if(inRange(absIyrange,ptty))then
        if(inRange(absIxrange,pttx).and.&
           inRange(absIzrange,pttz))then
	write(4,'(7E30.19)')n*dt,&
    ex(pttx-absIxrange%init+1,ptty-absIyrange%init+1,pttz-absIzrange%init+1),&
    ey(pttx-absIxrange%init+1,ptty-absIyrange%init+1,pttz-absIzrange%init+1),&
    ez(pttx-absIxrange%init+1,ptty-absIyrange%init+1,pttz-absIzrange%init+1),&
    hx(pttx-absIxrange%init+1,ptty-absIyrange%init+1,pttz-absIzrange%init+1),&
    hy(pttx-absIxrange%init+1,ptty-absIyrange%init+1,pttz-absIzrange%init+1),&
    hz(pttx-absIxrange%init+1,ptty-absIyrange%init+1,pttz-absIzrange%init+1)
        endif
        if(n.eq.1)then
        print *,"TFt: myrank=",myrank,ptty-absIyrange%init+1
        endif
    do m=1,freqsteps,1
	do k=1,n_z,1
	    do i=1,n_x,1
real_ext(i,k,m)=real_ext(i,k,m)+dcos(arg(m)*t)*ex(i,ptty-absIyrange%init+1,k)*dt
 img_ext(i,k,m)= img_ext(i,k,m)+dsin(arg(m)*t)*ex(i,ptty-absIyrange%init+1,k)*dt

real_ezt(i,k,m)=real_ezt(i,k,m)+dcos(arg(m)*t)*ez(i,ptty-absIyrange%init+1,k)*dt
 img_ezt(i,k,m)= img_ezt(i,k,m)+dsin(arg(m)*t)*ez(i,ptty-absIyrange%init+1,k)*dt

real_hxt(i,k,m)=real_hxt(i,k,m)+dcos(arg(m)*t)*hx(i,ptty-absIyrange%init+1,k)*dt
 img_hxt(i,k,m)= img_hxt(i,k,m)+dsin(arg(m)*t)*hx(i,ptty-absIyrange%init+1,k)*dt

real_hzt(i,k,m)=real_hzt(i,k,m)+dcos(arg(m)*t)*hz(i,ptty-absIyrange%init+1,k)*dt
 img_hzt(i,k,m)= img_hzt(i,k,m)+dsin(arg(m)*t)*hz(i,ptty-absIyrange%init+1,k)*dt
	    enddo
	enddo
    enddo
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	imprime muestra temporal
	if(nout.eq.n)then
	    nout=nout+dn
	    frame=frame+1
	    write(framestr,'(I4)') frame

	    if(printTimesVTK.and.((ycord.eq.0.or.ycord.eq.(Ymax-1))))then
		tmpFilename="timesVTK/r"//trim(adjustl(rankstr))&
		    //"Ex_time"//trim(adjustl(framestr))//".vtk"
		open(unit=100+frame,file=tmpFilename,status="unknown")
		call makeVTKheader(100+frame,'Ex',N_x,N_y,N_z)
		do k=1,N_z,1
		    do j=1,N_y,1
			do i=1,N_x,1
			    write(100+frame,"(E30.19)",advance="no") Hz(i,j,k)
			enddo
		    enddo
		    write(100+frame,*)
		enddo
	    endif

	    if(inRange(absIxrange,ic).and.&
		inRange(absIzrange,kc).and.&
		printTimesLine)then
	tmpFilename="timesLiney/r"//trim(adjustl(rankstr))&
	    //"E_time"//trim(adjustl(framestr))//".dat"
	    open(unit=400+frame,file=tmpFilename,status="unknown")
		do j=1,n_y,1
!		    write(400+frame,*)&
!	    (j+absIyrange%init-1)*ds/nano,&
!	    Ex(ic-absIxrange%init+1,j,kc-absIzrange%init+1),&
!	    Hz(ic-absIxrange%init+1,j,kc-absIzrange%init+1)
		    write(400+frame,'(7F10.5)')&
	    (j+absIyrange%init-1)*ds/nano,&
	    Ex(ic-absIxrange%init+1,j,kc-absIzrange%init+1),&
	    Ey(ic-absIxrange%init+1,j,kc-absIzrange%init+1),&
	    Ez(ic-absIxrange%init+1,j,kc-absIzrange%init+1),&
	    Hx(ic-absIxrange%init+1,j,kc-absIzrange%init+1),&
	    Hy(ic-absIxrange%init+1,j,kc-absIzrange%init+1),&
	    Hz(ic-absIxrange%init+1,j,kc-absIzrange%init+1)
		enddo
	    endif

	    if(inRange(absIzrange,kc).and.printTimesMap)then
	tmpFilename="timesGnuplotmap/r"//trim(adjustl(rankstr))&
	    //"Ekc_time"//trim(adjustl(framestr))//".dat"
	    open(unit=500+frame,file=tmpFilename,status="unknown")

		do j=1,N_y,pf
		    do i=1,N_x,pf
	write(500+frame,'(8E30.19)')&
	    (i+absIxrange%init-1)*ds/nano,&
	    (j+absIyrange%init-1)*ds/nano,&
	    Ex(i,j,kc-absIzrange%init+1),&
	    Ey(i,j,kc-absIzrange%init+1),&
	    Ez(i,j,kc-absIzrange%init+1),&
	    Hx(i,j,kc-absIzrange%init+1),&
	    hy(i,j,kc-absIzrange%init+1),&
	    Hz(i,j,kc-absIzrange%init+1)
		    enddo
		    write(500+frame,*)
		enddo
	    endif
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	enddo
!		End FDTD time loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		print TFs to file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(inRange(absIyrange,ptiny))then
        tmpFilename="transforms/r"//trim(adjustl(rankstr))//"Syin.dat"
        open(unit=10,file=tmpFilename,status="unknown")
        do m=1,freqsteps,1
	Syin=0.d0
	do k=1,n_z,1
	    do i=1,n_x,1
	    Syin=Syin+(&
	    (real_hxin(i,k,m)*real_ezin(i,k,m))+&
	    ( img_hxin(i,k,m)* img_ezin(i,k,m))-&
	    (real_hzin(i,k,m)*real_exin(i,k,m))-&
	    ( img_hzin(i,k,m)* img_exin(i,k,m))&
	    )*ds
	    enddo
	enddo
	lambda=2.d0*pi*c0/arg(m)
	write(10,*)lambda/nano,Syin
        enddo
    endif

    if(inRange(absIyrange,ptry))then
        tmpFilename="transforms/r"//trim(adjustl(rankstr))//"Syr.dat"
        open(unit=11,file=tmpFilename,status="unknown")

        tmpFilename="transforms/r"//trim(adjustl(rankstr))//"tf_r.dat"
        open(unit=21,file=tmpFilename,status="unknown")
        do m=1,freqsteps,1
	    Syr=0.d0
	    lambda=2.d0*pi*c0/arg(m)
	    do k=1,n_z,1
		do i=1,n_x,1
		    Syr=Syr+(&
		    (real_hxr(i,k,m)*real_ezr(i,k,m))+&
		    ( img_hxr(i,k,m)* img_ezr(i,k,m))-&
		    (real_hzr(i,k,m)*real_exr(i,k,m))-&
		    ( img_hzr(i,k,m)* img_exr(i,k,m))&
		    )*ds
		    if(k.eq.kc)then
			write(21,'(2F15.5,4E20.10)')lambda/nano,i*ds/nano,&
			real_hxr(i,k,m)*real_ezr(i,k,m),&
			img_hxr(i,k,m)*img_ezr(i,k,m),&
			real_hzr(i,k,m)*real_exr(i,k,m),&
			img_hzr(i,k,m)* img_exr(i,k,m)
		    endif
		enddo
		if(k.eq.kc)then
		    write(21,*)
		endif
	    enddo
	    write(11,*)lambda/nano,Syr

	enddo
    endif

    if(inRange(absIyrange,ptty))then
        tmpFilename="transforms/r"//trim(adjustl(rankstr))//"Syt.dat"
        open(unit=12,file=tmpFilename,status="unknown")
        do m=1,freqsteps,1
        Syt=0.d0
        do k=1,n_z,1
	do i=1,n_x,1
	    Syt=Syt+(&
	    (real_hxt(i,k,m)*real_ezt(i,k,m))+&
	    ( img_hxt(i,k,m)* img_ezt(i,k,m))-&
	    (real_hzt(i,k,m)*real_ext(i,k,m))-&
	    ( img_hzt(i,k,m)* img_ext(i,k,m))&
	    )*ds
	enddo
        enddo
    lambda=2.d0*pi*c0/arg(m)
    write(12,*)lambda/nano,Syt
        enddo
    endif

	call mpi_finalize(ierror)

	contains
	include 'include/myfdtd_func.f90'
	include 'include/mympi_func.f90'

	function inRange(haystack,needle) result(found)
	    type(Irange),intent(in) :: haystack
	    integer, intent(in) :: needle
	    logical :: found
	    found =.false.
	    if(haystack%init.le.needle.and.&
		haystack%ends.ge.needle)then
		found =.true.
	    endif
	end function inRange

!    FUNCTION Replace_Text (s,text,rep)  RESULT(outs)
!        CHARACTER(*)        :: s,text,rep
!        CHARACTER(LEN(s)+100) :: outs     ! provide outs with extra 100 char len
!        INTEGER             :: i, nt, nr
!
!        outs = s ; nt = LEN(text) ; nr = LEN(rep)
!        DO
!	i = INDEX(outs,text(:nt)) ; IF (i == 0) EXIT
!	outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
!        END DO
!    END FUNCTION Replace_Text

	end program my3dfdtd
