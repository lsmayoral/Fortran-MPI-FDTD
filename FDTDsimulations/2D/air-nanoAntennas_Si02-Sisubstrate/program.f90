	program my3dfdtd
	use constants_mod
	use mydatatypes
	use simpleDielectric_mod
	use gold_mod
	use BERpmly_HzExEy_mod
	implicit none
	include 'mpif.h'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	control variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	logical, parameter :: PBC_x=.true.,PBC_y=.false.
	logical, parameter :: printTimesMap=.true.
	logical, parameter :: printTimesLine=.true.
	logical, parameter :: printTimesVTK=.false.

	character(len=4) :: framestr
	character(len=6) :: rankstr,rankstr2
	character(len=600) :: tmpFilename,names
	integer, parameter :: pf=1
	integer :: nout,dn,n,frame,m,i,j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	MPI constants and variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer :: nproces,ncolumn,ierror,myrank,status(mpi_status_size)
	integer :: COMM_2D
	integer :: tagE1_WE,tagE2_WE,tagH1_WE,tagH2_WE
        integer :: tagE1_NS,tagE2_NS,tagH1_NS,tagH2_NS
        integer :: tagE3_WE,tagE4_WE,tagH3_WE,tagH4_WE
        integer :: tagE3_NS,tagE4_NS,tagH3_NS,tagH4_NS
        integer, parameter :: NORTH=1
        integer, parameter :: SOUTH=2
        integer, parameter :: EAST=3
	integer, parameter :: WEST=4
	integer:: neighborN,neighborS,neighborE,neighborW
	integer :: liney 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Grid process Layout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real*8, parameter :: lambda0=633.d0*nano
	real*8, parameter :: lambda1=lambda0/3.88d0
	real*8, parameter :: f0=c0/lambda0
	real*8, parameter :: stabilityFactor=0.99d0
	real*8, parameter :: Ku=0.5d0!stabilityFactor/dsqrt(2.d0)
	real*8, parameter :: ds=3.999d0*nano 
	real*8, parameter :: dt=ku*ds/c0
	integer, parameter :: lambda0ds=floor(lambda0/ds)
	integer, parameter :: lambda1ds=floor(lambda1/ds)
	real*8, parameter :: xreal=500.d0*nano
	real*8, parameter :: yreal=5064.d0*nano
	integer, parameter :: Xmax=1,Ymax=3	!number of processor in each direction
	integer, parameter :: nx=floor(xreal/ds)
	integer, parameter :: ny=floor(yreal/ds)
	integer, parameter :: ic=nx/2
	integer, parameter :: jc=ny/2
	integer, parameter :: n_x=nx/Xmax
	integer, parameter :: n_y=ny/Ymax
	integer, parameter :: i_c=n_x/2
	integer, parameter :: j_c=n_y/2
	integer, parameter :: Ndim(2)=(/ Xmax,Ymax /)
	integer, dimension(2) :: coords
	integer :: xcord,ycord
	logical, parameter :: periodic(2)=(/ .false.,.false. /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	field variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real*8 :: t
	real*8, dimension(0:N_x+1,0:N_y+1), target :: Hz,Ex,Ey,Dx,Dy
	real*8, dimension(:,:), pointer :: PEx,PEy,PDx,PDy,PHz
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	Material variables
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	type(ConstantDielectric) :: air,SiO2,Si
	type(Gold) :: Au
	type(Irange) :: absIxrange,absIyrange
	type(Irange) :: locIxrange,locIyrange
	integer, parameter,dimension(2) :: rank_dims=(/ n_x,n_y /)
	real*8, parameter :: n_SiO2=1.46d0
	real*8, parameter :: n_Si=3.88d0
	real*8, parameter :: n_air=1.d0
	real*8, parameter :: SiO2_ep=n_SiO2**2
	real*8, parameter :: Si_ep=n_Si**2
	real*8, parameter :: air_ep=n_air**2
	real*8 :: goldLength, goldThick
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		PML variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	type(pmly) :: mypmly
	integer, parameter :: npmlx=10,npmly=10
	integer, parameter, dimension(2) :: npmls=(/ npmlx,npmly/)
	real*8, parameter :: order=3.d0
	real*8, parameter :: R0=dexp(-16.d0) !16->(0.0025,0.01) 18->(0.003,0.01)
	real*8, parameter :: sigmam=(order+1.d0)*c0*log(1.d0/R0)/(2.d0*npmly*ds)
!	npml=5, R0=exp(-8) => r=0.0025
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
	real*8 :: gauss_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		fourier Transform variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, parameter :: freqsteps=200

    real*8, parameter :: lambdainit=500d-9
    real*8, parameter :: lambdaend=800d-9
    real*8, parameter :: finit=c0/lambdaend
    real*8, parameter :: fend=c0/lambdainit
    real*8, parameter :: df=(fend-finit)/freqsteps

    real*8,dimension(n_x,freqsteps) :: real_hzin,img_hzin
    real*8,dimension(n_x,freqsteps) :: real_hzt,img_hzt
    real*8,dimension(n_x,freqsteps) :: real_hzr,img_hzr

    real*8,dimension(n_x,freqsteps) :: real_exin,img_exin
    real*8,dimension(n_x,freqsteps) :: real_ext,img_ext
    real*8,dimension(n_x,freqsteps) :: real_exr,img_exr

    real*8,dimension(freqsteps) :: arg
    real*8 :: Syin,Syr,Syt,freq,lambda
!	incidence monitor
    integer, parameter :: ptinx=ic
    integer, parameter :: ptiny=sposy+2
!	reflection monitor
    integer, parameter :: ptrx=ptinx
    integer, parameter :: ptry=ptiny
!	transmission monitor
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

    real_hzin=0.d0
    img_hzin=0.d0
    real_hzr=0.d0
    img_hzr=0.d0
    real_hzt=0.d0
    img_hzt=0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Initialize MPI environment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call my_mpi_init()
	write(rankstr,'(I4)') myrank

!	Matrix blocks creations for MPI messaging
	call create_mpi_planes()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Initialize variables (set everything to zero)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Ex=0.d0
	Ey=0.d0
	Hz=0.d0

	PEx=>Ex
	PEy=>Ey

	PDx=>Dx
	PDy=>Dy

	PHz=>Hz
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	air initialization
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	air%xinit=1
	air%xend=Nx
	air%yinit=npmly
	air%yend=0.8*Ny!floor(2532.d0*nano/ds)
	air%constant_ep=air_ep
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	gold initialization
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	goldLength=100.d0*nano
	goldThick=20.d0*nano

	Au%xinit=ic-0.5*floor(goldLength/ds)+1
	Au%xend=ic+0.5*floor(goldLength/ds)
!	Au%xinit=1
!	Au%xend=Nx
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
	SiO2%constant_ep=SiO2_ep
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	Si initialization
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	Si%xinit=1
	Si%xend=Nx
	Si%yinit=SiO2%yend+1
	Si%yend=Ny-npmly
	Si%constant_ep=Si_ep
	ptty=0.5*(Si%yinit+Si%yend)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	time windowing
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    distance4incident=ds*real(air%yend-sposy,8)+beta*4.5d0*lambda0/n_air
    vel4incident=c0/n_air
    time4incident=distance4incident/vel4incident
    tfinStop=floor(time4incident/dt)
!    Nt=0.75*tfinStop
    Nt=4*tfinStop
    dn=aint(Nt/20.d0)	! printing rate per frame
    nout=dn
!    nout=4*aint(Nt/20.d0)	! printing rate per frame
!    dn=1
    frame=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	PML initiation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call mypmly%initialize(air_ep*ep0,Si_ep*ep0,sigmam,npmls,order,dt,ds,coords,rank_dims,Ndim)

	 call air%initialize_material(absIxrange,absIyrange,rank_dims)
	  call Si%initialize_material(absIxrange,absIyrange,rank_dims)
	call SiO2%initialize_material(absIxrange,absIyrange,rank_dims)
	  call Au%initialize_material(absIxrange,absIyrange,rank_dims)

	if(inRange(absIyrange,Au%yinit).or.inRange(absIyrange,Au%yend))then
	    tmpFilename="r"//trim(adjustl(rankstr))//"goldVars.dat"
	    open(unit=30,file=tmpFilename,status="unknown")
	    do j=1,n_y,1
		do i=1,n_x,1
		    write(30,'(2F15.5,12E20.10)')&
		(j+absIyrange%init-1)*ds/nano,(i+absIxrange%init-1)*ds/nano,&
		Au%ha(i,j),Au%hb(i,j),Au%ga(i,j),&
		Au%gb1(i,j),Au%gb2(i,j),&
		Au%gc1(i,j),Au%gc2(i,j)
		enddo
		write(30,*)
	    enddo
	endif
	if(myrank.eq.0)then
	    print *,"sigmam=",sigmam
	    print *,"npmls=",npmls
	    print *,"Nx=",Nx,"Ny=",Ny
	    print *,"N_x=",N_x,"N_y=",N_y
	    print *,"ds=",ds
	    print *,"dt=",dt
	    print *,"dn=",dn
	    print *,"Nt=",Nt
	    print *,"Nt/dn=",Nt/dn
	    print *,"R0=",R0!(-(4.d0*pi*ds*((order**npmly)-1.d0)/(c0*log(order)*dt*Nt)))
	    print *,"lambda0=",lambda0,floor(lambda0/ds)
	    print *,"lambda1=",lambda1,floor(lambda1/ds)
	    print *,"source_sposy=",sposy
	    print *,"Au:xi",Au%xinit,"Au:xf",Au%xend
	    print *,"Au:yi",Au%yinit,"Au:yf",Au%yend
	    print *,"Si:xi",Si%xinit,"Si:xf",Si%xend
	    print *,"Si:yi",Si%yinit,"Si:yf",Si%yend
	    print *,"Si:ep",Si%constant_ep,dsqrt(Si%constant_ep)
	    print *,"air:xi",air%xinit,"air:xf",air%xend
	    print *,"air:yi",air%yinit,"air:yf",air%yend
	    print *,"air:ep",air%constant_ep,dsqrt(air%constant_ep)
	    print *,"time4incident=",time4incident
	    print *,"distance4incident=",distance4incident
	    print *,"vel4incident=",vel4incident,vel4incident/c0
	    print *,"tfinStop=",tfinStop
	    print *,"xBDlow=",mypmly%xBDlow,"xBDhigh=",mypmly%xBDhigh
	    print *,"yBDlow=",mypmly%yBDlow,"yBDhigh=",mypmly%yBDhigh
	    print *,"ic=",ic,"jc=",jc
	    print *,"i_c=",i_c,"j_c=",j_c
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
       inRange(absIyrange,ptiny))then
    tmpFilename="transforms/r"//trim(adjustl(rankstr))//"infld.dat"
        open(unit=2,file=tmpFilename,status="unknown")
    endif

    if(inRange(absIxrange,ptrx).and.&
       inRange(absIyrange,ptry))then
        tmpFilename="transforms/r"//trim(adjustl(rankstr))//"rfld.dat"
        open(unit=3,file=tmpFilename,status="unknown")
    endif

    if(inRange(absIxrange,pttx).and.&
       inRange(absIyrange,ptty))then
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
		      locIyrange%init,locIyrange%ends)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	source definition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(inRange(absIyrange,sposy))then
	    if(n.eq.1)then
		print *,"source:",myrank
	    endif

!	plane wave
	    gauss_t=dexp(-0.5d0*((t-t0)/st)**2)*dsin(w0*t)
		do i=1,n_x,1
    Dx(i,sposy-absIyrange%init+1)=Dx(i,sposy-absIyrange%init+1)+gauss_t
		enddo
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Calculate E from D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call air%update_E(pex,pey,pdx,pdy)
        call Au%update_E(pex,pey,pdx,pdy)
	call SiO2%update_E(pex,pey,pdx,pdy)
        call Si%update_E(pex,pey,pdx,pdy)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Calculate E field in PML layers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call mypmly%updateE(PEx,PHz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	communication of the E-fields components at the subspace boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    call interchange_E()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Calculate H field inside simulation space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call updateH(locIxrange%init,locIxrange%ends,&
		      locIyrange%init,locIyrange%ends)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Calculate H field in PML layers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call mypmly%updateH(PHz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	communication of the H-fields components at the subspace boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call interchange_H()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Calculate TFs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(inRange(absIyrange,ptiny).and.&
    n.lt.tfinStop)then
	if(inRange(absIxrange,ptinx))then
	    write(2,*)n*dt,&
    ex(ptinx-absIxrange%init+1,ptiny-absIyrange%init+1),&
    hz(ptinx-absIxrange%init+1,ptiny-absIyrange%init+1)
	endif
        if(n.eq.1)then
	print *,"TFin: myrank=",myrank,ptiny-absIyrange%init+1
        endif
	do m=1,freqsteps,1
	    do i=1,n_x,1
real_exin(i,m)=real_exin(i,m)+dcos(arg(m)*t)*ex(i,ptiny-absIyrange%init+1)*dt
 img_exin(i,m)= img_exin(i,m)+dsin(arg(m)*t)*ex(i,ptiny-absIyrange%init+1)*dt

real_hzin(i,m)=real_hzin(i,m)+dcos(arg(m)*t)*hz(i,ptiny-absIyrange%init+1)*dt
 img_hzin(i,m)= img_hzin(i,m)+dsin(arg(m)*t)*hz(i,ptiny-absIyrange%init+1)*dt
	    enddo
	enddo
    endif

    if(inRange(absIyrange,ptry).and.&
    n.gt.tfinStop)then
        if(inRange(absIxrange,ptrx))then
	write(3,'(E30.19,4F30.19)')n*dt,&
    ex(ptrx-absIxrange%init+1,ptry-absIyrange%init+1),&
    hz(ptrx-absIxrange%init+1,ptry-absIyrange%init+1)
        endif
        if(n.eq.1)then
        print *,"TFr: myrank=",myrank,ptry-absIyrange%init+1
        endif
	do m=1,freqsteps,1
	    do i=1,n_x,1
real_exr(i,m)=real_exr(i,m)+dcos(arg(m)*t)*ex(i,ptry-absIyrange%init+1)*dt
 img_exr(i,m)= img_exr(i,m)+dsin(arg(m)*t)*ex(i,ptry-absIyrange%init+1)*dt

real_hzr(i,m)=real_hzr(i,m)+dcos(arg(m)*t)*hz(i,ptry-absIyrange%init+1)*dt
 img_hzr(i,m)= img_hzr(i,m)+dsin(arg(m)*t)*hz(i,ptry-absIyrange%init+1)*dt
	    enddo
	enddo
    endif

    if(inRange(absIyrange,ptty))then
        if(inRange(absIxrange,pttx))then
	write(4,'(E30.19,4F30.19)')n*dt,&
    ex(pttx-absIxrange%init+1,ptty-absIyrange%init+1),&
    hz(pttx-absIxrange%init+1,ptty-absIyrange%init+1)
        endif
        if(n.eq.1)then
        print *,"TFt: myrank=",myrank,ptty-absIyrange%init+1
        endif
	do m=1,freqsteps,1
	    do i=1,n_x,1
real_ext(i,m)=real_ext(i,m)+dcos(arg(m)*t)*ex(i,ptty-absIyrange%init+1)*dt
 img_ext(i,m)= img_ext(i,m)+dsin(arg(m)*t)*ex(i,ptty-absIyrange%init+1)*dt

real_hzt(i,m)=real_hzt(i,m)+dcos(arg(m)*t)*hz(i,ptty-absIyrange%init+1)*dt
 img_hzt(i,m)= img_hzt(i,m)+dsin(arg(m)*t)*hz(i,ptty-absIyrange%init+1)*dt
	    enddo
	enddo
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	imprime muestra temporal
	if(nout.eq.n)then
	    nout=nout+dn
	    frame=frame+1
	    write(framestr,'(I4)') frame

	    if(inRange(absIxrange,ic).and.&
		printTimesLine)then
	tmpFilename="timesLiney/r"//trim(adjustl(rankstr))&
	    //"E_time"//trim(adjustl(framestr))//".dat"
	    open(unit=400+frame,file=tmpFilename,status="unknown")
		do j=1,n_y,1
		    write(400+frame,*)&
	    (j+absIyrange%init-1)*ds/nano,&
	    Ex(ic-absIxrange%init+1,j),&
	    Hz(ic-absIxrange%init+1,j)
		enddo
	    endif

	    if(printTimesMap)then
	tmpFilename="timesGnuplotmap/r"//trim(adjustl(rankstr))&
	    //"E_time"//trim(adjustl(framestr))//".dat"
	    open(unit=500+frame,file=tmpFilename,status="unknown")

		do j=1,N_y,pf
		    do i=1,N_x,pf
	write(500+frame,*)&
	    (i+absIxrange%init-1)*1.d0,&
	    (j+absIyrange%init-1)*1.d0,&
	    Hz(i,j)
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
	    do i=1,n_x,1
		Syin=Syin-(&
	    (real_hzin(i,m)*real_exin(i,m))+&
	    ( img_hzin(i,m)* img_exin(i,m))&
	    )*ds
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
	    do i=1,n_x,1
		Syr=Syr-(&
	    (real_hzr(i,m)*real_exr(i,m))+&
	    ( img_hzr(i,m)* img_exr(i,m))&
	    )*ds
    write(21,'(2F15.5,2E20.10)')lambda/nano,i*ds/nano,&
	real_exr(i,m)*real_hzr(i,m),img_hzr(i,m)*img_exr(i,m)
	enddo
    write(11,*)lambda/nano,Syr
    write(21,*)
        enddo
    endif

    if(inRange(absIyrange,ptty).and.&
       n.gt.1)then
        tmpFilename="transforms/r"//trim(adjustl(rankstr))//"Syt.dat"
        open(unit=12,file=tmpFilename,status="unknown")
        do m=1,freqsteps,1
        Syt=0.d0
	    do i=1,n_x,1
		Syt=Syt-(&
	    (real_hzt(i,m)*real_ext(i,m))+&
	    ( img_hzt(i,m)* img_ext(i,m))&
	    )*ds
	enddo
    lambda=2.d0*pi*c0/arg(m)
    write(12,*)lambda/nano,Syt
        enddo
    endif

	call mpi_finalize(ierror)

	contains
	include 'include/myfdtd_2Dfunc.f90'
	include 'include/mympi_2Dfunc.f90'

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

	end program my3dfdtd
