    program test
    use constants_mod
    use mydatatypes
    use BERpmly_ExHz_mod
    use simpleDielectric_mod
    implicit none
    include 'mpif.h'
    character(len=4) :: framestr
    character(len=100) :: tmpFilename,rankstr,tmpstr,msg
    logical, parameter :: printLine=.false. ! if true, print 20 time frames of Ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	mpi variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: tagE_f,tagE_b=2,tagH_f=3,tagH_b=4
    integer :: nproces,ncolumn,ierror,myrank,status(mpi_status_size)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	simulation variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real*8, parameter :: lambda0=633.d0*nano
    real*8, parameter :: f0=c0/lambda0
    real*8, parameter :: ds=1.d0*nano
    real*8, parameter :: Ku=0.99d0
    real*8, parameter :: dt=Ku*ds/c0
    integer, parameter :: lambda0ds=floor(lambda0/ds)
    real*8, parameter :: yreal=20.d0*lambda0+344.d0*nano
    integer, parameter :: npos=floor(yreal/ds)
    integer, parameter :: NCmax=4 ! Number of cores
    integer, parameter :: n_pos=npos/NCmax !Number of points per core
    integer, parameter :: jc=0.5*npos ! center of sim. space
    integer, parameter :: j_c=0.5*n_pos
    real*8 :: t
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	PML variables
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    type(pmly_ExHz):: mypmly,mypmly_inc
    real*8 :: sigmam,order,R0
    integer, parameter :: npml=10
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	field variables
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    real*8, pointer, dimension(:) :: pdx,pex,phz
    real*8, target, dimension(0:n_pos+1) :: ex,hz
    real*8, target, dimension(n_pos) :: dx
!	incident field variables
    real*8, target, dimension(0:n_pos+1) :: ex_inc,hz_inc
    real*8, pointer, dimension(:) :: pex_inc,phz_inc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	material variables
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    type(Irange) :: absIrange,locIrange
    type(ConstantDielectric) :: air,SiO2,Si
    real*8, parameter :: air_n=1.d0
    real*8, parameter :: SiO2_n=1.46d0
    real*8, parameter :: Si_n=3.772d0
    real*8, parameter :: air_ep=air_n**2
    real*8, parameter :: SiO2_ep=SiO2_n**2
    real*8, parameter :: Si_ep=Si_n**2
    real*8, parameter :: z0=dsqrt(mu0/ep0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		fourier Transform variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, parameter :: freqsteps=100

    real*8, parameter :: linit=200d-9
    real*8, parameter :: lend=1000d-9
    real*8, parameter :: finit=c0/lend
    real*8, parameter :: fend=c0/linit
    real*8, parameter :: df=(fend-finit)/freqsteps

    real*8,dimension(freqsteps) :: pumpreal_hzin,pumpimg_hzin,pumpreal_hzr,pumpimg_hzr,&
	pumpreal_exin,pumpimg_exin,pumpreal_exr,pumpimg_exr

    real*8,dimension(freqsteps) :: Rreal_hzin,Rimg_hzin,Rreal_hzr,Rimg_hzr,&
	Rreal_exin,Rimg_exin,Rreal_exr,Rimg_exr

    real*8,dimension(freqsteps) :: arg
    real*8 :: pump_Syin,R_Syr,freq,lambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	source variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real*8, parameter :: periodoin=1.d0/f0
    real*8, parameter :: w0=2.d0*pi*f0
    real*8, parameter :: beta=2.d0
    real*8, parameter :: st=beta*periodoin
    real*8, parameter :: t0=4.d0*st
    integer, parameter :: spos=npml+floor(10.d0*nano/ds)
    integer, parameter :: TFSF_bd=spos+floor(5.d0*nano/ds)
!	incidence monitor
    integer, parameter :: monitor_pti=TFSF_bd-floor(5.d0*nano/ds)
!	reflection monitor
    integer, parameter :: monitor_ptr=TFSF_bd+floor(5.d0*nano/ds)
    real*8 :: gauss_t
!        10   r     s  i   B      SiO2
!|--npml--|---|-----|--|---|------|--------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	time duration variables
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    real*8 :: time4incident,distance4incident,vel4incident,timeSim
    integer :: n,m,j,dn,nout,frame,tfinStop,Nt
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	Fourier Transform initialization
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    pumpreal_exin=0.d0; pumpimg_exin=0.d0; pumpreal_exr=0.d0; pumpimg_exr=0.d0
    pumpreal_hzin=0.d0; pumpimg_hzin=0.d0; pumpreal_hzr=0.d0; pumpimg_hzr=0.d0

    Rreal_exin=0.d0; Rimg_exin=0.d0; Rreal_exr=0.d0; Rimg_exr=0.d0
    Rreal_hzin=0.d0; Rimg_hzin=0.d0; Rreal_hzr=0.d0; Rimg_hzr=0.d0

    freq=finit
    do m=1,freqsteps,1
	arg(m)=2.d0*pi*freq
	freq=freq+df
    enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	Field initialization
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ex=0.d0; hz=0.d0; dx=0.d0
    pdx=>Dx; pex=>Ex; phz=>Hz
    ex_inc=0.d0; hz_inc=0.d0
    pex_inc=>ex_inc; phz_inc=>hz_inc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Initialize MPI environment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call my_mpi_init()
    write(rankstr,'(I4)') myrank
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Initialize pml
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    order=3.d0
    R0=dexp(-16.d0)
    sigmam=(order+1.d0)*c0*log(1.d0/R0)/(2.d0*npml*ds)
    call mypmly%initialize(air_ep*ep0,Si_ep*ep0,sigmam,npml,order,dt,ds,myrank,n_pos,NCmax)
    call mypmly_inc%initialize(air_ep*ep0,air_ep*ep0,sigmam,npml,order,dt,ds,myrank,n_pos,NCmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	medium initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    air%posinit=npml+1
!    air%posend=air%posinit+floor(20.d0*nano/ds)-1
    air%posend=air%posinit+floor(20.d0*lambda0/ds)-1
    air%constant_ep=air_ep

    SiO2%posinit=air%posend+1
    SiO2%posend=SiO2%posinit+floor(284.d0*nano/ds)-1
    SiO2%constant_ep=SiO2_ep

    Si%posinit=SiO2%posend+1
    Si%posend=npos-npml
    Si%constant_ep=Si_ep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	time Window
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    distance4incident=ds*real(jc,8)!+beta*3.d0*lambda0/back_n
    vel4incident=c0/air_n
    time4incident=distance4incident/vel4incident
    tfinStop=floor(time4incident/dt)
    Nt=6*tfinStop+3*floor(t0/dt)
    dn=aint(Nt/20.d0)	! printing rate per frame
    nout=dn
    frame=0

    call air%initialize_material(absIrange,n_pos)
    call SiO2%initialize_material(absIrange,n_pos)
    call Si%initialize_material(absIrange,n_pos)

    if(myrank.eq.0)then
	tmpFilename="SimParameters.dat"
        open(unit=1,file=tmpFilename,status="unknown")
        write(1,*) "ds=",ds
        write(1,*) "beta=",beta

	print *,"c0=",c0,c0/1e8
	print *,"lambda0=",lambda0,lambda0/nano,"nm",floor(lambda0/ds),"cells"
	print *,"f0=",f0,c0/f0
	print *,"w0=",w0,2.d0*pi*c0/w0
	print *,"sigmam=",sigmam
        print *,"npml=",npml
        print *,"Npos=",Npos
        print *,"N_pos=",N_pos
        print *,"ds=",ds
        print *,"dt=",dt
        print *,"dn=",dn
        print *,"Nt=",Nt
        print *,"Nt/dn=",Nt/dn
        print *,"source_spos=",spos
        print *,"air:posi",air%posinit,"air:posf",air%posend
        print *,"air:ep",air%constant_ep,dsqrt(air%constant_ep)
        print *,"SiO2:posi",SiO2%posinit,"SiO2:posf",SiO2%posend
        print *,"SiO2:ep",SiO2%constant_ep,dsqrt(SiO2%constant_ep)
        print *,"Si:posi",Si%posinit,"Si:posf",Si%posend
	print *,"Si:ep",Si%constant_ep,dsqrt(Si%constant_ep)
        print *,"time4incident=",time4incident
        print *,"distance4incident=",distance4incident
        print *,"vel4incident=",vel4incident,vel4incident/c0
        print *,"tfinStop=",tfinStop
        print *,"yBDlow=",mypmly%yBDlow,"yBDhigh=",mypmly%yBDhigh
        print *,"monitor_pti=",monitor_pti,"monitor_ptr=",monitor_ptr
	print *,"spos=",spos
	print *,"TFSF_bd=",TFSF_bd
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		open file for input fields of the FT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(inRange(absIrange,monitor_pti))then
	tmpFilename="transforms/r"//trim(adjustl(rankstr))//"infld.dat"
        open(unit=2,file=tmpFilename,status="unknown")
    endif
    if(inRange(absIrange,monitor_ptr))then
	tmpFilename="transforms/r"//trim(adjustl(rankstr))//"rfld.dat"
        open(unit=3,file=tmpFilename,status="unknown")
    endif
    call mpi_barrier(mpi_comm_world,ierror)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Start FDTD time loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do n=1,Nt,1
!    do n=1,1,1
	t=real(n-1,8)*dt

	do j=locIrange%init,locIrange%ends,1
	    dx(j)=dx(j)+Ku*(hz(j)-hz(j-1))
	    ex_inc(j)=ex_inc(j)+Ku*(hz_inc(j)-hz_inc(j-1))
	enddo

!	source
	if(inRange(absIrange,spos))then
	    gauss_t=dexp(-0.5d0*((t-t0)/st)**2)*dsin(w0*t)
	    ex_inc(spos-absIrange%init+1)=ex_inc(spos-absIrange%init+1)+gauss_t
        endif
	call mypmly_inc%updateE(pex_inc)
        call interchange_E(pex_inc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Consistency condition at TFSF Boundary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(inRange(absIrange,TFSF_bd))then
	    dx(TFSF_bd-absIrange%init+1)=dx(TFSF_bd-absIrange%init+1)-Ku*hz_inc(TFSF_bd-absIrange%init)
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Calculate E from D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call air%update_E(pex,pdx)
	call SiO2%update_E(pex,pdx)
        call Si%update_E(pex,pdx)
	call mypmly%updateE(pex)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	communication of the E-fields components at the subspace boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call interchange_E(pex)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do j=locIrange%init,locIrange%ends,1
	    hz(j)=hz(j)+Ku*(ex(j+1)-ex(j))
	    hz_inc(j)=hz_inc(j)+Ku*(ex_inc(j+1)-ex_inc(j))
	enddo
	call mypmly_inc%updateH(phz_inc)
	call mypmly%updateH(phz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Consistency condition at TFSF Boundary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(inRange(absIrange,TFSF_bd))then
	    hz(TFSF_bd-absIrange%init)=hz(TFSF_bd-absIrange%init)-Ku*ex_inc(TFSF_bd-absIrange%init+1)
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	communication of the H-fields components at the subspace boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call interchange_H(phz_inc)
        call interchange_H(phz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Calculate TFs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(inRange(absIrange,monitor_pti))then
        write(2,'(3E30.19)')n*dt,&
	    ex_inc(monitor_pti-absIrange%init+1),hz_inc(monitor_pti-absIrange%init+1)
	do m=1,freqsteps,1
    pumpreal_exin(m)=pumpreal_exin(m)+dcos(arg(m)*t)*ex_inc(monitor_pti-absIrange%init+1)*dt
     pumpimg_exin(m)= pumpimg_exin(m)+dsin(arg(m)*t)*ex_inc(monitor_pti-absIrange%init+1)*dt
    pumpreal_hzin(m)=pumpreal_hzin(m)+dcos(arg(m)*t)*hz_inc(monitor_pti-absIrange%init+1)*dt
     pumpimg_hzin(m)= pumpimg_hzin(m)+dsin(arg(m)*t)*hz_inc(monitor_pti-absIrange%init+1)*dt
	enddo
    endif

    if(inRange(absIrange,monitor_ptr))then
        write(3,'(3E30.19)')n*dt,&
	    ex(monitor_ptr-absIrange%init+1),hz(monitor_ptr-absIrange%init+1)
	do m=1,freqsteps,1
    Rreal_exr(m)=Rreal_exr(m)+dcos(arg(m)*t)*ex(monitor_ptr-absIrange%init+1)*dt
     Rimg_exr(m)= Rimg_exr(m)+dsin(arg(m)*t)*ex(monitor_ptr-absIrange%init+1)*dt
    Rreal_hzr(m)=Rreal_hzr(m)+dcos(arg(m)*t)*hz(monitor_ptr-absIrange%init+1)*dt
     Rimg_hzr(m)= Rimg_hzr(m)+dsin(arg(m)*t)*hz(monitor_ptr-absIrange%init+1)*dt
	enddo
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	imprime muestra temporal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(nout.eq.n.and.printLine)then
        nout=nout+dn
        frame=frame+1
        write(framestr,'(I4)') frame
	tmpFilename="timesLine/r"//trim(adjustl(rankstr))//&
	"frame"//trim(adjustl(framestr))//".dat"
        open(unit=400+frame,file=tmpFilename,status="unknown")
	do j=1,n_pos,1
	    write(400+frame,*)(j+absIrange%init-1)*ds/nano,ex(j),hz(j)
	enddo
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	End FDTD time loop
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		print TFs to file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(inRange(absIrange,monitor_pti))then
        tmpFilename="transforms/r"//trim(adjustl(rankstr))//"pump_Syin.dat"
        open(unit=10,file=tmpFilename,status="unknown")
        do m=1,freqsteps,1
	    pump_Syin=(-(pumpreal_hzin(m)*pumpreal_exin(m))-(pumpimg_hzin(m)*pumpimg_exin(m)))*ds*z0

	    lambda=2.d0*pi*c0/arg(m)
	    write(10,'(4E30.19)')lambda/nano,pump_Syin
        enddo
    endif

    if(inRange(absIrange,monitor_ptr))then
        tmpFilename="transforms/r"//trim(adjustl(rankstr))//"R_Syr.dat"
        open(unit=16,file=tmpFilename,status="unknown")
        do m=1,freqsteps,1
	    R_Syr=(-(Rreal_hzr(m)*Rreal_exr(m))-(Rimg_hzr(m)*Rimg_exr(m)))*ds*z0

	    lambda=2.d0*pi*c0/arg(m)
	    write(16,'(4E30.19)')lambda/nano,R_Syr
	enddo
    endif
    call mpi_finalize(ierror)
contains
    include 'include/mympi_1Dfunc.f90'

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

    end program test
