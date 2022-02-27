    module constants_mod
	implicit none
	real*8,parameter :: pi=4.d0*datan(1.d0)
! physics constants
	real*8,parameter :: ep0=8.85418781762039d-12 ![in C^2/(Nxm^2)]
	real*8,parameter :: mu0=pi*4d-7 ![in H/m]
	real*8,parameter :: c0=1.d0/dsqrt(ep0*mu0)
!	   real*8,parameter :: h=4.135667662d-15 ![in eV]
	real*8,parameter :: h=6.62607015d-34![in Jxs]
	real*8,parameter :: hbar=h/(2.d0*pi)!1.054571817d-34 
	real*8,parameter :: boltzman_const=1.380649d-23 ![in J/K]
	real*8,parameter :: echarge=1.602176634d-19 ![in C]
! convertion factors
	real*8,parameter :: J2eV=6.241509d18 ! #A Joules= #B J2eV
	real*8,parameter :: eV2J=1.d0/J2eV	  ! #B Ev= #A J2eV
! prefix constants
	real*8,parameter :: nano=1.d-9
	real*8,parameter :: micro=1.d-6
	real*8,parameter :: femto=1d-15
	real*8,parameter :: Tera=1d12
	real*8,parameter :: Peta=1d15
! complex constants
	complex*16,parameter :: unoi=dcmplx(0.d0,1.d0)
	complex*16,parameter :: zero=dcmplx(0.d0,0.d0)

    contains 
    subroutine show_consts()
	print *,"pi=",pi
	print *,"ep0",ep0
	print *,"mu0",mu0
	print *,"c0",c0
	print *,"h",h
	print *,"hbar",hbar
	print *,"nano",nano
	print *,"femto",femto
	print *,"Tera",Tera
	print *,"Peta",Peta
	print *,"zero=",zero
	print *,"unoi=",unoi
	end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		create vtk header
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine makeVTKheader(fileid,msg,N_x,N_y,N_z)
        integer, intent(in) :: fileid,N_x,N_y,N_z
    character(*) :: msg
!    print *,"fileid=",fileid
!    print *,"msg=",msg
    write(fileid,"(a)") "# vtk DataFile Version 2.0"
    write(fileid,"(a)") msg
    write(fileid,"(a)") "ASCII"
    write(fileid,"(a)") "DATASET STRUCTURED_POINTS"
    write(fileid,"(a,5I5)") "DIMENSIONS", N_x,N_y,N_z
    write(fileid,"(a)") "ASPECT_RATIO 1 1 1"
    write(fileid,"(a)") "ORIGIN 0 0 0"
    write(fileid,"(a,I10)") "POINT_DATA",N_x*N_y*N_z
    write(fileid,"(a)") "SCALARS volume_Ez double 1"
    write(fileid,"(a)") "LOOKUP_TABLE default"

    end subroutine makeVTKheader
    end module 
