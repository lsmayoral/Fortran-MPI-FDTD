	module my_matfunctions_mod
	use constants_mod
	implicit none
	contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Drude+2CP model for GOLD dielectric function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    function gold_ep(w) result(epg)
		real*8, intent(in) :: w
		complex*16 :: epg,critpnts
		integer :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		Parameter for Gold was taken from
!		Appl Phys A (2011) 103: 849-853
!		DOI 10.1007/s00339-010-6224-9
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	real*8, parameter :: ep=1.d0,wp=5.76d15,rp=wp/100.d0
		real*8, parameter :: ep=1.1431d0,rp=1.0805d14,wp=1.3202d16
		real*8, parameter, dimension(2) :: al=(/ 0.26698d0,3.0834d0 /)
		real*8, parameter, dimension(2) :: phil=(/ -1.2371d0,-1.0968d0 /)
		real*8, parameter, dimension(2) :: oml=(/ 3.8711d15,4.1684d15 /)
		real*8, parameter, dimension(2) :: rl=(/ 4.4642d14,2.4555d15 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		critpnts=zero
		    do i=1,2,1
			critpnts=critpnts+al(i)*oml(i)*&
        ((cdexp(unoi*phil(i))/(oml(i)-w-unoi*rl(i)))+&
        (cdexp(-unoi*phil(i))/(oml(i)+w+unoi*rl(i))))
		    enddo
		epg=ep-((wp*wp)/(w*w+unoi*rp*w))+critpnts
	    end function gold_ep

	end module my_matfunctions_mod
