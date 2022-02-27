subroutine updateEs((xi,xf,yi,yf,zi,zf)
        integer, intent(in) :: xi,xf,yi,yf,zi,zf
        real*8 :: curlHzy,curlHyz,curlHxz,curlHzx,curlHyx,curlHxy
        integer :: i,j,k
        do i=xi,xf,1
	    do j=yi,yf,1
		do k=ki,kf,1
curlHzy=Hzx(i,  j,  k)  +Hzy(i,  j,  k)  -Hzx(i,  j-1,k)  -Hzy(i,  j-1,k)
curlHyz=Hyz(i,  j,  k-1)+Hyx(i,  j,  k-1)-Hyz(i,  j,  k)  -Hyx(i,  j,  k)
curlHxz=Hxy(i,  j,  k)  +Hxz(i,  j,  k)  -Hxy(i,  j,  k-1)-Hxz(i,  j,  k-1)
curlHzx=Hzx(i-1,j,  k)  +Hzy(i-1,j,  k)  -Hzx(i,  j,  k)  -Hzy(i,  j,  k)
curlHyx=Hyz(i,  j,  k)  +Hyx(i,  j,  k)  -Hyz(i-1,j,  k)  -Hyx(i-1,j,  k)
curlHxy=Hxy(i,  j-1,k)  +Hxz(i,  j-1,k)  -Hxy(i,  j,  k)  -Hxz(i,  j,  k)

	Exy(i,j,k)=Exy(i,j,k)+Ku*curlHzy
	Exz(i,j,k)=Exz(i,j,k)+Ku*curlHyz
	Eyz(i,j,k)=Eyz(i,j,k)+Ku*curlHxz
	Eyx(i,j,k)=Eyx(i,j,k)+Ku*curlHzx
	Ezx(i,j,k)=Ezx(i,j,k)+Ku*curlHyx
	Ezy(i,j,k)=Ezy(i,j,k)+Ku*curlHxy
		enddo
	    enddo
        enddo
    end subroutine updateEs

function totalEx(x,y,z) result(Ex_)
    integer, intent(in) :: x,y,z
    Ex_=Exy(x,y,z)+Exz(x,y,z)
end function totalEx

function totalEy(x,y,z) result(Ey_)
    integer, intent(in) :: x,y,z
    Ey_=Eyz(x,y,z)+Eyx(x,y,z)
end function totalEy

function totalEz(x,y,z) result(Ez_)
    integer, intent(in) :: x,y,z
    Ez_=Ezx(x,y,z)+Ezy(x,y,z)
end function totalEz

subroutine updateH(xi,xf,yi,yf,zi,zf)
        integer, intent(in) :: xi,xf,yi,yf,zi,zf
        real*8 :: curlEzy,curlEyz,curlExz,curlEzx,curlEyx,curlExy
        integer :: i,j,k
        do i=xi,xf,1
	    do j=yi,yf,1
		do k=ki,kf,1
!curlEzy=Ezx(i,  j+1,k)  +Ezy(i,  j+1,k)  -Ezx(i,  j,  k)  -Ezy(i,  j,  k)
!curlEyz=Eyz(i,  j,  k)  +Eyx(i,  j,  k)  -Eyz(i,  j,  k+1)-Eyx(i,  j,  k+1)
!curlExz=Exy(i,  j,  k+1)+Exz(i,  j,  k+1)-Exy(i,  j,  k)  -Exz(i,  j,  k)
!curlEzx=Ezx(i,  j,  k)  +Ezy(i,  j,  k)  -Ezx(i+1,j,  k)  -Ezy(i+1,j,  k)
!curlEyx=Eyz(i+1,j,  k)  +Eyx(i+1,j,  k)  -Eyz(i,  j,  k)  -Eyx(i,  j,  k)
!curlExy=Exy(i,  j,  k)  +Exz(i,  j,  k)  -Exy(i,  j+1,k)  -Exz(i,  j+1,k)

curlEzy=Ezx(i,  j,  k)  +Ezy(i,  j,  k)  -Ezx(i,  j-1,k)  -Ezy(i,  j-1,k)
curlEyz=Eyz(i,  j,  k-1)+Eyx(i,  j,  k-1)-Eyz(i,  j,  k)  -Eyx(i,  j,  k)
curlExz=Exy(i,  j,  k)  +Exz(i,  j,  k)  -Exy(i,  j,  k-1)-Exz(i,  j,  k-1)
curlEzx=Ezx(i-1,j,  k)  +Ezy(i-1,j,  k)  -Ezx(i,  j,  k)  -Ezy(i,  j,  k)
curlEyx=Eyz(i,  j,  k)  +Eyx(i,  j,  k)  -Eyz(i-1,j,  k)  -Eyx(i-1,j,  k)
curlExy=Exy(i,  j-1,k)  +Exz(i,  j-1,k)  -Exy(i,  j,  k)  -Exz(i,  j,  k)

	Hxy(i,j,k)=Hxy(i,j,k)+Ku*curlEzy
	Hxz(i,j,k)=Hxz(i,j,k)+Ku*curlEyz
	Hyz(i,j,k)=Hyz(i,j,k)+Ku*curlExz
	Hyx(i,j,k)=Hyx(i,j,k)+Ku*curlEzx
	Hzx(i,j,k)=Hzx(i,j,k)+Ku*curlEyx
	Hzy(i,j,k)=Hzy(i,j,k)+Ku*curlExy
		enddo
	    enddo
        enddo
    end subroutine updateHs
