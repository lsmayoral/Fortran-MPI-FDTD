!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		calculate E field in vaccum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine updateE_vaccum(x0,x1,y0,y1,z0,z1)
    integer, intent(in) :: x0,x1,y0,y1,z0,z1
	do k=z0,z1,1
	    do j=y0,y1,1
		do i=x0,x1,1
		    Ex(i,j,k)=Dx(i,j,k)
		    Ey(i,j,k)=Dy(i,j,k)
		    Ez(i,j,k)=Dz(i,j,k)
		enddo
	    enddo
        enddo
    end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		update ecuation for D field into vacuum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine updateD(x0,x1,y0,y1,z0,z1)
	integer, intent(in) :: x0,x1,y0,y1,z0,z1
	real*8 :: curlhzy,curlhyz,curlhxz,curlhzx,curlhyx,curlhxy
	    do k=z0,z1,1
		do j=y0,y1,1
		    do i=x0,x1,1
	curlhzy=hz(i,j,k)-hz(i,j-1,k)
	curlhyz=-hy(i,j,k)+hy(i,j,k-1)
	curlhxz=hx(i,j,k)-hx(i,j,k-1)
	curlhzx=-hz(i,j,k)+hz(i-1,j,k)
	curlhyx=hy(i,j,k)-hy(i-1,j,k)
	curlhxy=-hx(i,j,k)+hx(i,j-1,k)

	    dx(i,j,k)=dx(i,j,k)+Ku*(curlhzy+curlhyz)
	    dy(i,j,k)=dy(i,j,k)+Ku*(curlhxz+curlhzx)
	    dz(i,j,k)=dz(i,j,k)+Ku*(curlhyx+curlhxy)
!	    ex(i,j,k)=ex(i,j,k)+Ku*(curlhzy+curlhyz)
!	    ey(i,j,k)=ey(i,j,k)+Ku*(curlhxz+curlhzx)
!	    ez(i,j,k)=ez(i,j,k)+Ku*(curlhyx+curlhxy)
		    enddo
		enddo
	    enddo
	end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		update ecuation for H field into vacuum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine updateH(x0,x1,y0,y1,z0,z1)
	integer, intent(in) :: x0,x1,y0,y1,z0,z1
	real*8 :: curleyz,curlezy,curlezx,curlexz,curlexy,curleyx
	    do k=z0,z1,1
		do j=y0,y1,1
		    do i=x0,x1,1
	curleyz=ey(i,j,k+1)-ey(i,j,k)
	curlezy=-ez(i,j+1,k)+ez(i,j,k)
	curlezx=ez(i+1,j,k)-ez(i,j,k)
	curlexz=-ex(i,j,k+1)+ex(i,j,k)
	curlexy=ex(i,j+1,k)-ex(i,j,k)
	curleyx=-ey(i+1,j,k)+ey(i,j,k)

	    hx(i,j,k)=hx(i,j,k)+Ku*(curleyz+curlezy)
	    hy(i,j,k)=hy(i,j,k)+Ku*(curlezx+curlexz)
	    hz(i,j,k)=hz(i,j,k)+Ku*(curlexy+curleyx)
		    enddo
		enddo
	    enddo
	end subroutine
