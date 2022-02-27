!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c		calculate E field in vaccum
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updateE_vaccum(x0,x1,y0,y1)
	integer, intent(in) :: x0,x1,y0,y1
	    do i=x0,x1,1
		do j=y0,y1,1
		    Ex(i,j)=Dx(i,j)
		    Ey(i,j)=Dy(i,j)
		end do
	    end do
	end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c		update ecuation for D field
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updateD(x0,x1,y0,y1)
	integer, intent(in) :: x0,x1,y0,y1
	real*8 :: curlhx,curlhy
	    do i=x0,x1,1
		do j=y0,y1,1
!		    curlhx=hz(i+1,j)-hz(i,j)
!		    curlhy=hz(i,j+1)-hz(i,j)
		    curlhx=hz(i,j)-hz(i-1,j)
		    curlhy=hz(i,j)-hz(i,j-1)

		    dx(i,j)=dx(i,j)+Ku*curlhy
		    dy(i,j)=dy(i,j)-Ku*curlhx
		enddo
	    enddo
	end subroutine updateD
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c		update ecuation for H field into pmls
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updateH(x0,x1,y0,y1)
	integer, intent(in) :: x0,x1,y0,y1
	real*8 :: curlE
	    do i=x0,x1,1
		do j=y0,y1,1
!		    curlE=ey(i,j)-ey(i-1,j)-ex(i,j)+ex(i,j-1)
		    curlE=ey(i+1,j)-ey(i,j)-ex(i,j+1)+ex(i,j)
		    hz(i,j)=hz(i,j)-Ku*curlE
		enddo
	    enddo
	end subroutine updateH