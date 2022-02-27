!c*********************************************
!	set consitency indexes
!*********************************************
	subroutine set_consistency_indexes(absrangesx,absrangesy)
	type(Irange), intent(in) :: absrangesx,absrangesy
	integer :: xi,xf,yi,yf
	xi=absrangesx%init
	xf=absrangesx%ends
	yi=absrangesy%init
	yf=absrangesy%ends
	consis_xi=0
	consis_xf=-1
	if(xi.le.ib.and.xf.ge.ia)then
	    consis_xi=1
	    consis_xf=n_x
	    if(xi.le.ia)then
		consis_xi=ia-xi+1
	    end if
	    if(xf.ge.ib)then
		consis_xf=ib-xi+1
	    end if
	end if
	consis_yi=0
	consis_yf=-1
	if(yi.le.jb.and.yf.ge.ja)then
	    consis_yi=1
	    consis_yf=n_y
	    if(yi.le.ja)then
		consis_yi=ja-yi+1
	    end if
	    if(yf.ge.jb)then
		consis_yf=jb-yi+1
	    end if
	end if
!	print *,"consis_x:",myrank,consis_xi,consis_xf
!c	print *,"consis_y:",myrank,consis_yi,consis_yf
	end subroutine
!c*********************************************
!c	update hz_inc
!c*********************************************
	subroutine update_hzinc(y0,y1)
	integer, intent(in) :: y0,y1
	integer :: j
	    do j=y0,y1,1
		hz_inc(j)=hz_inc(j)+Ku*(ex_inc(j)-ex_inc(j-1))
	    end do
	end subroutine
!c*********************************************************
!c	update incident field 
!c*********************************************************
	subroutine update_exinc(y0,y1)
	integer, intent(in) :: y0,y1
	integer :: j
	real*8 :: Kuk
	    do j=y0,y1,1
		Kuk=Ku/ga_inc(j)
		ex_inc(j)=ex_inc(j)+Kuk*(hz_inc(j+1)-hz_inc(j))
	    end do
	end subroutine
!c*********************************************
!c		Boundary condition for the 1D incident buffer
!c*********************************************
	subroutine hzinc_boundaryCond()
	    if(ycord.eq.0)then
		hz_inc(1)=hz_inc_low_m2
		hz_inc_low_m2=hz_inc_low_m1
		hz_inc_low_m1=hz_inc(2)
	    end if
	    if(ycord.eq.(Ymax-1))then
		hz_inc(n_y-1)=hz_inc_high_m2
		hz_inc_high_m2=hz_inc_high_m1
		hz_inc_high_m1=hz_inc(n_y-2)
	    end if
	end subroutine
!c*********************************************************
!c	Consistency conditions for Hz
!c*********************************************************
	subroutine consistCond_hz(absrangesx,absrangesy)
	type(Irange), intent(in) :: absrangesx,absrangesy
	integer :: xi,xf,yi,yf
	xi=absrangesx%init
	xf=absrangesx%ends
	yi=absrangesy%init
	yf=absrangesy%ends
	    if(yi.le.ja.and.yf.ge.ja)then
		do i=consis_xi,consis_xf,1
		    hz(i,ja-yi+1)=hz(i,ja-yi+1)+Ku*ex_inc(ja-yi)
		end do
	    end if
	    if(yi.le.jb.and.yf.ge.jb)then
		do i=consis_xi,consis_xf,1
		    hz(i,jb-yi+1)=hz(i,jb-yi+1)-Ku*ex_inc(jb-yi+1)
		end do
	    end if
	end subroutine
!c*********************************************************
!c	consistency conditions for Dx
!c*********************************************************
	subroutine consistCond_dx(absrangesx,absrangesy)
	type(Irange), intent(in) :: absrangesx,absrangesy
	integer :: xi,xf,yi,yf
	xi=absrangesx%init
	xf=absrangesx%ends
	yi=absrangesy%init
	yf=absrangesy%ends
	    if(yi.le.ja.and.yf.ge.ja)then
		do i=consis_xi,consis_xf,1
		    dx(i,ja-yi)=dx(i,ja-yi)+Ku*hz_inc(ja-yi+1)
		end do
	    end if
	    if(yi.le.jb.and.yf.ge.jb)then
		do i=consis_xi,consis_xf,1
		    dx(i,jb-yi+1)=dx(i,jb-yi+1)-Ku*hz_inc(jb-yi+1)
		end do
	    end if
	end subroutine
!c*********************************************************
!c	consistency conditions for Dy
!c*********************************************************
	subroutine consistCond_dy(absrangesx,absrangesy)
	type(Irange), intent(in) :: absrangesx,absrangesy
	integer :: xi,xf,yi,yf
	xi=absrangesx%init
	xf=absrangesx%ends
	yi=absrangesy%init
	yf=absrangesy%ends
	    if(xi.le.ia.and.xf.ge.ia.and.ia.gt.1)then
		do j=consis_yi,consis_yf,1
		    dy(ia-xi,j)=dy(ia-xi,j)-Ku*hz_inc(j)
		end do
	    end if
	    if(xi.le.ib.and.xf.ge.ib.and.ib.lt.Nx)then
		do j=consis_yi,consis_yf,1
		    dy(ib-xi+1,j)=dy(ib-xi+1,j)+Ku*hz_inc(j)
		end do
	    end if
	end subroutine
!c*********************************************
!c	synchronization for the Hz_inc field
!c*********************************************
	subroutine synchonize_hzinc()
	    if(ycord.gt.0)then
!c	    if(n.eq.1)then
!c		print *,"send(H1):",myrank
!c	    end if
	call mpi_send(hz_inc(1),1,mpi_double_precision,neighborW,&
		tagCH1_WE,mpi_comm_world,ierror)
	    endif
	    if(ycord.lt.(Ymax-1))then
!c	    if(n.eq.1)then
!c		print *,"recv(H1):",myrank
!c	    end if
	call mpi_recv(hz_inc(n_y+1),1,mpi_double_precision,&
		neighborE,tagCH1_WE,mpi_comm_world,status,ierror)
	    end if
	    if(ycord.lt.(Ymax-1)) then
!c	    if(n.eq.1)then
!c		print *,"send(H2):",myrank
!c	    end if
	call mpi_send(hz_inc(n_y),1,mpi_double_precision,neighborE,&
		tagCH2_WE,mpi_comm_world,ierror)
	    endif
	    if(ycord.gt.0) then
!c	    if(n.eq.1)then
!c		print *,"recv(H2):",myrank
!c	    end if
	call mpi_recv(hz_inc(0),1,mpi_double_precision,&
		neighborW,tagCH2_WE,mpi_comm_world,status,ierror)
	    endif
	end subroutine

!c*********************************************
!c	synchronization for the ex_inc field
!c*********************************************
	subroutine synchonize_exinc()
	    if(ycord.gt.0)then
!c	    if(n.eq.1)then
!c		print *,"send(Ex1):",myrank
!c	    end if
	call mpi_send(ex_inc(1),1,mpi_double_precision,neighborW,&
		tagCE1_WE,mpi_comm_world,ierror)
	    endif
	    if(ycord.lt.(Ymax-1))then
!c	    if(n.eq.1)then
!c		print *,"recv(Ex1):",myrank
!c	    end if
	call mpi_recv(ex_inc(n_y+1),1,mpi_double_precision,&
		neighborE,tagCE1_WE,mpi_comm_world,status,ierror)
	    end if
	    if(ycord.lt.(Ymax-1)) then
!c	    if(n.eq.1)then
!c		print *,"send(Ex2):",myrank
!c	    end if
	call mpi_send(ex_inc(n_y),1,mpi_double_precision,neighborE,&
		tagCE2_WE,mpi_comm_world,ierror)
	endif
	    if(ycord.gt.0) then
!c	    if(n.eq.1)then
!c		print *,"recv(Ex2):",myrank
!c	    end if
	call mpi_recv(ex_inc(0),1,mpi_double_precision,&
		neighborW,tagCE2_WE,mpi_comm_world,status,ierror)
	    endif
	end subroutine