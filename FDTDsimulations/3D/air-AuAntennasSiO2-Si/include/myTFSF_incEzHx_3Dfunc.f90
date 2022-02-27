!c*********************************************
!	set consitency indexes
!*********************************************
	subroutine set_consistency_indexes(absrangesx,absrangesy,absrangesz)
	type(Irange), intent(in) :: absrangesx,absrangesy,absrangesz
	integer :: xi,xf,yi,yf,zi,zf
	xi=absrangesx%init
	xf=absrangesx%ends
	yi=absrangesy%init
	yf=absrangesy%ends
	zi=absrangesz%init
	zf=absrangesz%ends

	consisD_xi=0
	consisD_xf=-1
	consisH_xi=0
	consisH_xf=-1
	if(xi.le.ib.and.xf.ge.ia)then
	    consisD_xi=1
	    consisD_xf=n_x
	    consisH_xi=1
	    consisH_xf=n_x
	    if(xi.le.ia)then
		consisD_xi=ia-xi+1
		consisH_xi=ia-xi+1
	    endif
	    if(xf.ge.ib)then
		consisD_xf=ib-xi+1
		consisH_xf=ib-xi+1
	    endif
	endif
	consisD_yi=0
	consisD_yf=-1
	consisH_yi=0
	consisH_yf=-1
	if(yi.le.jb.and.yf.ge.ja)then
	    consisD_yi=1
	    consisD_yf=n_y
	    consisH_yi=1
	    consisH_yf=n_y
	    if(yi.le.ja)then
		consisD_yi=ja-yi+1
		consisH_yi=ja-yi+1
	    endif
	    if(yf.ge.jb)then
		consisD_yf=jb-yi
		consisH_yf=jb-yi+1
	    endif
	endif
	consisD_zi=0
	consisD_zf=-1
	consisH_zi=0
	consisH_zf=-1
	if(zi.le.kb.and.zf.ge.ka)then
	    consisD_zi=1
	    consisD_zf=n_z
	    consisH_zi=1
	    consisH_zf=n_z
	    if(zi.le.ka)then
		consisD_zi=ka-zi+1
		consisH_zi=ka-zi+1
	    endif
	    if(zf.ge.kb)then
		consisD_zf=kb-zi+1
		consisH_zf=kb-zi+1
	    endif
	endif
	write(*,'(I2,A11,I3,A11,I3,A11,I3,A11,I3,A11,I3,A11,I3,A11,I3,A11,I3)') myrank,&
			"consisD_xi=",consisD_xi,"consisD_xf=",consisD_xf,&
			"consisD_yi=",consisD_yi,"consisD_yf=",consisD_yf,&
			"consisD_zi=",consisD_zi,"consisD_zf=",consisD_zf
	write(*,'(I2,A11,I3,A11,I3,A11,I3,A11,I3,A11,I3,A11,I3,A11,I3,A11,I3)') myrank,&
			"consisH_xi=",consisH_xi,"consisH_xf=",consisH_xf,&
			"consisH_yi=",consisH_yi,"consisH_yf=",consisH_yf,&
			"consisH_zi=",consisH_zi,"consisH_zf=",consisH_zf
	end subroutine
!c*********************************************
!c	update hx_inc (propagating along Y)
!c*********************************************
	subroutine update_hxincy(y0,y1)
	integer, intent(in) :: y0,y1
	integer :: j
	    do j=y0,y1,1
		hx_inc(j)=hx_inc(j)+Ku*(ez_inc(j)-ez_inc(j+1))
	    enddo
	end subroutine
!c*********************************************************
!c	update Ez field (propagating along Y)
!c*********************************************************
	subroutine update_ezincz(y0,y1)
	integer, intent(in) :: y0,y1
	integer :: j
	    do j=y0,y1,1
		ez_inc(j)=ez_inc(j)+tfsf_cb(j)*(hx_inc(j-1)-hx_inc(j))
	    enddo
	end subroutine
!c*********************************************
!c		Boundary condition for the 1D incident buffer
!c*********************************************
! for Ex source propagating along Z
	subroutine ezincy_boundaryCond()
	    if(ycord.eq.0)then
		ez_inc(2)=    ez_inc_low_m2
		ez_inc_low_m2=ez_inc_low_m1
		ez_inc_low_m1=ez_inc(3)
	    endif
	    if(ycord.eq.(Ymax-1))then
		ez_inc(n_y-2)=ez_inc_high_m2
		ez_inc_high_m2=ez_inc_high_m1
		ez_inc_high_m1=ez_inc(n_y-3)
	    endif
	end subroutine
!c*********************************************************
!c	consistency conditions for D
!c*********************************************************
!	Dy
	subroutine consistCond_dy(absrangesx,absrangesy,absrangesz)
	type(Irange), intent(in) :: absrangesx,absrangesy,absrangesz
	integer :: zi,zf,i,j
!	xi=absrangesx%init
!	xf=absrangesx%ends
!	yi=absrangesy%init
!	yf=absrangesy%ends
	zi=absrangesz%init
	zf=absrangesz%ends
	    if(zi.le.ka.and.zf.ge.ka)then
		do i=consisD_xi,consisD_xf,1
		    do j=consisD_yi,consisD_yf,1
			dy(i,j,ka-zi+1)=dy(i,j,ka-zi+1)-Ku*hx_inc(j)
		    enddo
		enddo
	    endif

	    if(zi.le.(kb+1).and.zf.ge.(kb+1))then
		do i=consisD_xi,consisD_xf,1
		    do j=consisD_yi,consisD_yf,1
			dy(i,j,kb-zi+2)=dy(i,j,kb-zi+2)+Ku*hx_inc(j)
		    enddo
		enddo
	    endif
	end subroutine

!	Dz
	subroutine consistCond_dz(absrangesx,absrangesy,absrangesz)
	type(Irange), intent(in) :: absrangesx,absrangesy,absrangesz
	integer :: yi,yf,i,k
!	xi=absrangesx%init
!	xf=absrangesx%ends
	yi=absrangesy%init
	yf=absrangesy%ends
!	zi=absrangesz%init
!	zf=absrangesz%ends
	    if(yi.le.ja.and.yf.ge.ja)then
		do i=consisD_xi,consisD_xf,1
		    do k=consisD_zi,consisD_zf,1
			dz(i,ja-yi+1,k)=dz(i,ja-yi+1,k)+Ku*hx_inc(ja-1)
		    enddo
		enddo
	    endif
	    if(yi.le.jb.and.yf.ge.jb)then
		do i=consisD_xi,consisD_xf,1
		    do k=consisD_zi,consisD_zf,1
			dz(i,jb-yi+1,k)=dz(i,jb-yi+1,k)-Ku*hx_inc(jb)
		    enddo
		enddo
	    endif
	end subroutine
!c*********************************************************
!c	Consistency conditions for H
!c*********************************************************
!	Hx
	subroutine consistCond_hx(absrangesx,absrangesy,absrangesz)
	type(Irange), intent(in) :: absrangesx,absrangesy,absrangesz
	integer :: yi,yf,i,k
!	xi=absrangesx%init
!	xf=absrangesx%ends
	yi=absrangesy%init
	yf=absrangesy%ends
!	zi=absrangesz%init
!	zf=absrangesz%ends
	    if(yi.le.(ja-1).and.yf.ge.(ja-1))then
		do i=consisH_xi,consisH_xf,1
		    do k=consisH_zi,consisH_zf,1
			hx(i,ja-yi,k)=hx(i,ja-yi,k)+Ku*ez_inc(ja)
		    enddo
		enddo
	    endif

	    if(yi.le.jb.and.yf.ge.jb)then
		do i=consisH_xi,consisH_xf,1
		    do k=consisH_zi,consisH_zf,1
			hx(i,jb-yi+1,k)=hx(i,jb-yi+1,k)-Ku*ez_inc(jb)
		    enddo
		enddo
	    endif
	end subroutine
!	Hy
	subroutine consistCond_hy(absrangesx,absrangesy,absrangesz)
	type(Irange), intent(in) :: absrangesx,absrangesy,absrangesz
	integer :: xi,xf,j,k
	xi=absrangesx%init
	xf=absrangesx%ends
!	yi=absrangesy%init
!	yf=absrangesy%ends
!	zi=absrangesz%init
!	zf=absrangesz%ends
	    if(xi.le.(ia-1).and.xf.ge.(ia-1))then
		do j=consisH_yi,consisH_yf,1
		    do k=consisH_zi,consisH_zf,1
			hy(ia-xi,j,k)=hy(ia-xi,j,k)-Ku*ez_inc(j)
		    enddo
		enddo
	    endif

	    if(xi.le.ib.and.xf.ge.ib)then
		do j=consisH_yi,consisH_yf,1
		    do k=consish_zi,consisH_zf,1
			hy(ib-xi+1,j,k)=hy(ib-xi+1,j,k)+Ku*ez_inc(j)
		    enddo
		enddo
	    endif
	end subroutine
!c*********************************************
!c	synchronization for the Hz_inc field
!c*********************************************
!	subroutine synchonize_hzinc()
!	    if(ycord.gt.0)then
!!c	    if(n.eq.1)then
!!!c		print *,"send(H1):",myrank
!!c	    end if
!		call mpi_send(hz_inc(1),1,mpi_double_precision,neighborN,&
!		    tagCH1_NS,mpi_comm_world,ierror)
!	    endif
!	    if(ycord.lt.(Ymax-1))then
!!c	    if(n.eq.1)then
!!c		print *,"recv(H1):",myrank
!!c	    end if
!		call mpi_recv(hz_inc(n_y+1),1,mpi_double_precision,&
!		    neighborS,tagCH1_NS,mpi_comm_world,status,ierror)
!	    endif
!	    if(ycord.lt.(Ymax-1)) then
!!c	    if(n.eq.1)then
!!c		print *,"send(H2):",myrank
!!c	    end if
!	call mpi_send(hz_inc(n_y),1,mpi_double_precision,neighborS,&
!		tagCH2_NS,mpi_comm_world,ierror)
!	    endif
!	    if(ycord.gt.0) then
!!c	    if(n.eq.1)then
!!c		print *,"recv(H2):",myrank
!!c	    end if
!	call mpi_recv(hz_inc(0),1,mpi_double_precision,&
!		neighborN,tagCH2_NS,mpi_comm_world,status,ierror)
!	    endif
!	end subroutine
!
!!c*********************************************
!!c	synchronization for the ex_inc field
!!c*********************************************
!	subroutine synchonize_exinc()
!	    if(ycord.gt.0)then
!!c	    if(n.eq.1)then
!!c		print *,"send(Ex1):",myrank
!!c	    end if
!	call mpi_send(ex_inc(1),1,mpi_double_precision,neighborN,&
!		tagCE1_NS,mpi_comm_world,ierror)
!	    endif
!	    if(ycord.lt.(Ymax-1))then
!!c	    if(n.eq.1)then
!!c		print *,"recv(Ex1):",myrank
!!c	    end if
!	call mpi_recv(ex_inc(n_y+1),1,mpi_double_precision,&
!		neighborS,tagCE1_NS,mpi_comm_world,status,ierror)
!	    end if
!	    if(ycord.lt.(Ymax-1)) then
!!c	    if(n.eq.1)then
!!c		print *,"send(Ex2):",myrank
!!c	    end if
!	call mpi_send(ex_inc(n_y),1,mpi_double_precision,neighborS,&
!		tagCE2_NS,mpi_comm_world,ierror)
!	endif
!	    if(ycord.gt.0) then
!!c	    if(n.eq.1)then
!!c		print *,"recv(Ex2):",myrank
!!c	    end if
!	call mpi_recv(ex_inc(0),1,mpi_double_precision,&
!		neighborN,tagCE2_NS,mpi_comm_world,status,ierror)
!	    endif
!!!	end subroutine