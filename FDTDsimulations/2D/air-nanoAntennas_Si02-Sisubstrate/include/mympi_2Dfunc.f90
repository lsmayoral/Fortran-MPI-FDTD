!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c		determine the neightboring ranks
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_neighbor(coord,dir,Xmax,Ymax,comm,neighbor)
	implicit none
	integer,dimension(2), intent(in) :: coord
	integer, intent(in) :: dir,comm,Xmax,Ymax
	integer, dimension(2) :: newcoord
	integer, parameter :: NORTH=1
	integer, parameter :: SOUTH=2
	integer, parameter :: EAST=3
	integer, parameter :: WEST=4
	integer :: x,y,error,neighbor
	if(dir.eq.NORTH)then
	    if(coord(1).eq.(Xmax-1).and.PBC_x)then
		x=0
		y=coord(2)
	    else
		x=coord(1)+1
		y=coord(2)
	    end if
	else if(dir.eq.SOUTH)then
	    if(coord(1).eq.0.and.PBC_x)then
		x=Xmax-1
		y=coord(2)
	    else
		x=coord(1)-1
		y=coord(2)
	    end if
	else if(dir.eq.EAST)then
	    if(coord(2).eq.(Ymax-1).and.PBC_y)then
		x=coord(1)
		y=0
	    else
		x=coord(1)
		y=coord(2)+1
	    end if
	else if(dir.eq.WEST)then
	    if(coord(2).eq.0.and.PBC_y)then
		x=coord(1)
		y=Ymax-1
	    else
		x=coord(1)
		y=coord(2)-1
	    end if
	else
	    x=-1
	    y=-1
	end if
	newcoord(1)=x
	newcoord(2)=y
	if(x.lt.0.or.y.lt.0.or.x.gt.(Xmax-1).or.y.gt.&
		(Ymax-1))then
	neighbor=-1
	else
	    CALL MPI_Cart_rank(comm,newcoord,neighbor,error)
	end if
	end subroutine get_neighbor
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
!c		create Matrix blocks creations
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
	subroutine create_mpi_planes()
!c	for y line
	    CALL MPI_TYPE_VECTOR(n_y+2,1,n_x+2,MPI_double_precision,&
		liney,ierror)
	    CALL MPI_TYPE_COMMIT(liney,ierror)
!c	for x line is no necesary (contiguous data)
	end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c	Initialize MPI environment
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine my_mpi_init()
	    call mpi_init(ierror)
	    call mpi_comm_size(mpi_comm_world,nproces,ierror)
	    call mpi_comm_rank(mpi_comm_world,myrank,ierror)

!c		Define parallelization grid layout
	    call mpi_cart_create(mpi_comm_world,2,Ndim,periodic,&
		.false.,COMM_2D,ierror)
	    call mpi_cart_coords(COMM_2D,myrank,2,coords,ierror)
	    xcord=coords(1)
	    ycord=coords(2)
!c	    write(*,"(a7,3I6)")"myrank:",myrank,xcord,ycord
	    call mpi_barrier(mpi_comm_world,ierror)
	    call get_neighbor(coords,NORTH,Xmax,Ymax,comm_2D,neighborN)
	    call get_neighbor(coords,SOUTH,Xmax,Ymax,comm_2D,neighborS)
	    call get_neighbor(coords,EAST,Xmax,Ymax,comm_2D,neighborE)
	    call get_neighbor(coords,WEST,Xmax,Ymax,comm_2D,neighborW)
!c	    write(*,"(a10,5I6)")"neighbors:",myrank,neighborN,neighborS,
!c     *		neighborE,neighborW
	    tagE1_WE=1
	    tagE2_WE=2
	    tagE1_NS=3
	    tagE2_NS=4
	    tagH1_WE=5
	    tagH2_WE=6
	    tagH1_NS=7
	    tagH2_NS=8
	    call set_index()

	end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c		set rank indexs
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine set_index()
	    absIxrange%init=xcord*n_x+1
	    absIxrange%ends=absIxrange%init+n_x-1
	    absIyrange%init=ycord*N_y+1
	    absIyrange%ends=absIyrange%init+N_y-1

	    locIxrange%init=1
	    locIxrange%ends=N_x
	    locIyrange%init=1
	    locIyrange%ends=N_y

	    if(.not.PBC_x)then
		if(xcord.eq.0)then
		    locIxrange%init=1
		end if
		if(xcord.eq.(Xmax-1))then
		    locIxrange%ends=N_x
		end if
	    end if
	    if(.not.PBC_y)then
		if(ycord.eq.0)then
		    locIyrange%init=npmly
		end if
		if(ycord.eq.(Ymax-1))then
		    locIyrange%ends=N_y-npmly
		end if
	    end if

	write(*,"(a10,I4,7I6)")"absRange:",myrank,absIxrange%init,&
		absIxrange%ends,absIyrange%init,absIyrange%ends
	write(*,"(a10,I4,4I6)")"locRange:",myrank,locIxrange%init,&
	    locIxrange%ends,locIyrange%init,locIyrange%ends
	end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c		E field synchonization calls
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine interchange_E()
!c	First step: Ex(1,N_y+1)<-Ex(1,1) and Hz(1,N_y)->Hz(1,0)
	    if(PBC_y)then
	call mpi_sendrecv(ex(0,1),N_x+2,mpi_double_precision,neighborW,&
	tagE1_WE,ex(0,N_y+1),N_x+2,mpi_double_precision,&
		neighborE,tagE1_WE,mpi_comm_world,status,ierror)

	call mpi_sendrecv(ex(0,N_y),N_x+2,mpi_double_precision,neighborE,&
		tagE2_WE,ex(0,0),N_x+2,mpi_double_precision,&
		neighborW,tagE2_WE,mpi_comm_world,status,ierror)
	else
	    if(ycord.gt.0) then
	call mpi_send(ex(0,1),N_x+2,mpi_double_precision,neighborW,&
		tagE1_WE,mpi_comm_world,ierror)
	    endif
	    if(ycord.lt.(Ymax-1)) then
	call mpi_recv(ex(0,N_y+1),N_x+2,mpi_double_precision,&
		neighborE,tagE1_WE,mpi_comm_world,status,ierror)
	    endif
	    if(ycord.lt.(Ymax-1)) then
	call mpi_send(ex(0,N_y),N_x+2,mpi_double_precision,neighborE,&
		tagE2_WE,mpi_comm_world,ierror)
	    endif
	    if(ycord.gt.0) then
	call mpi_recv(ex(0,0),N_x+2,mpi_double_precision,&
		neighborW,tagE2_WE,mpi_comm_world,status,ierror)
	    endif
	end if

!c	second step: Ey(N_x+1,1)<-Ey(1,1) and Ey(N_x,1)->Ey(0,1)
	if(PBC_x)then
!	    call mpi_sendrecv(ey(1,npmly+1),1,liney,neighborS,&
!		tagE1_NS,ey(N_x+1,npmly+1),1,liney,&
!		neighborN,tagE1_NS,mpi_comm_world,status,ierror)
!	    call mpi_sendrecv(ey(N_x,npmly+1),1,liney,neighborN,&
!		tagE2_NS,ey(0,npmly+1),1,liney,&
!		neighborS,tagE2_NS,mpi_comm_world,status,ierror)
	    call mpi_sendrecv(ey(1,0),1,liney,neighborS,&
		tagE1_NS,ey(N_x+1,0),1,liney,&
		neighborN,tagE1_NS,mpi_comm_world,status,ierror)

	    call mpi_sendrecv(ey(N_x,0),1,liney,neighborN,&
		tagE2_NS,ey(0,0),1,liney,&
		neighborS,tagE2_NS,mpi_comm_world,status,ierror)
	else
	    if(xcord.gt.0)then
		call mpi_send(ey(1,0),1,liney,neighborS,&
		tagE1_NS,mpi_comm_world,ierror)
	    end if
	    if(xcord.lt.(Xmax-1))then
		call mpi_recv(ey(N_x+1,0),1,liney,&
		neighborN,tagE1_NS,mpi_comm_world,status,ierror)
	    end if
	    if(xcord.lt.(Xmax-1))then
		call mpi_send(ey(N_x,0),1,liney,neighborN,&
		tagE2_NS,mpi_comm_world,ierror)
	    end if
	    if(xcord.gt.0)then
		call mpi_recv(ey(0,0),1,liney,&
		neighborS,tagE2_NS,mpi_comm_world,status,ierror)
	    end if
	end if
	call mpi_barrier(mpi_comm_world,ierror)
	end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c		H field synchonization calls
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine interchange_H()
!c	First step: Hz(1,N_y+1)<-Hz(1,1) and Hz(1,N_y)->Hz(1,0)
	if(PBC_y)then
	call mpi_sendrecv(hz(0,1),N_x+2,mpi_double_precision,neighborW,&
		tagH1_WE,hz(0,N_y+1),N_x+2,mpi_double_precision,&
		neighborE,tagH1_WE,mpi_comm_world,status,ierror)

	call mpi_sendrecv(hz(0,N_y),N_x+2,mpi_double_precision,neighborE,&
		tagH2_WE,hz(0,0),N_x+2,mpi_double_precision,&
		neighborW,tagH2_WE,mpi_comm_world,status,ierror)
	else
	    if(ycord.gt.0) then
	call mpi_send(hz(0,1),N_x+2,mpi_double_precision,neighborW,&
		tagH1_WE,mpi_comm_world,ierror)
	    endif
	    if(ycord.lt.(Ymax-1)) then
	call mpi_recv(hz(0,N_y+1),N_x+2,mpi_double_precision,&
		neighborE,tagH1_WE,mpi_comm_world,status,ierror)
	    endif
	    if(ycord.lt.(Ymax-1)) then
	call mpi_send(hz(0,N_y),N_x+2,mpi_double_precision,neighborE,&
		tagH2_WE,mpi_comm_world,ierror)
	    endif
	    if(ycord.gt.0) then
	call mpi_recv(hz(0,0),N_x+2,mpi_double_precision,&
		neighborW,tagH2_WE,mpi_comm_world,status,ierror)
	    endif
	end if
!c	second step: Hz(N_x+1,1)<-Hz(1,1) and Hz(N_x,1)->Hz(0,1)
	if(PBC_x)then
!	    call mpi_sendrecv(hz(1,npmly+1),1,liney,neighborS,&
!		tagH1_NS,hz(N_x+1,npmly+1),1,liney,&
!		neighborN,tagH1_NS,mpi_comm_world,status,ierror)
!	    call mpi_sendrecv(hz(N_x,npmly+1),1,liney,neighborN,&
!		tagH2_NS,hz(0,npmly+1),1,liney,&
!		neighborS,tagH2_NS,mpi_comm_world,status,ierror)
	    call mpi_sendrecv(hz(1,0),1,liney,neighborS,&
		tagH1_NS,hz(N_x+1,0),1,liney,&
		neighborN,tagH1_NS,mpi_comm_world,status,ierror)

	    call mpi_sendrecv(hz(N_x,0),1,liney,neighborN,&
		tagH2_NS,hz(0,0),1,liney,&
		neighborS,tagH2_NS,mpi_comm_world,status,ierror)
	else
	    if(xcord.gt.0)then
		call mpi_send(hz(1,0),1,liney,neighborS,&
		tagH1_NS,mpi_comm_world,ierror)
	    end if
	    if(xcord.lt.(Xmax-1))then
		call mpi_recv(hz(N_x+1,0),1,liney,&
		neighborN,tagH1_NS,mpi_comm_world,status,ierror)
	    end if

	    if(xcord.lt.(Xmax-1))then
		call mpi_send(hz(N_x,0),1,liney,neighborN,&
		tagH2_NS,mpi_comm_world,ierror)
	    end if
	    if(xcord.gt.0)then
		call mpi_recv(hz(0,0),1,liney,&
		neighborS,tagH2_NS,mpi_comm_world,status,ierror)
	    end if
	end if

	call mpi_barrier(mpi_comm_world,ierror)
	end subroutine

