!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		determine the neightboring ranks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine get_neigbor(coord,dir,Xmax,Ymax,Zmax,comm,neigbor)
	implicit none
	integer,dimension(3), intent(in) :: coord
	integer, intent(in) :: dir,comm,Xmax,Ymax,Zmax
	integer, dimension(3) :: newcoord
	integer, parameter :: NORTH=1
	integer, parameter :: SOUTH=2
	integer, parameter :: EAST=3
	integer, parameter :: WEST=4
	integer, parameter :: UP=5
	integer, parameter :: DOWN=6
	integer :: x,y,z,error,me,neigbor
	if(dir.eq.NORTH)then
	    x=coord(1)
	    z=coord(3)

	    if(coord(2).eq.(0).and.PBC_y)then
		y=Ymax-1
	    else
		y=coord(2)-1
	    endif

	else if(dir.eq.SOUTH)then
	    x=coord(1)
	    z=coord(3)

	    if(coord(2).eq.(Ymax-1).and.PBC_y)then
		y=0
	    else
		y=coord(2)+1
	    endif

	else if(dir.eq.EAST)then
	    x=coord(1)
	    y=coord(2)
	    if(coord(3).eq.(Zmax-1).and.PBC_z)then
		z=0
	    else
		z=coord(3)+1
	    endif
	else if(dir.eq.WEST)then
	    x=coord(1)
	    y=coord(2)
	    if(coord(3).eq.(0).and.PBC_z)then
		z=Zmax-1
	    else
		z=coord(3)-1
	    endif
	else if(dir.eq.UP)then
	    y=coord(2)
	    z=coord(3)
	    if(coord(1).eq.(Xmax-1).and.PBC_x)then
		x=0
	    else
		x=coord(1)+1
	    endif
	else if(dir.eq.DOWN)then
	    y=coord(2)
	    z=coord(3)
	    if(coord(1).eq.(0).and.PBC_x)then
		x=Xmax-1
	    else
		x=coord(1)-1
	    endif
	else
	    x=-1
	    y=-1
	    z=-1
	endif
	newcoord(1)=x
	newcoord(2)=y
	newcoord(3)=z
	if(x.lt.0.or.y.lt.0.or.z.lt.0.or.x.gt.(Xmax-1).or.y.gt.(Ymax-1)&
	.or.z.gt.(Zmax-1))then
	neigbor=-1
	else
	    CALL MPI_Cart_rank(comm,newcoord,neigbor,error)
	endif
	end subroutine get_neigbor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Initialize MPI environment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine my_mpi_init()
	call mpi_init(ierror)
	call mpi_comm_size(mpi_comm_world,nproces,ierror)
	call mpi_comm_rank(mpi_comm_world,myrank,ierror)
!		Define parallelization grid layout
	call mpi_cart_create(mpi_comm_world,3,Ndim,periodic,&
		.false.,COMM_3D,ierror)
	call mpi_cart_coords(COMM_3D,myrank,3,coords,ierror)
	xcord=coords(1)
	ycord=coords(2)
	zcord=coords(3)
	write(*,"(a7,I4,3I2)")"myrank=",myrank,xcord,ycord,zcord
	call mpi_barrier(mpi_comm_world,ierror)

	call get_neigbor(coords,NORTH,Xmax,Ymax,Zmax,comm_3D,neigborN)
	call get_neigbor(coords,SOUTH,Xmax,Ymax,Zmax,comm_3D,neigborS)
	call get_neigbor(coords,EAST,Xmax,Ymax,Zmax,comm_3D,neigborE)
	call get_neigbor(coords,WEST,Xmax,Ymax,Zmax,comm_3D,neigborW)
	call get_neigbor(coords,UP,Xmax,Ymax,Zmax,comm_3D,neigborU)
	call get_neigbor(coords,DOWN,Xmax,Ymax,Zmax,comm_3D,neigborD)
	tagE1_WE=1
	tagE2_WE=2
	tagE3_WE=3
	tagE4_WE=4
	tagE1_NS=5
	tagE2_NS=6
	tagE3_NS=7
	tagE4_NS=8
	tagH1_WE=9
	tagH2_WE=10
	tagH3_WE=11
	tagH4_WE=12
	tagH1_NS=13
	tagH2_NS=14
	tagH3_NS=15
	tagH4_NS=16
	tagE1_UD=17
	tagE2_UD=18
	tagE3_UD=19
	tagE4_UD=20
	tagH1_UD=21
	tagH2_UD=22
	tagH3_UD=23
	tagH4_UD=24

	call set_index()

	END subroutine my_mpi_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Set local and absolute indixes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine set_index()

        absIxrange%init=xcord*n_x+1
        absIxrange%ends=absIxrange%init+n_x-1
        absIyrange%init=ycord*N_y+1
        absIyrange%ends=absIyrange%init+N_y-1
        absIzrange%init=zcord*N_z+1
        absIzrange%ends=absIzrange%init+N_z-1

        locIxrange%init=1
        locIxrange%ends=N_x
        locIyrange%init=1
        locIyrange%ends=N_y
        locIzrange%init=1
        locIzrange%ends=N_z

        if(.not.PBC_x)then
	    if(xcord.eq.0)then
		locIxrange%init=1!npmlx+1
	    endif
	    if(xcord.eq.(Xmax-1))then
		locIxrange%ends=N_x!-npmlx
	    endif
        endif

        if(.not.PBC_y)then
	    if(ycord.eq.0)then
		locIyrange%init=npmly+1
	    endif
	    if(ycord.eq.(Ymax-1))then
		locIyrange%ends=N_y-npmly
	    endif
        endif

        if(.not.PBC_z)then
	    if(zcord.eq.0)then
		locIzrange%init=1!npmlz+1
	    endif
	    if(zcord.eq.(Zmax-1))then
		locIzrange%ends=N_z!-npmlz
	    endif
        endif

        if(myrank.le.9)then
	    write(*,"(a10,I2,6I6)")"locRange:",myrank,&
	    locIxrange%init,locIxrange%ends,&
	    locIyrange%init,locIyrange%ends,&
	    locIzrange%init,locIzrange%ends
	    write(*,"(a10,I2,6I9)")"neighbors:",myrank,&
	    neigborN,neigborS,neigborE,neigborW,neigborU,neigborD
	    write(*,"(a10,I2,6I6)")"absRange:",myrank,&
		absIxrange%init,absIxrange%ends,&
		absIyrange%init,absIyrange%ends,&
		absIzrange%init,absIzrange%ends
        else if(myrank.ge.10.and.myrank.le.99)then
	    write(*,"(a10,I4,6I6)")"locRange:",myrank,&
	    locIxrange%init,locIxrange%ends,&
	    locIyrange%init,locIyrange%ends,&
	    locIzrange%init,locIzrange%ends
	    write(*,"(a10,I3,6I9)")"neighbors:",myrank,&
	    neigborN,neigborS,neigborE,neigborW,neigborU,neigborD
	    write(*,"(a10,I3,6I6)")"absRange:",myrank,&
	    	absIxrange%init,absIxrange%ends,&
		absIyrange%init,absIyrange%ends,&
		absIzrange%init,absIzrange%ends
        else
	    write(*,"(a10,I4,6I6)")"locRange:",myrank,&
	    locIxrange%init,locIxrange%ends,&
	    locIyrange%init,locIyrange%ends,&
	    locIzrange%init,locIzrange%ends
	    write(*,"(a10,I4,6I9)")"neighbors:",myrank,&
	    neigborN,neigborS,neigborE,neigborW,neigborU,neigborD
	    write(*,"(a10,I4,6I6)")"absRange:",myrank,&
		absIxrange%init,absIxrange%ends,&
		absIyrange%init,absIyrange%ends,&
		absIzrange%init,absIzrange%ends
        endif
	end subroutine set_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		create Matrix blocks creations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine create_mpi_planes()
!	for XZ plane
!	    CALL MPI_TYPE_VECTOR(N_z+2,N_x+2,(N_x+2)*(N_y+2),&
!		MPI_double_precision,TypeMatrXZ,ierror)

	    CALL MPI_TYPE_VECTOR(N_z+2,N_x+2,(N_x+2)*(N_y+2),&
		MPI_double_precision,TypeMatrXZ,ierror)
	    CALL MPI_TYPE_COMMIT(TypeMatrXZ,ierror)

!	for YZ plane
	    CALL MPI_TYPE_VECTOR((N_y+2)*(N_z+2),1,N_x+2,&
		MPI_double_precision,TypeMatrYZ,ierror)
	    CALL MPI_TYPE_COMMIT(TypeMatrYZ,ierror)
!	for XY plane
	    countXY=(N_x+2)*(N_y+2)
	end subroutine create_mpi_planes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		E field synchonization calls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c		E field synchonization calls
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine interchange_E()
	if(n.eq.1)then
	    print *,"starting E interchange",myrank
	endif
!		firs step of syncronization [send(F(1)), recv(F(N_i+1))
	    if(PBC_x)then
		if(n.eq.1)then
		    print *,"PBCX_Ey1 interchange",myrank
		endif
		call MPI_SENDRECV(Ey(1,0,0),1,TypeMatrYZ,neigborD,tagE1_UD,&
		    Ey(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagE1_UD,comm_3D,status,&
		    ierror)

		if(n.eq.1)then
		    print *,"PBCX_Ez1 interchange",myrank
		endif

		call MPI_SENDRECV(Ez(1,0,0),1,TypeMatrYZ,neigborD,tagE2_UD,&
		    Ez(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagE2_UD,comm_3D,status,&
		    ierror)
	    else
		if(xcord.gt.0)then
		if(n.eq.1)then
		    print *,"SendEy_D1",myrank
		endif
		    call MPI_Send(Ey(1,0,0),1,TypeMatrYZ,neigborD,tagE1_UD,&
			comm_3D,status,ierror)
		endif
		if(xcord.lt.(Xmax-1))then
		if(n.eq.1)then
		    print *,"RecvEy_D1",myrank
		endif
		    call MPI_Recv(Ey(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagE1_UD,comm_3D,&
			status,ierror)
		endif
		if(xcord.gt.0)then
		if(n.eq.1)then
		    print *,"SendEz_D1",myrank
		endif
		    call MPI_Send(Ez(1,0,0),1,TypeMatrYZ,neigborD,tagE2_UD,&
			comm_3D,status,ierror)
		endif
		if(xcord.lt.(Xmax-1))then
		if(n.eq.1)then
		    print *,"RecvEy_D1",myrank
		endif
		    call MPI_Recv(Ez(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagE2_UD,comm_3D,&
			status,ierror)
		endif
	    endif

	    if(PBC_y)then
		if(n.eq.1)then
		    print *,"PBCY_Ex1 interchange",myrank
		endif
		call MPI_SENDRECV(Ex(0,1,0),1,TypeMatrXZ,neigborN,tagE1_NS,&
		    Ex(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagE1_NS,comm_3D,status,&
		    ierror)
		if(n.eq.1)then
		    print *,"PBCY_Ez1 interchange",myrank
		endif
		call MPI_SENDRECV(Ez(0,1,0),1,TypeMatrXZ,neigborN,tagE2_NS,&
		    Ez(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagE2_NS,comm_3D,status,&
		    ierror)
	    else
		if(ycord.gt.0)then
		if(n.eq.1)then
		    print *,"SendEx_S1",myrank
		endif
		    call MPI_Send(Ex(0,1,0),1,TypeMatrXZ,neigborN,tagE1_NS,&
 			comm_3D,status,ierror)
		endif
		if(ycord.lt.(Ymax-1))then
		if(n.eq.1)then
		    print *,"RecvEx_S1",myrank
		endif
		    call MPI_Recv(Ex(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagE1_NS,comm_3D,&
			status,ierror)
		endif
		if(ycord.gt.0)then
		if(n.eq.1)then
		    print *,"SendEz_S1",myrank
		endif
		    call MPI_Send(Ez(0,1,0),1,TypeMatrXZ,neigborN,tagE2_NS,&
			comm_3D,status,ierror)
		endif
		if(ycord.lt.(Ymax-1))then
		if(n.eq.1)then
		    print *,"RecvEz_S1",myrank
		endif
		    call MPI_Recv(Ez(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagE2_NS,comm_3D,&
			status,ierror)
		endif
	    endif

	    if(PBC_z)then
		if(n.eq.1)then
		    print *,"PBCZ_Ex1 interchange",myrank
		endif
		call MPI_SENDRECV(Ex(0,0,1),countXY,MPI_double_precision,&
		    neigborW,tagE1_WE,Ex(0,0,N_z+1),countXY,MPI_double_precision,&
		    neigborE,tagE1_WE,comm_3D,status,ierror)

		if(n.eq.1)then
		    print *,"PBCY_Ey1 interchange",myrank
		endif

		call MPI_SENDRECV(Ey(0,0,1),countXY,MPI_double_precision,&
		    neigborW,tagE2_WE,Ey(0,0,N_z+1),countXY,MPI_double_precision,&
		    neigborE,tagE2_WE,comm_3D,status,ierror)
	    else
		if(zcord.gt.0)then
		if(n.eq.1)then
		    print *,"SendEx_E1",myrank
		endif
		    call MPI_Send(Ex(0,0,1),countXY,MPI_double_precision,&
			neigborW,tagE1_WE,comm_3D,status,ierror)
		endif
		if(zcord.lt.Zmax-1)then
		if(n.eq.1)then
		    print *,"RecvEx_E1",myrank
		endif
		    call MPI_Recv(Ex(0,0,N_z+1),countXY,MPI_double_precision,&
			neigborE,tagE1_WE,comm_3D,status,ierror)
		endif

		if(zcord.gt.0)then
		if(n.eq.1)then
		    print *,"SendEy_E1",myrank
		endif
		    call MPI_Send(Ey(0,0,1),countXY,MPI_double_precision,&
			neigborW,tagE2_WE,comm_3D,status,ierror)
		endif
		if(zcord.lt.Zmax-1)then
		if(n.eq.1)then
		    print *,"RecvEy_E1",myrank
		endif
		    call MPI_Recv(Ey(0,0,N_z+1),countXY,MPI_double_precision,neigborE,&
			tagE2_WE,comm_3D,status,ierror)
		endif
	    endif
!	second step of syncronization [send(F(N_i)), recv(F(0))
	    if(PBC_x)then
		if(n.eq.1)then
		    print *,"PBCX_Ey2 interchange",myrank
		endif
		call MPI_SENDRECV(Ey(N_x,0,0),1,TypeMatrYZ,neigborU,tagE3_UD,&
		    Ey(0,0,0),1,TypeMatrYZ,neigborD,tagE3_UD,comm_3D,status,&
		    ierror)
		if(n.eq.1)then
		    print *,"PBCX_Ez2 interchange",myrank
		endif
		call MPI_SENDRECV(Ez(N_x,0,0),1,TypeMatrYZ,neigborU,tagE4_UD,&
		    Ez(0,0,0),1,TypeMatrYZ,neigborD,tagE4_UD,comm_3D,status,&
		    ierror)
	    else
		if(xcord.gt.0)then
		if(n.eq.1)then
		    print *,"RecvEy_D2",myrank
		endif
		    call MPI_Recv(Ey(0,0,0),1,TypeMatrYZ,neigborD,tagE3_UD,&
			comm_3D,status,ierror)
		endif
		if(xcord.lt.(Xmax-1))then
		if(n.eq.1)then
		    print *,"SendvEy_D2",myrank
		endif
		    call MPI_Send(Ey(N_x,0,0),1,TypeMatrYZ,neigborU,tagE3_UD,comm_3D,&
			status,ierror)
		endif
		if(xcord.gt.0)then
		if(n.eq.1)then
		    print *,"RecvEz_D2",myrank
		endif
		    call MPI_Recv(Ez(0,0,0),1,TypeMatrYZ,neigborD,tagE4_UD,&
			comm_3D,status,ierror)
		endif
		if(xcord.lt.(Xmax-1))then
		if(n.eq.1)then
		    print *,"SendEz_D2",myrank
		endif
		    call MPI_Send(Ez(N_x,0,0),1,TypeMatrYZ,neigborU,tagE4_UD,comm_3D,&
			status,ierror)
		endif
	    endif
	    if(PBC_y)then
		if(n.eq.1)then
		    print *,"PBCY_Ex2 interchange",myrank
		endif
		call MPI_SENDRECV(Ex(0,N_y,0),1,TypeMatrXZ,neigborS,tagE3_NS,&
		    Ex(0,0,0),1,TypeMatrXZ,neigborN,tagE3_NS,comm_3D,status,&
		    ierror)
		if(n.eq.1)then
		    print *,"PBCY_Ez2 interchange",myrank
		endif
		call MPI_SENDRECV(Ez(0,N_y,0),1,TypeMatrXZ,neigborS,tagE4_NS,&
		    Ez(0,0,0),1,TypeMatrXZ,neigborN,tagE4_NS,comm_3D,status,&
		    ierror)
	    else
		if(ycord.gt.0)then
		if(n.eq.1)then
		    print *,"RecvEx_N2",myrank
		endif
		    call MPI_Recv(Ex(0,0,0),1,TypeMatrXZ,neigborN,tagE3_NS,&
 			comm_3D,status,ierror)
		endif
		if(ycord.lt.(Ymax-1))then
		if(n.eq.1)then
		    print *,"SendEx_N2",myrank
		endif
		    call MPI_Send(Ex(0,N_y,0),1,TypeMatrXZ,neigborS,tagE3_NS,comm_3D,&
			status,ierror)
		endif
		if(ycord.gt.0)then
		if(n.eq.1)then
		    print *,"RecvEz_N2",myrank
		endif
		    call MPI_Recv(Ez(0,0,0),1,TypeMatrXZ,neigborN,tagE4_NS,&
			comm_3D,status,ierror)
		endif
		if(ycord.lt.(Ymax-1))then
		if(n.eq.1)then
		    print *,"RecvEz_N2",myrank
		endif
		    call MPI_Send(Ez(0,N_y,0),1,TypeMatrXZ,neigborS,tagE4_NS,comm_3D,&
			status,ierror)
		endif
	    endif

	    if(PBC_z)then
		if(n.eq.1)then
		    print *,"PBCZ_Ex2 interchange",myrank
		endif
		call MPI_SENDRECV(Ex(0,0,N_z),countXY,MPI_double_precision,&
		    neigborE,tagE3_WE,Ex(0,0,0),countXY,MPI_double_precision,&
		    neigborW,tagE3_WE,comm_3D,status,ierror)
		if(n.eq.1)then
		    print *,"PBCZ_Ey2 interchange",myrank
		endif
		call MPI_SENDRECV(Ey(0,0,N_z),countXY,MPI_double_precision,&
		    neigborE,tagE4_WE,Ey(0,0,0),countXY,MPI_double_precision,&
		    neigborW,tagE4_WE,comm_3D,status,ierror)
	    else
		if(zcord.gt.0)then
		if(n.eq.1)then
		    print *,"RecvEx_E2",myrank
		endif
		    call MPI_Recv(Ex(0,0,0),countXY,MPI_double_precision,&
			neigborW,tagE3_WE,comm_3D,status,ierror)
		endif
		if(zcord.lt.(Zmax-1))then
		if(n.eq.1)then
		    print *,"SendEx_E2",myrank
		endif
		    call MPI_Send(Ex(0,0,N_z),countXY,MPI_double_precision,&
			neigborE,tagE3_WE,comm_3D,status,ierror)
		endif

		if(zcord.gt.0)then
		    if(n.eq.1)then
			print *,"RecvEy_E2",myrank
		    endif
		    call MPI_Recv(Ey(0,0,0),countXY,MPI_double_precision,&
			neigborW,tagE4_WE,comm_3D,status,ierror)
		endif
		if(zcord.lt.(Zmax-1))then
		    if(n.eq.1)then
			print *,"SendEy_E2",myrank
		    endif
		    call MPI_Send(Ey(0,0,N_z),countXY,MPI_double_precision,neigborE,&
			tagE4_WE,comm_3D,status,ierror)
		endif
	    endif
	if(n.eq.1)then
	    print *,"ending E interchange",myrank
	endif
	call mpi_barrier(mpi_comm_world,ierror)
	end subroutine interchange_E
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		H field synchonization calls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine interchange_H()
	if(n.eq.1)then
	    print *,"starting H interchange",myrank
	endif
!	first step of syncronization [send(F(N_i)), recv(F(0))
	    if(PBC_x)then
		call MPI_SENDRECV(Hy(N_x,0,0),1,TypeMatrYZ,neigborU,tagH3_UD,&
		    Hy(0,0,0),1,TypeMatrYZ,neigborD,tagH3_UD,comm_3D,status,&
		    ierror)

		call MPI_SENDRECV(Hz(N_x,0,0),1,TypeMatrYZ,neigborU,tagH4_UD,&
		    Hz(0,0,0),1,TypeMatrYZ,neigborD,tagH4_UD,comm_3D,status,&
		    ierror)
	    else
		if(xcord.gt.0)then
		    call MPI_Recv(Hy(0,0,0),1,TypeMatrYZ,neigborD,tagH3_UD,&
			comm_3D,status,ierror)
		endif
		if(xcord.lt.(Xmax-1))then
		    call MPI_Send(Hy(N_x,0,0),1,TypeMatrYZ,neigborU,tagH3_UD,comm_3D,&
			status,ierror)
		endif
		if(xcord.gt.0)then
		    call MPI_Recv(Hz(0,0,0),1,TypeMatrYZ,neigborD,tagH4_UD,&
			comm_3D,status,ierror)
		endif
		if(xcord.lt.(Xmax-1))then
		    call MPI_Send(Hz(N_x,0,0),1,TypeMatrYZ,neigborU,tagH4_UD,comm_3D,&
			status,ierror)
		endif
	    endif
	    if(PBC_y)then
		call MPI_SENDRECV(Hx(0,N_y,0),1,TypeMatrXZ,neigborS,tagH3_NS,&
		    Hx(0,0,0),1,TypeMatrXZ,neigborN,tagH3_NS,comm_3D,status,&
		    ierror)
		call MPI_SENDRECV(Hz(0,N_y,0),1,TypeMatrXZ,neigborS,tagH4_NS,&
		    Hz(0,0,0),1,TypeMatrXZ,neigborN,tagH4_NS,comm_3D,status,&
		    ierror)
	    else
		if(ycord.gt.0)then
		    call MPI_Recv(Hx(0,0,0),1,TypeMatrXZ,neigborN,tagH3_NS,&
 			comm_3D,status,ierror)
		endif
		if(ycord.lt.(Ymax-1))then
		    call MPI_Send(Hx(0,N_y,0),1,TypeMatrXZ,neigborS,tagH3_NS,comm_3D,&
			status,ierror)
		endif
		if(ycord.gt.0)then
		    call MPI_Recv(Hz(0,0,0),1,TypeMatrXZ,neigborN,tagH4_NS,&
			comm_3D,status,ierror)
		endif
		if(ycord.lt.(Ymax-1))then
		    call MPI_Send(Hz(0,N_y,0),1,TypeMatrXZ,neigborS,tagH4_NS,comm_3D,&
			status,ierror)
		endif
	    endif

	    if(PBC_z)then
		call MPI_SENDRECV(Hx(0,0,N_z),countXY,MPI_double_precision,&
		    neigborE,tagH3_WE,Hx(0,0,0),countXY,MPI_double_precision,&
		    neigborW,tagH3_WE,comm_3D,status,ierror)
		call MPI_SENDRECV(Hy(0,0,N_z),countXY,MPI_double_precision,&
		    neigborE,tagH4_WE,Hy(0,0,0),countXY,MPI_double_precision,&
		    neigborW,tagH4_WE,comm_3D,status,ierror)
	    else
		if(zcord.gt.0)then
		    call MPI_Recv(Hx(0,0,0),countXY,MPI_double_precision,&
			neigborW,tagH3_WE,comm_3D,status,ierror)
		endif
		if(zcord.lt.Zmax-1)then
		    call MPI_Send(Hx(0,0,N_z),countXY,MPI_double_precision,&
			neigborE,tagH3_WE,comm_3D,status,ierror)
		endif

		if(zcord.gt.0)then
		    call MPI_Recv(Hy(0,0,0),countXY,MPI_double_precision,&
			neigborW,tagH4_WE,comm_3D,status,ierror)
		endif
		if(zcord.lt.Zmax-1)then
		    call MPI_Send(Hy(0,0,N_z),countXY,MPI_double_precision,neigborE,&
			tagH4_WE,comm_3D,status,ierror)
		endif
	    endif
!		second step of syncronization [send(F(1)), recv(F(N_i+1))
	    if(PBC_x)then
		call MPI_SENDRECV(Hy(1,0,0),1,TypeMatrYZ,neigborD,tagH1_UD,&
		    Hy(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagH1_UD,comm_3D,status,&
		    ierror)

		call MPI_SENDRECV(Hz(1,0,0),1,TypeMatrYZ,neigborD,tagH2_UD,&
		    Hz(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagH2_UD,comm_3D,status,&
		    ierror)
	    else
		if(xcord.gt.0)then
		    call MPI_Send(Hy(1,0,0),1,TypeMatrYZ,neigborD,tagH1_UD,&
			comm_3D,status,ierror)
		endif
		if(xcord.lt.(Xmax-1))then
		    call MPI_Recv(Hy(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagH1_UD,comm_3D,&
			status,ierror)
		endif
		if(xcord.gt.0)then
		    call MPI_Send(Hz(1,0,0),1,TypeMatrYZ,neigborD,tagH2_UD,&
			comm_3D,status,ierror)
		endif
		if(xcord.lt.(Xmax-1))then
		    call MPI_Recv(Hz(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagH2_UD,comm_3D,&
			status,ierror)
		endif
	    endif

	    if(PBC_y)then
		call MPI_SENDRECV(Hx(0,1,0),1,TypeMatrXZ,neigborN,tagH1_NS,&
		    Hx(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagH1_NS,comm_3D,status,&
		    ierror)
		call MPI_SENDRECV(Hz(0,1,0),1,TypeMatrXZ,neigborN,tagH2_NS,&
		    Hz(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagH2_NS,comm_3D,status,&
		    ierror)
	    else
		if(ycord.gt.0)then
		    call MPI_Send(Hx(0,1,0),1,TypeMatrXZ,neigborN,tagH1_NS,&
 			comm_3D,status,ierror)
		endif
		if(ycord.lt.(Ymax-1))then
		    call MPI_Recv(Hx(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagH1_NS,comm_3D,&
			status,ierror)
		endif
		if(ycord.gt.0)then
		    call MPI_Send(Hz(0,1,0),1,TypeMatrXZ,neigborN,tagH2_NS,&
			comm_3D,status,ierror)
		endif
		if(ycord.lt.(Ymax-1))then
		    call MPI_Recv(Hz(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagH2_NS,comm_3D,&
			status,ierror)
		endif
	    endif

	    if(PBC_z)then
		call MPI_SENDRECV(Hx(0,0,1),countXY,MPI_double_precision,&
		    neigborW,tagH1_WE,Hx(0,0,N_z+1),countXY,MPI_double_precision,&
		    neigborE,tagH1_WE,comm_3D,status,ierror)

		call MPI_SENDRECV(Hy(0,0,1),countXY,MPI_double_precision,&
		    neigborW,tagH2_WE,Hy(0,0,N_z+1),countXY,MPI_double_precision,&
		    neigborE,tagH2_WE,comm_3D,status,ierror)
	    else
		if(zcord.gt.0)then
		    call MPI_Send(Hx(0,0,1),countXY,MPI_double_precision,&
			neigborW,tagH1_WE,comm_3D,status,ierror)
		endif
		if(zcord.lt.Zmax-1)then
		    call MPI_Recv(Hx(0,0,N_z+1),countXY,MPI_double_precision,&
			neigborE,tagH1_WE,comm_3D,status,ierror)
		endif

		if(zcord.gt.0)then
		    call MPI_Send(Hy(0,0,1),countXY,MPI_double_precision,&
			neigborW,tagH2_WE,comm_3D,status,ierror)
		endif
		if(zcord.lt.Zmax-1)then
		    call MPI_Recv(Hy(0,0,N_z+1),countXY,MPI_double_precision,neigborE,&
			tagH2_WE,comm_3D,status,ierror)
		endif
	    endif
	if(n.eq.1)then
	    print *,"ending H interchange",myrank
	endif
	call mpi_barrier(mpi_comm_world,ierror)
	end subroutine interchange_H
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c		alternative E field synchonization calls
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine interchange_E2()
!c		firs step of syncronization [send(F(1)), recv(F(N_i+1))
!ccccccccccccccccccccccc
!c	for the first ranks along z
	if(zcord.eq.0)then
	    if(n.eq.1)then
	    write(*,"(a16,2I3)")"mpi(E) recv Exy:",myrank,neigborE
	    end if
	call MPI_Recv(Ex(0,0,N_z+1),countXY,MPI_double_precision,&
		neigborE,tagE1_WE,comm_3D,status,ierror)
	call MPI_Recv(Ey(0,0,N_z+1),countXY,MPI_double_precision,&
		neigborE,tagE2_WE,comm_3D,status,ierror)
	endif
!c	for the first ranks along y
	if(ycord.eq.0)then
	    if(n.eq.1)then
	    write(*,"(a16,2I3)")"mpi(E) recv Exz:",myrank,neigborS
	    end if
	call MPI_Recv(Ex(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagE1_NS,&
		comm_3D,status,ierror)
	call MPI_Recv(Ez(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagE2_NS,&
	comm_3D,status,ierror)
	end if
!c	for the first ranks along x
	if(xcord.eq.0)then
	    if(n.eq.1)then
	    write(*,"(a16,2I3)")"mip(E) recv Eyz:",myrank,neigborU
	    end if
	call MPI_Recv(Ey(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagE1_UD,&
		comm_3D,status,ierror)
	call MPI_Recv(Ez(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagE2_UD,&
		comm_3D,status,ierror)
	end if

!c	for the middle ranks
	if(neigborE.ge.0.and.neigborW.ge.0)then
	    if(n.eq.1)then
	    write(*,"(a20,2I3)")"mpi(E) sendRecv E-W:",myrank,neigborE
	    end if
	call MPI_SENDRECV(Ex(0,0,1),countXY,MPI_double_precision,&
	neigborW,tagE1_WE,Ex(0,0,N_z+1),countXY,MPI_double_precision,&
	neigborE,tagE1_WE,comm_3D,status,ierror)
	call MPI_SENDRECV(Ey(0,0,1),countXY,MPI_double_precision,&
	neigborW,tagE2_WE,Ey(0,0,N_z+1),countXY,MPI_double_precision,&
	neigborE,tagE2_WE,comm_3D,status,ierror)
	end if

	if(neigborN.ge.0.and.neigborS.ge.0)then
	    if(n.eq.1)then
	    print *,"mpi(E) sendRecv N-S:",myrank,neigborW,neigborE
	    end if
	call MPI_SENDRECV(Ex(0,1,0),1,TypeMatrXZ,neigborN,tagE1_NS,&
	Ex(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagE1_NS,comm_3D,status&
	,ierror)
	call MPI_SENDRECV(Ez(0,1,0),1,TypeMatrXZ,neigborN,tagE2_NS,&
	Ez(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagE2_NS,comm_3D,status,&
	ierror)
	end if

	if(neigborU.ge.0.and.neigborD.ge.0)then
	    if(n.eq.1)then
	    print *,"mpi(E) sendRecv U-D:",myrank,neigborU,neigborD
	    end if
	call MPI_SENDRECV(Ey(1,0,0),1,TypeMatrYZ,neigborD,tagE1_UD,&
	Ey(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagE1_UD,comm_3D,status,&
	ierror)
	call MPI_SENDRECV(Ez(1,0,0),1,TypeMatrYZ,neigborD,tagE2_UD,&
	Ez(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagE2_UD,comm_3D,status,&
	ierror)
	end if
!c	for the last ranks along z
	if(zcord.eq.Zmax-1)then
	    if(n.eq.1)then
	    write(*,"(a16,2I3)")"mpi(E) send Exy:",myrank,neigborW
	    end if
	call MPI_Send(Ex(0,0,1),countXY,MPI_double_precision,&
	neigborW,tagE1_WE,comm_3D,status,ierror)
	call MPI_Send(Ey(0,0,1),countXY,MPI_double_precision,neigborW,&
	tagE2_WE,comm_3D,status,ierror)
	end if
!c	for the last ranks along y
	if(ycord.eq.Ymax-1)then
	    if(n.eq.1)then
	    write(*,"(a16,2I3)")"mpi(E) send Exz:",myrank,neigborN
	    end if
	call MPI_Send(Ex(0,1,0),1,TypeMatrXZ,neigborN,tagE1_NS,comm_3D,&
		status,ierror)
	call MPI_Send(Ez(0,1,0),1,TypeMatrXZ,neigborN,tagE2_NS,comm_3D,&
		status,ierror)
	end if
!c	for the last ranks along x
	if(xcord.eq.Xmax-1)then
	    if(n.eq.1)then
	    write(*,"(a16,2I3)")"mpi(E) send Eyz:",myrank,neigborD
	    end if
	call MPI_Send(Ey(1,0,0),1,TypeMatrYZ,neigborD,tagE1_UD,comm_3D,&
		status,ierror)
	call MPI_Send(Ez(1,0,0),1,TypeMatrYZ,neigborD,tagE2_UD,comm_3D,&
		status,ierror)
	end if
	call mpi_barrier(mpi_comm_world,ierror)
	    if(n.eq.1)then
	    print *,"Esincronyzed1",myrank
	    end if
!ccccccccccccccccccccccc
!c	c	second step of syncronization [send(F(N_i)), recv(F(0))
!ccccccccccccccccccccccc
!c	for the first ranks along z
	if(zcord.eq.0)then
	    if(n.eq.1)then
	    write(*,"(a17,2I3)")"mpi(E2) send Exy:",myrank,neigborE
	    end if
	call MPI_Send(Ex(0,0,N_z),countXY,MPI_double_precision,&
	neigborE,tagE3_WE,comm_3D,status,ierror)
	call MPI_Send(Ey(0,0,N_z),countXY,MPI_double_precision,&
	neigborE,tagE4_WE,comm_3D,status,ierror)
	end if
!c	for the first ranks along y
	if(ycord.eq.0)then
	    if(n.eq.1)then
	    write(*,"(a17,2I3)")"mpi(E2) send Exz:",myrank,neigborS
	    end if
	call MPI_Send(Ex(0,N_y,0),1,TypeMatrXZ,neigborS,tagE3_NS,&
		comm_3D,status,ierror)
	call MPI_Send(Ez(0,N_y,0),1,TypeMatrXZ,neigborS,tagE4_NS,&
		comm_3D,status,ierror)
	end if
!c	for the first ranks along x
	if(xcord.eq.0)then
	    if(n.eq.1)then
	    write(*,"(a16,2I3)")"mip(E2) send Eyz:",myrank,neigborU
	    end if
	call MPI_Send(Ey(N_x,0,0),1,TypeMatrYZ,neigborU,tagE3_UD,&
		comm_3D,status,ierror)
	call MPI_Send(Ez(N_x,0,0),1,TypeMatrYZ,neigborU,tagE4_UD,&
		comm_3D,status,ierror)
	end if

!c	for the middle ranks
	if(neigborE.ge.0.and.neigborW.ge.0)then
	    if(n.eq.1)then
	    print *,"mpi(E2) sendRecv E-W:",myrank,neigborE,neigborW
	    end if
	call MPI_SENDRECV(Ex(0,0,N_z),countXY,MPI_double_precision,&
	neigborE,tagE3_WE,Ex(0,0,0),countXY,MPI_double_precision,&
	neigborW,tagE3_WE,comm_3D,status,ierror)

	call MPI_SENDRECV(Ey(0,0,N_z),countXY,MPI_double_precision,&
	neigborE,tagE4_WE,Ey(0,0,0),countXY,MPI_double_precision,&
	neigborW,tagE4_WE,comm_3D,status,ierror)
	end if

	if(neigborN.ge.0.and.neigborS.ge.0)then
	    if(n.eq.1)then
	    print *,"mpi(E2) sendRecv N-S:",myrank,neigborN,neigborS
	    end if
	call MPI_SENDRECV(Ex(0,N_y,0),1,TypeMatrXZ,neigborS,tagE3_NS,&
	Ex(0,0,0),1,TypeMatrXZ,neigborN,tagE3_NS,comm_3D,status&
	,ierror)

	call MPI_SENDRECV(Ez(0,N_y,0),1,TypeMatrXZ,neigborS,tagE4_NS,&
	Ez(0,0,0),1,TypeMatrXZ,neigborN,tagE4_NS,comm_3D,status,&
	ierror)
	end if

	if(neigborU.ge.0.and.neigborD.ge.0)then
	    if(n.eq.1)then
	    print *,"mpi(E2) sendRecv U-D:",myrank,neigborU,neigborD
	    end if
	call MPI_SENDRECV(Ey(N_x,0,0),1,TypeMatrYZ,neigborU,tagE3_UD,&
	Ey(0,0,0),1,TypeMatrYZ,neigborD,tagE3_UD,comm_3D,status,&
	ierror)

	call MPI_SENDRECV(Ez(N_x,0,0),1,TypeMatrYZ,neigborU,tagE4_UD,&
	Ez(0,0,0),1,TypeMatrYZ,neigborD,tagE4_UD,comm_3D,status,&
	ierror)
	end if
!c	for the last ranks along z
	if(zcord.eq.(Zmax-1))then
	    if(n.eq.1)then
	    write(*,"(a17,2I3)")"mpi(E2) recv Exy:",myrank,neigborW
	    end if
	call MPI_Recv(Ex(0,0,0),countXY,MPI_double_precision,&
	neigborW,tagE3_WE,comm_3D,status,ierror)
	call MPI_Recv(Ey(0,0,0),countXY,MPI_double_precision,&
	neigborW,tagE4_WE,comm_3D,status,ierror)
	end if
!c	for the last ranks along y
	if(ycord.eq.(Ymax-1))then
	    if(n.eq.1)then
	    write(*,"(a17,2I3)")"mpi(E2) recv Exz:",myrank,neigborN
	    end if
	call MPI_Recv(Ex(0,0,0),1,TypeMatrXZ,neigborN,tagE3_NS,comm_3D,&
		status,ierror)
	call MPI_Recv(Ez(0,0,0),1,TypeMatrXZ,neigborN,tagE4_NS,comm_3D,&
		status,ierror)
	end if
!	for the last ranks along x
	if(xcord.eq.(Xmax-1))then
	    if(n.eq.1)then
	    write(*,"(a17,2I3)")"mpi(E2) recv Eyz:",myrank,neigborD
	    end if
	call MPI_Recv(Ey(0,0,0),1,TypeMatrYZ,neigborD,tagE3_UD,comm_3D,&
		status,ierror)
	call MPI_Recv(Ez(0,0,0),1,TypeMatrYZ,neigborD,tagE4_UD,comm_3D,&
		status,ierror)
	end if
	call mpi_barrier(mpi_comm_world,ierror)
	if(n.eq.1)then
	    print *,"Esincronyzed2",myrank
	end if
	end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c		H field synchonization calls
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine interchange_H2()
!ccccccccccccccccccccccc
!c	firs step of syncronization  [send(F(N_i)), recv(F(0))
!ccccccccccccccccccccccc
!c	for the last ranks along z
	if(zcord.eq.Zmax-1)then
	    if(n.eq.1)then
	    print *,"mpi(H) recv Hxy:",neigborW,tagH1_WE,tagH2_WE
	    end if
	call MPI_Recv(Hx(0,0,0),countXY,MPI_double_precision,neigborW,&
		tagH1_WE,comm_3D,status,ierror)
	call MPI_Recv(Hy(0,0,0),countXY,MPI_double_precision,neigborW,&
		tagH2_WE,comm_3D,status,ierror)
	end if
!c	for the last ranks along y
	if(ycord.eq.Ymax-1)then
	    if(n.eq.1)then
	    print *,"mpi(H) recv Hxz:",neigborN,tagH1_NS,tagH2_NS
	    end if
	call MPI_Recv(Hx(0,0,0),1,TypeMatrXZ,neigborN,tagH1_NS,comm_3D,&
		status,ierror)
	call MPI_Recv(Hz(0,0,0),1,TypeMatrXZ,neigborN,tagH2_NS,comm_3D,&
		status,ierror)
	end if
!c	for the last ranks along x
	if(xcord.eq.Xmax-1)then
	    if(n.eq.1)then
	    print *,"mpi(H) recv Hyz:",neigborD,tagH1_UD,tagH2_UD
	    end if
	call MPI_Recv(Hy(0,0,0),1,TypeMatrYZ,neigborD,tagH1_UD,comm_3D,&
		status,ierror)
	call MPI_Recv(Hz(0,0,0),1,TypeMatrYZ,neigborD,tagH2_UD,comm_3D,&
		status,ierror)
	end if

!c	for the middle ranks
	if(neigborE.ge.0.and.neigborW.ge.0)then
	    if(n.eq.1)then
	    print *,"mpi(H) sendrecv Hxy:",neigborE,neigborW
	    end if
	call MPI_SENDRECV(Hx(0,0,N_z),countXY,MPI_double_precision,&
	neigborE,tagH1_WE,Hx(0,0,0),countXY,MPI_double_precision,&
	neigborW,tagH1_WE,comm_3D,status,ierror)

	call MPI_SENDRECV(Hy(0,0,N_z),countXY,MPI_double_precision,&
	neigborE,tagH2_WE,Hy(0,0,0),countXY,MPI_double_precision,&
	neigborW,tagH2_WE,comm_3D,status,ierror)
	end if

	if(neigborN.ge.0.and.neigborS.ge.0)then
	    if(n.eq.1)then
	    print *,"mpi(H) sendrecv Hxz:",neigborN,neigborS
	    end if
	call MPI_SENDRECV(Hx(0,N_y,0),1,TypeMatrXZ,neigborS,tagH1_NS,&
	Hx(0,0,0),1,TypeMatrXZ,neigborN,tagH1_NS,comm_3D,status,ierror)

	call MPI_SENDRECV(Hz(0,N_y,0),1,TypeMatrXZ,neigborS,tagH2_NS,&
	Hz(0,0,0),1,TypeMatrXZ,neigborN,tagH2_NS,comm_3D,status,ierror)
	end if

	if(neigborU.ge.0.and.neigborD.ge.0)then
	    if(n.eq.1)then
	    print *,"mpi(H) sendrecv Hyz:",neigborU,neigborD
	    end if
	call MPI_SENDRECV(Hy(N_x,0,0),1,TypeMatrYZ,neigborU,tagH1_UD,&
	Hy(0,0,0),1,TypeMatrYZ,neigborD,tagH1_UD,comm_3D,status,ierror)

	call MPI_SENDRECV(Hz(N_x,0,0),1,TypeMatrYZ,neigborU,tagH2_UD,&
	Hz(0,0,0),1,TypeMatrYZ,neigborD,tagH2_UD,comm_3D,status,ierror)
	end if

!c	for the first ranks along z
	if(zcord.eq.0)then
	    if(n.eq.1)then
	    print *,"mpi(H) send Hxy:",neigborE,tagH1_WE,tagH2_WE
	    end if
	call MPI_Send(Hx(0,0,N_z),countXY,MPI_double_precision,neigborE,&
	tagH1_WE,comm_3D,status,ierror)
	call MPI_Send(Hy(0,0,N_z),countXY,MPI_double_precision,neigborE,&
	tagH2_WE,comm_3D,status,ierror)
	end if
!c	for the first ranks along y
	if(ycord.eq.0)then
	    if(n.eq.1)then
	    print *,"mpi(H) send Hxz:",neigborS,tagH1_NS,tagH2_NS
	    end if
	call MPI_Send(Hx(0,N_y,0),1,TypeMatrXZ,neigborS,tagH1_NS,comm_3D,&
		status,ierror)
	call MPI_Send(Hz(0,N_y,0),1,TypeMatrXZ,neigborS,tagH2_NS,comm_3D,&
		status,ierror)
	end if
!c	for the first ranks along x
	if(xcord.eq.0)then
	    if(n.eq.1)then
	    print *,"mpi(H) send Hyz:",neigborU,tagH1_UD,tagH2_UD
	    end if
	call MPI_Send(Hy(N_x,0,0),1,TypeMatrYZ,neigborU,tagH1_UD,comm_3D,&
		status,ierror)
	call MPI_Send(Hz(N_x,0,0),1,TypeMatrYZ,neigborU,tagH2_UD,comm_3D,&
		status,ierror)
	end if
!	call mpi_barrier(mpi_comm_world,ierror)
	if(n.eq.1)then
	print *,"H syncronized1:",myrank
	end if
!ccccccccccccccccccccccc
!c	c	second step of syncronization [send(F(1)), recv(F(N_i+1))
!ccccccccccccccccccccccc
!c	for the last ranks along z
	if(zcord.eq.Zmax-1)then
	    if(n.eq.1)then
	    print *,"mpi(H2) send Hxy:",neigborW,tagH3_WE,tagH4_WE
	    end if
	call MPI_Send(Hx(0,0,1),countXY,MPI_double_precision,neigborW,&
	tagH3_WE,comm_3D,status,ierror)
	call MPI_Send(Hy(0,0,1),countXY,MPI_double_precision,neigborW,&
	tagH4_WE,comm_3D,status,ierror)
	end if
!c	for the last ranks along y
	if(ycord.eq.Ymax-1)then
	    if(n.eq.1)then
	    print *,"mpi(H2) send Hxz:",neigborN,tagH3_NS,tagH4_NS
	    end if
	call MPI_Send(Hx(0,1,0),1,TypeMatrXZ,neigborN,tagH3_NS,&
		comm_3D,status,ierror)
	call MPI_Send(Hz(0,1,0),1,TypeMatrXZ,neigborN,tagH4_NS,&
		comm_3D,status,ierror)
	end if
!c	for the last ranks along x
	if(xcord.eq.Xmax-1)then
	    if(n.eq.1)then
	    print *,"mpi(H2) send Hyz:",neigborD,tagH3_UD,tagH4_UD
	    end if
	call MPI_Send(Hy(1,0,0),1,TypeMatrYZ,neigborD,tagH3_UD,comm_3D,&
		status,ierror)
	call MPI_Send(Hz(1,0,0),1,TypeMatrYZ,neigborD,tagH4_UD,comm_3D,&
		status,ierror)
	end if

!c	for the middle ranks
	if(neigborE.ge.0.and.neigborW.ge.0)then
	    if(n.eq.1)then
	    print *,"mpi(H2) sendrecv Hxy:",neigborE,neigborW
	    end if
	call MPI_SENDRECV(Hx(0,0,1),countXY,MPI_double_precision,&
	neigborW,tagH3_WE,Hx(0,0,N_z+1),countXY,MPI_double_precision,&
	neigborE,tagH3_WE,comm_3D,status,ierror)

	call MPI_SENDRECV(Hy(0,0,1),countXY,MPI_double_precision,&
	neigborW,tagH4_WE,Hy(0,0,N_z+1),countXY,MPI_double_precision,&
	neigborE,tagH4_WE,comm_3D,status,ierror)
	end if

	if(neigborN.ge.0.and.neigborS.ge.0)then
	    if(n.eq.1)then
	    print *,"mpi(H2) sendrecv Hxz:",neigborN,neigborS
	    end if
	call MPI_SENDRECV(Hx(0,1,0),1,TypeMatrXZ,neigborN,tagH3_NS,&
	Hx(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagH3_NS,comm_3D,status,&
	ierror)

	call MPI_SENDRECV(Hz(0,1,0),1,TypeMatrXZ,neigborN,tagH4_NS,&
	Hz(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagH4_NS,comm_3D,status,&
	ierror)
	end if

	if(neigborU.ge.0.and.neigborD.ge.0)then
	    if(n.eq.1)then
	    print *,"mpi(H2) sendrecv Hyz:",neigborU,neigborD
	    end if
	call MPI_SENDRECV(Hy(1,0,0),1,TypeMatrYZ,neigborD,tagH3_UD,&
	Hy(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagH3_UD,comm_3D,status,&
	ierror)

	call MPI_SENDRECV(Hz(1,0,0),1,TypeMatrYZ,neigborD,tagH4_UD,&
	Hz(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagH4_UD,comm_3D,status,&
	ierror)
	end if

!c	for the first ranks along z
	if(zcord.eq.0)then
	    if(n.eq.1)then
	    print *,"mpi(H2) recv Hxy:",neigborE,tagH3_WE,tagH4_WE
	    end if
	call MPI_Recv(Hx(0,0,N_z+1),countXY,MPI_double_precision,&
	neigborE,tagH3_WE,comm_3D,status,ierror)
	call MPI_Recv(Hy(0,0,N_z+1),countXY,MPI_double_precision,&
	neigborE,tagH4_WE,comm_3D,status,ierror)
	end if
!c	for the first ranks along y
	if(ycord.eq.0)then
	    if(n.eq.1)then
	    print *,"mpi(H2) recv Hxz:",neigborS,tagH3_NS,tagH4_NS
	    end if
	call MPI_Recv(Hx(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagH3_NS,&
		comm_3D,status,ierror)
	call MPI_Recv(Hz(0,N_y+1,0),1,TypeMatrXZ,neigborS,tagH4_NS,&
		comm_3D,status,ierror)
	end if
!c	for the first ranks along x
	if(xcord.eq.0)then
	    if(n.eq.1)then
	    print *,"mpi(H) recv Hyz:",neigborU,tagH3_UD,tagH4_UD
	    end if
	call MPI_Recv(Hy(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagH3_UD,&
		comm_3D,status,ierror)
	call MPI_Recv(Hz(N_x+1,0,0),1,TypeMatrYZ,neigborU,tagH4_UD,&
		comm_3D,status,ierror)
	end if
	call mpi_barrier(mpi_comm_world,ierror)
	end subroutine
