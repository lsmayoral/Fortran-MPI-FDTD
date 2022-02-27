!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Initialize MPI environment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine my_mpi_init()
	call mpi_init(ierror)
	call mpi_comm_size(mpi_comm_world,nproces,ierror)
	call mpi_comm_rank(mpi_comm_world,myrank,ierror)

	tagE_f=1
	tagE_b=2
	tagH_f=3
	tagH_b=4
	call set_index()
	END subroutine my_mpi_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Set local and absolute indixes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine set_index()

        absIrange%init=myrank*n_pos+1
        absIrange%ends=absIrange%init+n_pos-1

        locIrange%init=1
        locIrange%ends=n_pos
	if(myrank.eq.0)then
	    locIrange%init=npml+1
	endif
	if(myrank.eq.(NCmax-1))then
	    locIrange%ends=N_pos-npml
	endif

        if(myrank.le.9)then
	    write(*,"(a10,I2,6I6)")"locRange:",myrank,locIrange%init,locIrange%ends
	    write(*,"(a10,I2,6I6)")"absRange:",myrank,absIrange%init,absIrange%ends

        else if(myrank.ge.10.and.myrank.le.99)then
	    write(*,"(a10,I4,6I6)")"locRange:",locIrange%init,locIrange%ends
	    write(*,"(a10,I3,6I6)")"absRange:",absIrange%init,absIrange%ends
        else
	    write(*,"(a10,I4,6I6)")"locRange:",myrank,locIrange%init,locIrange%ends
	    write(*,"(a10,I4,6I6)")"absRange:",myrank,absIrange%init,absIrange%ends
        endif
	end subroutine set_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		E field synchonization calls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c		E field synchonization calls
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine interchange_E(E)
	real*8, intent(in), pointer :: E(:)
	if(n.eq.1)then
	    print *,"starting E interchange",myrank
	endif
!		firs step of syncronization [send(F(1)), recv(F(N_i+1))
	if(myrank.gt.0)then
	    if(n.eq.1)then
		print *,"SendE_b",myrank
	    endif
	    call MPI_Send(E(1),1,MPI_double_precision,myrank-1,tagE_b,&
		mpi_comm_world,status,ierror)
	endif
	if(myrank.lt.(NCmax-1))then
	    if(n.eq.1)then
		print *,"RecvE_b",myrank
	    endif
	    call MPI_Recv(E(N_pos+1),1,MPI_double_precision,myrank+1,tagE_b,&
		mpi_comm_world,status,ierror)
	endif
!	second step of syncronization [send(F(N_i)), recv(F(0))
	if(myrank.gt.0)then
	    if(n.eq.1)then
	        print *,"RecvE_f",myrank
	    endif
	    call MPI_Recv(E(0),1,MPI_double_precision,myrank-1,tagE_f,&
		mpi_comm_world,status,ierror)
	endif
	if(myrank.lt.(NCmax-1))then
	    if(n.eq.1)then
		print *,"SendE_f",myrank
	    endif
	    call MPI_Send(E(N_pos),1,MPI_double_precision,myrank+1,tagE_f,&
		mpi_comm_world,status,ierror)
	endif
	call mpi_barrier(mpi_comm_world,ierror)
	end subroutine interchange_E
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		H field synchonization calls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine interchange_H(H)
	real*8, intent(in), pointer :: H(:)
	if(n.eq.1)then
	    print *,"starting H interchange",myrank
	endif
!	first step of syncronization [send(F(N_i)), recv(F(0))
	if(myrank.gt.0)then
	    call MPI_Recv(H(0),1,MPI_double_precision,myrank-1,tagH_b,&
		mpi_comm_world,status,ierror)
	endif
	if(myrank.lt.(NCmax-1))then
	    call MPI_Send(H(N_pos),1,MPI_double_precision,myrank+1,tagH_b,&
		mpi_comm_world,status,ierror)
	endif
!		second step of syncronization [send(F(1)), recv(F(N_i+1))
	if(myrank.gt.0)then
	    call MPI_Send(H(1),1,MPI_double_precision,myrank-1,tagH_f,&
		    mpi_comm_world,status,ierror)
	endif
	if(myrank.lt.(NCmax-1))then
	    call MPI_Recv(H(N_pos+1),1,MPI_double_precision,myrank+1,tagH_f,&
	    mpi_comm_world,status,ierror)
	endif
	call mpi_barrier(mpi_comm_world,ierror)
	end subroutine interchange_H
