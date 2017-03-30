
	subroutine task_start(rank,numtasks)
	integer rank,numtasks	
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end
	
c----------------------------------------------------------------------
	
	subroutine task_abort()
	print*,'Aborting the program...'
	stop
	end
c----------------------------------------------------------------------
	
	subroutine task_finish()
	print*,'program is finished...'
	stop
	end

c----------------------------------------------------------------------
        subroutine task_barrier()
	return
        end

c----------------------------------------------------------------------

	subroutine task_send_float(rank_to,buffer,length,tag,request)
	implicit none
	integer rank_to		! receiving task's rank
	real buffer(*)		! buffer of data
	integer length		! buffers' length
	integer tag		! tag of the message
	integer request		! request id
	print*, 'MPIsndf call from a single task program! Exiting...'
	stop
	end

c----------------------------------------------------------------------

	subroutine task_send_integer(rank_to,buffer,length,tag,request)
	implicit none
	integer rank_to		! receiving task's rank
	integer buffer(*)	! buffer of data
	integer length		! buffers' length
	integer tag		! tag of the message
	integer request
	print*, 'MPIsndi call from a single task program! Exiting...'
	stop
	end
	
c----------------------------------------------------------------------

	subroutine task_send_character
     &                            (rank_to,buffer,length,tag,request)
	implicit none
	integer rank_to		! receiving task's rank
	character*1 buffer(*)	! buffer of data
	integer length		! buffers' length
	integer tag		! tag of the message
	integer request
	print*, 'MPIsndi call from a single task program! Exiting...'
	stop
	end
	
c----------------------------------------------------------------------

        subroutine task_receive_float(buffer,length,request)
	real buffer(*)		! buffer of data
	integer length		! buffers' length
	integer request
	print*, 'MPIrcvf call from a single task program! Exiting...'
	stop
	end

c----------------------------------------------------------------------

        subroutine task_receive_charcater(buffer,length,request)
	character*1 buffer(*)	! buffer of data
	integer length		! buffers' length
	integer request
	print*, 'MPIrcvi call from a single task program! Exiting...'
	stop
	end

c----------------------------------------------------------------------

        subroutine task_receive_integer(buffer,length,request)
	integer buffer(*)	! buffer of data
	integer length		! buffers' length
	integer request
	print*, 'MPIrcvi call from a single task program! Exiting...'
	stop
	end

c----------------------------------------------------------------------
        subroutine task_wait(request,rank,tag)
	integer request
	integer rank, tag
	return
	end

c----------------------------------------------------------------------
        
        subroutine task_waitall(count,reqs,ranks,tags)
 	integer count,reqs(count)
	integer ranks(count),tags(count)
	return
	end

c----------------------------------------------------------------------
        subroutine task_test(request,flag,rank,tag)
	integer request
	integer rank, tag
	logical flag
	print*, 'MPItst call from a single task program! Exiting...'
	stop
	end

c----------------------------------------------------------------------

        subroutine task_sum_real(buffer1,buffer2,length)
	real buffer1(*)	! buffer of data
	real buffer2(*)	! buffer of data
	integer length		! buffers' length
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end

c----------------------------------------------------------------------

        subroutine task_sum_real8(buffer1,buffer2,length)
	real buffer1(*)	! buffer of data
	real buffer2(*)	! buffer of data
	integer length		! buffers' length
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end
c----------------------------------------------------------------------

        subroutine task_sum_integer(buffer1,buffer2,length)
	real buffer1(*)	! buffer of data
	real buffer2(*)	! buffer of data
	integer length		! buffers' length
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end
c----------------------------------------------------------------------

        subroutine task_max_real(buffer1,buffer2,length)
	real buffer1(*)	! buffer of data
	real buffer2(*)	! buffer of data
	integer length		! buffers' length
	return
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end
c----------------------------------------------------------------------

        subroutine task_max_integer(buffer1,buffer2,length)
	real buffer1(*)	! buffer of data
	real buffer2(*)	! buffer of data
	integer length		! buffers' length
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end
c----------------------------------------------------------------------

        subroutine task_min_real(buffer1,buffer2,length)
	real buffer1(*)	! buffer of data
	real buffer2(*)	! buffer of data
	integer length		! buffers' length
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end
c----------------------------------------------------------------------

        subroutine task_min_integer(buffer1,buffer2,length)
	real buffer1(*)	! buffer of data
	real buffer2(*)	! buffer of data
	integer length		! buffers' length
	print*, 'MPI call from a single task program! Exiting...'
	stop
	end
c------------------------------------------------------------------
       subroutine task_rank_to_index (rank,i,j)
        integer rank, i, j
	i=0
	j=0
	return
	end
