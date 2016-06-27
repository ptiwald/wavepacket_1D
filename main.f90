program paul_dynamik

  use data_grid
  use data_au
  
  implicit none
  integer:: I
  
  call input				! read parameters from input file
  allocate(PR(NR), Pot(NR,N_states), diab_coupling_matrix(NR,N_states,N_states))      ! allocate memory for momentum grid, diabatic potential energy surfaces
                                                                                 ! and the offdiagonal part of the diabatic matrix (i.e. the couplings)
  call p_grid				! initialize momentum grid

  call potential			! initialize diabatic potential energy surfaces and the couplings
  
  write(*,*)
  write(*,*) 'INITIAL KINETIC ENERGY [a.u.] : ', P_init**2/Mass*0.5d0
  write(*,*) 'REDUCED MASS [a.u.] : ', Mass
  write(*,*)
  
  call propagation_1D	! Propagation
  
   deallocate (pr, Pot, diab_coupling_matrix)
  
  stop	  
end program paul_dynamik


! _______________ Subroutines __________________________________________________


subroutine input

use data_grid
use data_au

implicit none
integer :: getpid

 open(10,file='input',status='old')
 
  read(10,*) dt			  ! dt = time step (a.u.).
  read(10,*) Nt			  ! Nt = number of time steps.	
  read(10,*) R_init		  ! RI = center of initial Gaussian (a.u.)
  read(10,*) kappa		  ! kappa = FWHM of initial Gaussian in momentum space (a.u.)
  read(10,*) P_init               ! PI = center of inial Gaussian in momentum space (a.u.)
  read(10,*) NR                   ! NR = number of grid points
  read(10,*) R0                   ! R0 = starting point of grid (a.u.)
  read(10,*) Rend                 ! Rend = end point of grid (a.u.)
  read(10,*) initial_state        ! index of initial state

  dR = (Rend - R0) / (NR - 1)
 
  R_init = R_init
  
  dpr = (2.d0 * pi) / (dR * NR)  

  print*
  print*, "The current process ID is ", getpid()
  print* 
  print*,'_________________________'
  print*
  print*,'Parameters'
  print*,'_________________________'
  print*
  print*,'dt  = ', SNGL(dt), 'a.u.'
  print*,'NR  = ', NR
  print*,'L   = ', SNGL(Rend - R0), 'a.u.'
  print*,'dpR = ', SNGL(dpR), 'a.u.'
  print*,'dR  = ', SNGL(dR), 'a.u.'
  print*
  print*,'__________________________'
  print*		  
                 
 end subroutine
 
!...................... Impulsgrid......................
 
subroutine p_grid

use data_grid
use data_au
 implicit none 
 integer:: I 
      
  dpr = (2.d0 * pi) / (dR * NR)  

  do I = 1, NR  
    if (I.le.(NR / 2)) then    
    PR(I) = (I - 1) * dpR    
    else    
    PR(I) = - (NR + 1 - I) * dpR    
    end if    
  end do
           
  return
end subroutine  
