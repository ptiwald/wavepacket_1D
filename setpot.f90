  subroutine potential
         
  use data_grid
  use data_au
  use interpol
 
  implicit none 
  
  integer :: i,j,k
  double precision :: R

  call read_input

  ! set Pot, the diabatic potential energy surfaces
  do j = 1, N_states
     do i = 1, NR
        r = r0 + (i-1) *dr
        pot(i,j) = diab_energy(r,j)
     end do
  end do

  ! diab_coupling_matrix, the off-diagonal part of the diabatic potential matrix
  diab_coupling_matrix(:,:,:) = 0.d0
  do j = 1, N_states-1
     do k = j+1, N_states
        do i = 1, nr
           r = r0 + (i-1) *dr
           diab_coupling_matrix(i,j,k) = diab_coupling(r,j,k)
           diab_coupling_matrix(i,k,j) = diab_coupling_matrix(i,j,k)
        end do
     end do
  end do

  open(10,file='potentials.out',status='unknown')
  open(23,file='coupling_1X.out',status='unknown')

  do I = 1, NR	  
    R = R0 + (I-1) * dR     			  
    write(10,"(F15.8,F15.8,F15.8)") R, pot(i,:)
    write(23,"(F15.8,F15.8)") & 
         R, diab_coupling_matrix(i,1,2)
  end do 
  
  close(10,status='keep')   
  close(23,status='keep')

return
end subroutine
