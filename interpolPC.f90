module interpol
  use data_grid

  implicit none

  public :: read_input, diab_energy, diab_coupling

  private
  integer, parameter :: N_input=212 ! number of data points for total energies and the nonadiabatic coupling vector
  double precision :: R_grid(N_input)
  double precision :: Vdia(N_input,N_states),Vdia_NEW(N_input,N_states) ! total energy curves in the diabatic basis
  double precision :: coupling(N_input,N_states,N_states),coupling_NEW(N_input,N_states,N_states) ! diabatic couplings
contains

!---------- reads input files for total energies
  subroutine read_input
    implicit none

    double precision :: dydR1,dydR2
    integer :: i,j

    Vdia(:,:) = 0.d0
    coupling(:,:,:) = 0.d0

    !read diabatic potential energy surfaces
    open(98,file='Vdiabatic.inp',status='old')
    do i=1,N_input
       read(98,*) R_grid(i),Vdia(i,:)
    end do
    close(98)

    !read diabatic couplings
    open(98,file='diab_coupling.inp',status='old')
    do i=1,N_input
       read(98,*) R_grid(i), &
                  coupling(i,1,2)
    end do
    close(98)

    !preprocess the total energy curves in the diabatic basis
    do i=1,N_states
       dydR1=(Vdia(2,i)-Vdia(1,i))/(R_grid(2)-R_grid(1))
       dydR2=(Vdia(N_input,i)-Vdia(N_input-1,i))/(R_grid(N_input)-R_grid(N_input-1))
       call spline(R_grid,Vdia(:,i),N_input,dydR1,dydR2,Vdia_NEW(:,i))
    end do

    !preprocess the total energy curves in the diabatic basis
    do i=1,N_states-1
       do j=i+1,N_states
          dydR1=(coupling(2,i,j)-coupling(1,i,j))/(R_grid(2)-R_grid(1))
          dydR2=(coupling(N_input,i,j)-coupling(N_input-1,i,j))/(R_grid(N_input)-R_grid(N_input-1))
          call spline(R_grid,coupling(:,i,j),N_input,dydR1,dydR2,coupling_NEW(:,i,j))
       end do
    end do

  end subroutine read_input


!---------- interpolates the diabatic potential energy surfaces
  function diab_energy(r,state)
    implicit none
    
    integer :: state
    double precision :: r,diab_energy,help

    if (r.lt.R_grid(1)) then
       print*, "WARNING: r is lower than smallest point in R_grid"
       stop
    end if

    if (r.gt.R_grid(N_input)) then
       diab_energy=Vdia(N_input,state)
       return
    end if

    call splint(R_grid,Vdia(:,state),Vdia_NEW(:,state),N_input,r,help)
    diab_energy=help

  end function diab_energy



!---------- interpolates the diabatic coupling matrix elements
  function diab_coupling(r,from,to)
    implicit none
    
    integer :: from,to
    double precision :: r,diab_coupling,help

    if (from.ge.to) then
       print*, "WARNING: coupling_index 1 is larger then index 2"
       stop
    end if

    if (r.lt.R_grid(1)) then
       print*, "WARNING: r is lower than smallest point in R_grid"
       diab_coupling=coupling(1,from,to)
       return
!       stop
    end if

    if (r.gt.R_grid(N_input)) then
       diab_coupling=coupling(N_input,from,to)
       return
    end if

    call splint(R_grid,coupling(:,from,to),coupling_NEW(:,from,to),N_input,r,help)
    diab_coupling=help

  end function diab_coupling



!---------- preprocesses the input data for interpolation with subroutine splint
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500000)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do
      return
      end SUBROUTINE spline


!---------- interpolation routine
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      return
      end SUBROUTINE splint


end module interpol
