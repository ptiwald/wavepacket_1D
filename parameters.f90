module data_grid
  public
  integer,parameter :: N_states = 2             ! Anzahl der elektron. Zustaende
  integer :: Nt        				! Zeitschritte, Anzahl der vib. Zustaende
  integer :: NR                                 ! Anzahl der Gridpunkte im Ortsraum
  integer :: initial_state                      ! Index des Anfangszustandes
  double precision:: R_init      		! Anfangskernabstand
  double precision:: P_init                     ! Anfangsimpuls
  double precision:: kappa			! FWHM des Anfangswellenpakets im Impulsraum
  double precision:: dt				! Zeitschritt
  double precision:: dpr, dr			! Aufloesung des Impuls- und Ortsgrids
  double precision:: R0                         ! Anfangspunkt des Ortsgrids
  double precision:: Rend                       ! Endpunkt des Ortsgrids
  double precision, allocatable, dimension(:) :: PR     ! Impulsgrid
  double precision, allocatable, dimension(:,:) :: Pot  ! Potential
  double precision, allocatable, dimension(:,:,:) :: diab_coupling_matrix ! off-diagonal part of diabatic potential matrix
end module data_grid


module data_au
  public
  double precision,parameter:: pi=3.141592653589793d0       ! that's just pi
  double precision,parameter:: Mass=6382.9491d0
  complex*16,parameter:: im=(0.d0,1.d0)   
end module data_au
