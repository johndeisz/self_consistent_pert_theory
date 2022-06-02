MODULE Hamiltonian

  integer :: N_dim
  double precision :: a(3,3)
  
  integer :: N_a
  double precision, dimension(:,:), allocatable :: r_atom

  double precision, dimension(:), allocatable :: ed
  integer, dimension(:), allocatable :: N_o

end MODULE Hamiltonian
