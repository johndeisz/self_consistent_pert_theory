MODULE lattice

  integer :: N_dim
  double precision :: a(3,3)
  integer :: Nl(3)
  
  integer :: N_a
  double precision, dimension(:,:), allocatable :: r_atom

end MODULE lattice
