program litio_variacional_nr_repl
  use nr_powell_compat
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: ndim = 2
  real(dp), parameter :: ftol = 1.0d-8
  real(dp), parameter :: hartree_to_ev = 27.211386245988d0

  integer :: i, iter
  real(dp) :: fret, p(ndim), xi(ndim, ndim), z
  real(dp), external :: ip

  common /core_z/ z

  xi = 0.0d0
  xi(1,1) = 1.0d0
  xi(2,2) = 1.0d0
  p = [1.5d0, 1.5d0]

  print *, 'CUAL ES EL VALOR DE Z '
  read(*,*) z

  call powell(p, xi, ndim, ftol, iter, fret, func_vec)

  print *, 'ITERACIONES:        = ', iter
  print *, 'ALFA, BETA          = ', (p(i), i=1, ndim)
  print *, 'ENER. ESTADO BASE   = ', hartree_to_ev * fret, ' eV'
  print *, 'ENER. IONIZACION    = ', hartree_to_ev * ip(p(1), p(2)), ' eV'

contains

  real(dp) function func_vec(x, n)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp) :: a, b, y, z

    common /core_z/ z

    a = x(1)
    b = x(2)
    y = 2.0d0 * a * b / ((2.0d0 * a + b)**5)

    func_vec = a*a - 2.0d0*z*a + (5.0d0/8.0d0)*a + (1.0d0/8.0d0)*b*b - (z*b/4.0d0) + &
               y*(8.0d0*a**4 + 20.0d0*a**3*b + 12.0d0*a*a*b*b + 10.0d0*a*b**3 + b**4)
  end function func_vec

end program litio_variacional_nr_repl

real(kind(1.0d0)) function ip(alfa, beta)
  implicit none

  real(kind(1.0d0)), intent(in) :: alfa, beta
  real(kind(1.0d0)) :: a, b, a2, b2, y, z

  common /core_z/ z

  a = alfa
  b = beta
  a2 = a*a
  b2 = b*b
  y = 2.0d0 * a * b / ((2.0d0*a + b)**5)

  ip = (1.0d0/8.0d0)*b*b - (z*b/4.0d0) + &
       y*(8.0d0*a2*a2 + 20.0d0*a2*a*b + 12.0d0*a2*b2 + 10.0d0*a*b*b2 + b2*b2)
end function ip
