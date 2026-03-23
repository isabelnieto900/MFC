module li_spectrum_fd
  implicit none
  private
  public :: potential_li, determinant_poly, bracket_roots, bisect_root, inverse_iteration

  integer, parameter :: dp = kind(1.0d0)

contains

  pure function potential_li(r, l, alpha) result(v)
    real(dp), intent(in) :: r, alpha
    integer, intent(in) :: l
    real(dp) :: v

    v = -1.0_dp / r - 2.0_dp / r * (1.0_dp + alpha * r) * exp(-2.0_dp * alpha * r) + &
      real(l * (l + 1), dp) / (2.0_dp * r * r)
  end function potential_li

  function determinant_poly(e, v, h) result(pn)
    real(dp), intent(in) :: e, h
    real(dp), intent(in) :: v(:)
    real(dp) :: pn

    integer :: i, n
    real(dp) :: h2, ai, p_im2, p_im1, p_i

    n = size(v)
    h2 = h * h

    ai = 2.0_dp + 2.0_dp * h2 * v(1) - 2.0_dp * h2 * e
    p_im2 = 0.0_dp
    p_im1 = ai

    if (n == 1) then
      pn = p_im1
      return
    end if

    do i = 2, n
      ai = 2.0_dp + 2.0_dp * h2 * v(i) - 2.0_dp * h2 * e
      p_i = ai * p_im1 - p_im2

      if (abs(p_i) > 1.0d100) then
        p_i = p_i * 1.0d-100
        p_im1 = p_im1 * 1.0d-100
        p_im2 = p_im2 * 1.0d-100
      end if

      p_im2 = p_im1
      p_im1 = p_i
    end do

    pn = p_im1
  end function determinant_poly

  subroutine bracket_roots(v, h, emin, emax, nscan, xb1, xb2, nroots)
    real(dp), intent(in) :: v(:), h, emin, emax
    integer, intent(in) :: nscan
    real(dp), intent(out) :: xb1(:), xb2(:)
    integer, intent(out) :: nroots

    integer :: i, nbmax
    real(dp) :: e_prev, e_curr, f_prev, f_curr, de

    nroots = 0
    nbmax = min(size(xb1), size(xb2))

    de = (emax - emin) / real(nscan, dp)
    e_prev = emin
    f_prev = determinant_poly(e_prev, v, h)

    do i = 1, nscan
      e_curr = emin + de * real(i, dp)
      f_curr = determinant_poly(e_curr, v, h)

      if (f_prev == 0.0_dp) then
        if (nroots < nbmax) then
          nroots = nroots + 1
          xb1(nroots) = e_prev - de
          xb2(nroots) = e_prev + de
        end if
      else if (f_prev * f_curr < 0.0_dp) then
        if (nroots < nbmax) then
          nroots = nroots + 1
          xb1(nroots) = e_prev
          xb2(nroots) = e_curr
        end if
      end if

      e_prev = e_curr
      f_prev = f_curr
    end do
  end subroutine bracket_roots

  function bisect_root(v, h, a, b, tol, max_iter) result(root)
    real(dp), intent(in) :: v(:), h, a, b, tol
    integer, intent(in) :: max_iter
    real(dp) :: root

    integer :: it
    real(dp) :: left, right, mid, f_left, f_mid

    left = a
    right = b
    f_left = determinant_poly(left, v, h)

    do it = 1, max_iter
      mid = 0.5_dp * (left + right)
      f_mid = determinant_poly(mid, v, h)

      if (abs(f_mid) < tol .or. abs(right - left) < tol) then
        root = mid
        return
      end if

      if (f_left * f_mid <= 0.0_dp) then
        right = mid
      else
        left = mid
        f_left = f_mid
      end if
    end do

    root = 0.5_dp * (left + right)
  end function bisect_root

  subroutine solve_tridiagonal(a, b, c, rhs, x)
    real(dp), intent(in) :: a(:), b(:), c(:), rhs(:)
    real(dp), intent(out) :: x(:)

    integer :: n, i
    real(dp), allocatable :: cp(:), dpv(:)
    real(dp) :: denom

    n = size(b)
    allocate(cp(n-1), dpv(n))

    denom = b(1)
    if (abs(denom) < 1.0d-14) denom = sign(1.0d-14, denom + 1.0d-30)
    cp(1) = c(1) / denom
    dpv(1) = rhs(1) / denom

    do i = 2, n - 1
      denom = b(i) - a(i-1) * cp(i-1)
      if (abs(denom) < 1.0d-14) denom = sign(1.0d-14, denom + 1.0d-30)
      cp(i) = c(i) / denom
      dpv(i) = (rhs(i) - a(i-1) * dpv(i-1)) / denom
    end do

    denom = b(n) - a(n-1) * cp(n-1)
    if (abs(denom) < 1.0d-14) denom = sign(1.0d-14, denom + 1.0d-30)
    dpv(n) = (rhs(n) - a(n-1) * dpv(n-1)) / denom

    x(n) = dpv(n)
    do i = n - 1, 1, -1
      x(i) = dpv(i) - cp(i) * x(i + 1)
    end do

    deallocate(cp, dpv)
  end subroutine solve_tridiagonal

  subroutine normalize_radial(h, vec)
    real(dp), intent(in) :: h
    real(dp), intent(inout) :: vec(:)

    real(dp) :: nrm

    nrm = sqrt(sum(vec * vec) * h)
    if (nrm > 0.0_dp) vec = vec / nrm
  end subroutine normalize_radial

  subroutine inverse_iteration(v, eig, h, n_iter, vec)
    real(dp), intent(in) :: v(:), eig, h
    integer, intent(in) :: n_iter
    real(dp), intent(out) :: vec(:)

    integer :: n, k
    real(dp), allocatable :: a(:), b(:), c(:), rhs(:), x(:)
    real(dp) :: shift

    n = size(v)
    allocate(a(n-1), b(n), c(n-1), rhs(n), x(n))

    a = -1.0_dp
    c = -1.0_dp

    rhs = 0.0_dp
    rhs(1) = 1.0_dp
    call normalize_radial(h, rhs)

    shift = 1.0d-10
    do k = 1, n_iter
      b = 2.0_dp + 2.0_dp * h * h * v - 2.0_dp * h * h * (eig + shift)
      call solve_tridiagonal(a, b, c, rhs, x)
      call normalize_radial(h, x)
      rhs = x
    end do

    vec = rhs
    deallocate(a, b, c, rhs, x)
  end subroutine inverse_iteration

end module li_spectrum_fd

program litio_espectro_fd_actualizado
  use li_spectrum_fd
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: hartree_to_ev = 27.211386245988_dp

  integer :: npt, l, i, j, nscan, nroots, max_states
  real(dp) :: dmax, h, alpha, emin, emax, tol, eig
  real(dp), allocatable :: r(:), v(:), eigvec(:), xb1(:), xb2(:)

  npt = 2500
  dmax = 120.0_dp
  alpha = 2.535930_dp
  nscan = 15000
  max_states = 40
  emin = -0.40_dp
  emax = -1.0d-4
  tol = 1.0d-12

  print *, 'Momentum angular l (0=s, 1=p, 2=d, ...):'
  read(*, *) l

  h = dmax / real(npt + 1, dp)

  allocate(r(npt), v(npt), eigvec(npt))
  allocate(xb1(max_states), xb2(max_states))

  do i = 1, npt
    r(i) = h * real(i, dp)
    v(i) = potential_li(r(i), l, alpha)
  end do

  call bracket_roots(v, h, emin, emax, nscan, xb1, xb2, nroots)

  if (nroots == 0) then
    print *, 'No se encontraron autovalores en el rango dado.'
    stop
  end if

  open(unit=10, file='litio_energias.dat', status='replace', action='write')
  open(unit=11, file='litio_funciones_onda.dat', status='replace', action='write')

  write(10, '(a)') '# idx  l   E(Ha)            E(eV)'

  do j = 1, nroots
    eig = bisect_root(v, h, xb1(j), xb2(j), tol, 200)
    write(10, '(i4,1x,i2,1x,es20.12,1x,es20.12)') j, l, eig, eig * hartree_to_ev

    call inverse_iteration(v, eig, h, 12, eigvec)

    write(11, '(a,i0,a,i0,a,es20.12,a)') '# estado=', j, ' l=', l, ' E(Ha)=', eig, ''
    do i = 1, npt, 5
      write(11, '(f14.8,1x,es20.12)') r(i), eigvec(i)
    end do
    write(11, *)
  end do

  close(10)
  close(11)

  print *, 'Autovalores encontrados: ', nroots
  print *, 'Archivo de energias: litio_energias.dat'
  print *, 'Archivo de funciones: litio_funciones_onda.dat'

  deallocate(r, v, eigvec, xb1, xb2)
end program litio_espectro_fd_actualizado
