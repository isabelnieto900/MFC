program litio_espectro_nr_repl
  use nr_compat_litio
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: npt_const = 6000, nbmax = 1000
  real(dp), parameter :: hartree_to_ev = 27.211386245988d0

  integer :: i, j, k, kk, jj, nscan, nb, nn
  real(dp) :: xb1(nbmax), xb2(nbmax), c(npt_const)
  real(dp) :: energia(nbmax), diag(npt_const), superd(npt_const), subd(npt_const)
  real(dp) :: vk(npt_const), vkm1(npt_const)
  real(dp) :: h, ele, d, x1, x2, tol, raiz, vp

  real(dp), external :: determ, pot

  common /spec_params/ h, ele, nn

  print *, 'MOMENTUM ANGULAR (l) = '
  read(*,*) ele
  print *, 'DISTANCIA MAXIMA D (u.a.) = '
  read(*,*) d

  nn = npt_const
  h = d / real(nn, dp)

  x1 = -10.0d0
  x2 = 0.0d0
  nscan = 2000
  nb = nbmax

  call zbrak(determ, x1, x2, nscan, xb1, xb2, nb)

  jj = 0
  do i = 1, nb
    tol = 1.0d-12 * abs(0.5d0 * (xb1(i) + xb2(i)))
    if (tol == 0.0d0) tol = 1.0d-12

    raiz = zbrent(determ, xb1(i), xb2(i), tol)
    jj = jj + 1
    energia(jj) = raiz
    print *, 'VALOR PROPIO', jj, hartree_to_ev * energia(jj), determ(raiz)
  end do

  open(unit=10, file='LITIO.DAT', status='replace', action='write')

  do i = 1, nn
    c(i) = pot(h * real(i, dp))
  end do

  do j = 1, jj
    vp = energia(j)

    do k = 1, nn
      diag(k) = 2.0d0 + 2.0d0*h*h*c(k) - 2.0d0*h*h*vp
    end do

    do k = 1, nn - 1
      superd(k) = -1.0d0
    end do

    subd(1) = 0.0d0
    do k = 2, nn
      subd(k) = -1.0d0
    end do

    vk = 0.0d0
    vk(1) = 1.0d0

    do kk = 1, 10
      call tridag(subd, diag, superd, vk, vkm1, nn)
      call norma(nn, vkm1, vk)
    end do

    do i = 1, nn, 10
      write(10,*) real(i, dp) * h, vk(i)
    end do
    write(10,*)
  end do

  close(10)

contains

  subroutine norma(n, vin, vout)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: vin(n)
    real(dp), intent(out) :: vout(n)
    real(dp) :: nrm

    nrm = sqrt(sum(vin*vin))
    if (nrm > 0.0d0) then
      vout = vin / nrm
    else
      vout = vin
    end if
  end subroutine norma

end program litio_espectro_nr_repl

real(kind(1.0d0)) function pot(x)
  implicit none

  real(kind(1.0d0)), intent(in) :: x
  real(kind(1.0d0)) :: h, ele
  integer :: nn
  real(kind(1.0d0)), parameter :: alfa = 2.535930d0

  common /spec_params/ h, ele, nn

  pot = -(1.0d0/x) - (2.0d0/x)*(1.0d0 + alfa*x)*exp(-2.0d0*alfa*x)
  pot = pot + ele*(ele + 1.0d0)/(2.0d0*x*x)
end function pot

real(kind(1.0d0)) function determ(x)
  implicit none

  real(kind(1.0d0)), intent(in) :: x
  real(kind(1.0d0)) :: h, ele
  integer :: nn, i
  real(kind(1.0d0)) :: h2, p0, p1, p2
  real(kind(1.0d0)), external :: pot

  common /spec_params/ h, ele, nn

  h2 = h*h
  p0 = 0.0d0
  p1 = 2.0d0 + 2.0d0*h2*pot(h) - 2.0d0*x*h2

  do i = 2, nn
    p2 = (2.0d0 + 2.0d0*h2*pot(h*real(i,kind(1.0d0))) - 2.0d0*h2*x)*p1 - p0
    p0 = p1
    p1 = p2
  end do

  determ = p2
end function determ
