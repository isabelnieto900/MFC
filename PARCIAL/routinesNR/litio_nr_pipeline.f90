program litio_nr_pipeline
  use nr_powell_compat
  use nr_compat_litio
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: npt = 16000, nbmax = 1000
  real(dp), parameter :: hartree_to_ev = 27.211386245988d0

  integer :: i, nfound_s, nfound_p, nfound_d
  real(dp) :: z
  real(dp) :: alpha0, beta0, alpha_opt, beta_opt, e0, ei
  integer :: it

  real(dp) :: r(npt), vc(npt), veff(npt)
  real(dp) :: rs2(npt), rs3(npt), rs4(npt), rp2(npt), rp3(npt), rd3(npt)
  real(dp) :: e_s(nbmax), e_p(nbmax), e_d(nbmax)
  real(dp), save :: h_spec, ele_spec, z_spec, alpha_spec
  integer, save :: nn_spec

  common /var_z/ z

  print *, 'Valor de Z (recomendado 3 para Li):'
  read(*,*) z

  alpha0 = 1.5d0
  beta0 = 1.5d0
  call optimize_variational(z, alpha0, beta0, alpha_opt, beta_opt, e0, it)
  ei = ionization_energy_li(alpha_opt, beta_opt, z)

  print *, 'Iteraciones (variacional) = ', it
  print *, 'Alpha optimo             = ', alpha_opt
  print *, 'Beta optimo              = ', beta_opt
  print *, 'E estado base [eV]       = ', e0 * hartree_to_ev
  print *, 'E ionizacion [eV]        = ', ei * hartree_to_ev

  call build_potential_table(z, 2.535930d0, npt, r, vc, veff)
  call solve_l_spectrum(z, 2.535930d0, 0.0d0, npt, -10.0d0, 0.0d0, 2500, e_s, nfound_s)
  call solve_l_spectrum(z, 2.535930d0, 1.0d0, npt, -10.0d0, 0.0d0, 2500, e_p, nfound_p)
  call solve_l_spectrum(z, 2.535930d0, 2.0d0, npt, -10.0d0, 0.0d0, 2500, e_d, nfound_d)

  if (nfound_s < 4) stop 'No se encontraron suficientes estados s para 2s/3s/4s'
  if (nfound_p < 2) stop 'No se encontraron suficientes estados p para 2p/3p'
  if (nfound_d < 1) stop 'No se encontro estado d para 3d'

  call wave_for_energy(z, 2.535930d0, 0.0d0, npt, e_s(2), rs2)
  call wave_for_energy(z, 2.535930d0, 0.0d0, npt, e_s(3), rs3)
  call wave_for_energy(z, 2.535930d0, 0.0d0, npt, e_s(4), rs4)
  call wave_for_energy(z, 2.535930d0, 1.0d0, npt, e_p(1), rp2)
  call wave_for_energy(z, 2.535930d0, 1.0d0, npt, e_p(2), rp3)
  call wave_for_energy(z, 2.535930d0, 2.0d0, npt, e_d(1), rd3)

  open(10, file='potential_litio.dat', status='replace', action='write')
  write(10,'(a)') '# r(a.u.)  Vcoulomb(Ha)  Veff(Ha)'
  do i = 1, npt
    write(10,'(3(es22.14,1x))') r(i), vc(i), veff(i)
  end do
  close(10)

  open(11, file='funciones_radiales_litio.dat', status='replace', action='write')
  write(11,'(a)') '# r(a.u.)  R_2s  R_3s  R_4s  R_2p  R_3p  R_3d'
  do i = 1, npt
    write(11,'(7(es22.14,1x))') r(i), rs2(i), rs3(i), rs4(i), rp2(i), rp3(i), rd3(i)
  end do
  close(11)

  open(12, file='energias_litio_nr.dat', status='replace', action='write')
  write(12,'(a)') '# estado  l  E(Ha)  E(eV)'
  write(12,'(a,1x,i1,1x,es22.14,1x,es22.14)') '2s', 0, e_s(2), e_s(2)*hartree_to_ev
  write(12,'(a,1x,i1,1x,es22.14,1x,es22.14)') '3s', 0, e_s(3), e_s(3)*hartree_to_ev
  write(12,'(a,1x,i1,1x,es22.14,1x,es22.14)') '4s', 0, e_s(4), e_s(4)*hartree_to_ev
  write(12,'(a,1x,i1,1x,es22.14,1x,es22.14)') '2p', 1, e_p(1), e_p(1)*hartree_to_ev
  write(12,'(a,1x,i1,1x,es22.14,1x,es22.14)') '3p', 1, e_p(2), e_p(2)*hartree_to_ev
  write(12,'(a,1x,i1,1x,es22.14,1x,es22.14)') '3d', 2, e_d(1), e_d(1)*hartree_to_ev
  close(12)

  print *, 'Archivos generados:'
  print *, '  potential_litio.dat'
  print *, '  funciones_radiales_litio.dat'
  print *, '  energias_litio_nr.dat'

contains

  subroutine optimize_variational(zin, a0, b0, aopt, bopt, eopt, iterations)
    implicit none
    real(dp), intent(in) :: zin, a0, b0
    real(dp), intent(out) :: aopt, bopt, eopt
    integer, intent(out) :: iterations

    real(dp) :: p(2), xi(2,2), fret

    z = zin
    p = [a0, b0]
    xi = 0.0d0
    xi(1,1) = 1.0d0
    xi(2,2) = 1.0d0

    call powell(p, xi, 2, 1.0d-8, iterations, fret, func_vec)

    aopt = p(1)
    bopt = p(2)
    eopt = fret
  end subroutine optimize_variational

  real(dp) function func_vec(x, n)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp) :: a, b, y

    a = x(1)
    b = x(2)
    y = 2.0d0 * a * b / ((2.0d0 * a + b)**5)

    func_vec = a*a - 2.0d0*z*a + (5.0d0/8.0d0)*a + (1.0d0/8.0d0)*b*b - (z*b/4.0d0) + &
               y*(8.0d0*a**4 + 20.0d0*a**3*b + 12.0d0*a*a*b*b + 10.0d0*a*b**3 + b**4)
  end function func_vec

  real(dp) function ionization_energy_li(alpha, beta, zloc)
    implicit none
    real(dp), intent(in) :: alpha, beta, zloc
    real(dp) :: y

    y = 2.0d0 * alpha * beta / (2.0d0 * alpha + beta)**5
    ionization_energy_li = beta**2 / 8.0d0 - zloc * beta / 4.0d0 + &
                           y * (8.0d0 * alpha**4 + 20.0d0 * alpha**3 * beta + 12.0d0 * alpha**2 * beta**2 + &
                                10.0d0 * alpha * beta**3 + beta**4)
  end function ionization_energy_li

  subroutine build_potential_table(zloc, aloc, n, rr, vcoul, veffo)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: zloc, aloc
    real(dp), intent(out) :: rr(n), vcoul(n), veffo(n)
    integer :: k
    real(dp) :: rmin, rmax, h

    rmin = 1.0d-4
    rmax = 120.0d0
    h = (rmax - rmin) / real(n - 1, dp)

    do k = 1, n
      rr(k) = rmin + real(k - 1, dp) * h
      vcoul(k) = -(zloc - 2.0d0) / rr(k)
      veffo(k) = -zloc / rr(k) + 2.0d0 / rr(k) * (1.0d0 - (1.0d0 + aloc*rr(k)) * exp(-2.0d0*aloc*rr(k)))
    end do
  end subroutine build_potential_table

  subroutine solve_l_spectrum(zloc, aloc, lloc, n, x1, x2, nscan, energies, nfound)
    implicit none
    integer, intent(in) :: n, nscan
    real(dp), intent(in) :: zloc, aloc, lloc, x1, x2
    real(dp), intent(out) :: energies(nbmax)
    integer, intent(out) :: nfound

    integer :: i, nb
    real(dp) :: xb1(nbmax), xb2(nbmax), tol, root

    call set_spec_params(zloc, aloc, lloc, n)

    nb = nbmax
    call zbrak(determ_l, x1, x2, nscan, xb1, xb2, nb)

    nfound = nb
    energies = 0.0d0

    do i = 1, nb
      tol = 1.0d-12 * abs(0.5d0 * (xb1(i) + xb2(i)))
      if (tol == 0.0d0) tol = 1.0d-12
      root = zbrent(determ_l, xb1(i), xb2(i), tol)
      energies(i) = root
    end do
  end subroutine solve_l_spectrum

  subroutine wave_for_energy(zloc, aloc, lloc, n, energy, wave)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: zloc, aloc, lloc, energy
    real(dp), intent(out) :: wave(n)

    integer :: k, kk
    real(dp) :: h, c(n), diag(n), superd(n), subd(n), vk(n), vkm1(n), nrm

    call set_spec_params(zloc, aloc, lloc, n)
    h = h_spec

    do k = 1, n
      c(k) = pot_l(h * real(k, dp))
      diag(k) = 2.0d0 + 2.0d0*h*h*c(k) - 2.0d0*h*h*energy
    end do

    do k = 1, n - 1
      superd(k) = -1.0d0
    end do
    superd(n) = 0.0d0

    subd(1) = 0.0d0
    do k = 2, n
      subd(k) = -1.0d0
    end do

    vk = 0.0d0
    vk(1) = 1.0d0

    do kk = 1, 10
      call tridag(subd, diag, superd, vk, vkm1, n)
      nrm = sqrt(sum(vkm1*vkm1) * h)
      if (nrm > 0.0d0) then
        vk = vkm1 / nrm
      else
        vk = vkm1
      end if
    end do

    wave = vk
  end subroutine wave_for_energy

  subroutine set_spec_params(zloc, aloc, lloc, n)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: zloc, aloc, lloc
    real(dp) :: d

    d = 120.0d0
    nn_spec = n
    h_spec = d / real(nn_spec, dp)
    ele_spec = lloc
    z_spec = zloc
    alpha_spec = aloc
  end subroutine set_spec_params

  real(dp) function pot_l(x)
    implicit none
    real(dp), intent(in) :: x

    pot_l = -z_spec / x + 2.0d0 / x * (1.0d0 - (1.0d0 + alpha_spec*x) * exp(-2.0d0*alpha_spec*x))
    pot_l = pot_l + ele_spec * (ele_spec + 1.0d0) / (2.0d0 * x * x)
  end function pot_l

  real(dp) function determ_l(x)
    implicit none
    real(dp), intent(in) :: x

    integer :: i
    real(dp) :: h2, p0, p1, p2

    h2 = h_spec * h_spec
    p0 = 0.0d0
    p1 = 2.0d0 + 2.0d0*h2*pot_l(h_spec) - 2.0d0*x*h2

    do i = 2, nn_spec
      p2 = (2.0d0 + 2.0d0*h2*pot_l(h_spec*real(i,dp)) - 2.0d0*h2*x) * p1 - p0
      p0 = p1
      p1 = p2
    end do

    determ_l = p2
  end function determ_l

end program litio_nr_pipeline
