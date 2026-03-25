program litio_numerov
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: n_states = 6
  integer, parameter :: n_grid = 12000

  real(dp), parameter :: alpha = 2.535930_dp
  real(dp), parameter :: r_min = 1.0d-4
  real(dp), parameter :: r_max = 35.0_dp
  real(dp), parameter :: e_max = -1.0d-4
  real(dp), parameter :: e_tol = 1.0d-10
  real(dp), parameter :: ev_per_au = 27.211386_dp

  character(len=8), dimension(n_states) :: labels
  integer, dimension(n_states) :: nvals, lvals, target_nodes

  real(dp), dimension(n_states) :: energies_au, energies_ev
  real(dp), dimension(n_states) :: e_lo_guess, e_hi_guess
  real(dp), dimension(n_grid) :: r_grid, r_tmp
  real(dp), dimension(n_grid, n_states) :: waves

  integer :: i
  logical :: ok

  labels = [character(len=8) :: '2s', '3s', '4s', '2p', '3p', '3d']
  nvals = [2, 3, 4, 2, 3, 3]
  lvals = [0, 0, 0, 1, 1, 2]

  do i = 1, n_states
    target_nodes(i) = nvals(i) - lvals(i) - 1
  end do

  ! Ventanas fisicas de busqueda para estados de valencia del Li.
  e_lo_guess = [ -0.22_dp, -0.09_dp, -0.05_dp, -0.16_dp, -0.07_dp, -0.07_dp ]
  e_hi_guess = [ -0.12_dp, -0.05_dp, -0.02_dp, -0.10_dp, -0.045_dp, -0.045_dp ]

  call build_grid(r_grid)

  write(*, '(a)') 'Resolviendo espectro del Litio con Numerov (Fortran)...'
  do i = 1, n_states
    call find_state_energy_numerov(lvals(i), target_nodes(i), e_lo_guess(i), e_hi_guess(i), e_tol, energies_au(i), ok)
    if (.not. ok) then
      write(*, '(a,a)') 'No se encontro estado para: ', trim(labels(i))
      stop 1
    end if

    energies_ev(i) = ev_per_au * energies_au(i)
    call integrate_wave_numerov(energies_au(i), lvals(i), r_tmp, ok)
    if (.not. ok) then
      write(*, '(a,a)') 'Fallo integracion de onda para: ', trim(labels(i))
      stop 1
    end if

    call normalize_wave(r_grid, r_tmp)
    waves(:, i) = r_tmp

    write(*, '(a,1x,a,2x,a,f12.6,2x,a,f12.6)') 'Estado', trim(labels(i)), 'E(a.u.)=', energies_au(i), 'E(eV)=', energies_ev(i)
  end do

  call write_energies('energias_numerov_fortran.dat', labels, lvals, target_nodes, energies_au, energies_ev)
  call write_waves('funciones_numerov_fortran.dat', labels, r_grid, waves)

  write(*, '(a)') 'Archivos generados:'
  write(*, '(a)') '  energias_numerov_fortran.dat'
  write(*, '(a)') '  funciones_numerov_fortran.dat'

contains

  subroutine build_grid(r)
    real(dp), intent(out) :: r(n_grid)
    real(dp) :: h
    integer :: j

    h = (r_max - r_min) / real(n_grid - 1, dp)
    do j = 1, n_grid
      r(j) = r_min + real(j - 1, dp) * h
    end do
  end subroutine build_grid

  real(dp) function v_eff(r, l)
    real(dp), intent(in) :: r
    integer, intent(in) :: l
    real(dp) :: centrifugal

    centrifugal = real(l * (l + 1), dp) / (2.0_dp * r * r)
    v_eff = (-1.0_dp / r) + (2.0_dp / r) * (1.0_dp + alpha * r) * exp(-2.0_dp * alpha * r) + centrifugal
  end function v_eff

  real(dp) function g_fun(r, e, l)
    real(dp), intent(in) :: r, e
    integer, intent(in) :: l
    g_fun = 2.0_dp * (v_eff(r, l) - e)
  end function g_fun

  subroutine endpoint_nodes_numerov(e, l, endpoint, nodes)
    real(dp), intent(in) :: e
    integer, intent(in) :: l
    real(dp), intent(out) :: endpoint
    integer, intent(out) :: nodes

    real(dp) :: h, g0, g1, g2, c_im1, c_i, c_ip1
    real(dp) :: u_im1, u_i, u_ip1, r0, r1, r2
    integer :: i

    h = (r_max - r_min) / real(n_grid - 1, dp)

    r0 = r_min
    r1 = r_min + h
    u_im1 = r0 ** real(l + 1, dp)
    u_i = r1 ** real(l + 1, dp)
    nodes = 0

    do i = 2, n_grid - 1
      r2 = r_min + real(i, dp) * h

      g0 = g_fun(r0, e, l)
      g1 = g_fun(r1, e, l)
      g2 = g_fun(r2, e, l)

      c_im1 = 1.0_dp + (h * h / 12.0_dp) * g0
      c_i = 1.0_dp - (5.0_dp * h * h / 12.0_dp) * g1
      c_ip1 = 1.0_dp + (h * h / 12.0_dp) * g2

      if (abs(c_ip1) < 1.0d-16) then
        endpoint = huge(1.0_dp)
        return
      end if

      u_ip1 = (2.0_dp * c_i * u_i - c_im1 * u_im1) / c_ip1

      if (.not. ieee_is_finite(u_ip1)) then
        endpoint = huge(1.0_dp)
        return
      end if

      if (u_i * u_ip1 < 0.0_dp) nodes = nodes + 1

      if (abs(u_ip1) > 1.0d80) then
        u_im1 = u_im1 * 1.0d-80
        u_i = u_i * 1.0d-80
        u_ip1 = u_ip1 * 1.0d-80
      end if

      u_im1 = u_i
      u_i = u_ip1
      r0 = r1
      r1 = r2
    end do

    endpoint = u_i
  end subroutine endpoint_nodes_numerov

  subroutine integrate_wave_numerov(e, l, wave, ok)
    real(dp), intent(in) :: e
    integer, intent(in) :: l
    real(dp), intent(out) :: wave(n_grid)
    logical, intent(out) :: ok

    real(dp) :: h, g0, g1, g2, c_im1, c_i, c_ip1
    real(dp) :: u_im1, u_i, u_ip1, r0, r1, r2
    integer :: i

    ok = .false.
    h = (r_max - r_min) / real(n_grid - 1, dp)

    r0 = r_min
    r1 = r_min + h
    u_im1 = r0 ** real(l + 1, dp)
    u_i = r1 ** real(l + 1, dp)

    wave(1) = u_im1
    wave(2) = u_i

    do i = 2, n_grid - 1
      r2 = r_min + real(i, dp) * h

      g0 = g_fun(r0, e, l)
      g1 = g_fun(r1, e, l)
      g2 = g_fun(r2, e, l)

      c_im1 = 1.0_dp + (h * h / 12.0_dp) * g0
      c_i = 1.0_dp - (5.0_dp * h * h / 12.0_dp) * g1
      c_ip1 = 1.0_dp + (h * h / 12.0_dp) * g2

      if (abs(c_ip1) < 1.0d-16) return

      u_ip1 = (2.0_dp * c_i * u_i - c_im1 * u_im1) / c_ip1
      if (.not. ieee_is_finite(u_ip1)) return

      if (abs(u_ip1) > 1.0d80) then
        u_im1 = u_im1 * 1.0d-80
        u_i = u_i * 1.0d-80
        u_ip1 = u_ip1 * 1.0d-80
      end if

      wave(i + 1) = u_ip1

      u_im1 = u_i
      u_i = u_ip1
      r0 = r1
      r1 = r2
    end do

    ok = .true.
  end subroutine integrate_wave_numerov

  subroutine refine_bisection_numerov(l, target_nodes, ea, eb, tol, e_star, ok)
    integer, intent(in) :: l, target_nodes
    real(dp), intent(in) :: ea, eb, tol
    real(dp), intent(out) :: e_star
    logical, intent(out) :: ok

    integer, parameter :: itmax = 220
    integer :: iter, na, nb, nm
    real(dp) :: a, b, m, fa, fb, fm

    ok = .false.
    a = ea
    b = eb

    call endpoint_nodes_numerov(a, l, fa, na)
    call endpoint_nodes_numerov(b, l, fb, nb)
    if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fb)) return
    if (na /= target_nodes .or. nb /= target_nodes) return
    if (fa * fb > 0.0_dp) return

    do iter = 1, itmax
      m = 0.5_dp * (a + b)
      call endpoint_nodes_numerov(m, l, fm, nm)
      if (.not. ieee_is_finite(fm)) return

      if (abs(b - a) < tol .or. abs(fm) < 1.0d-12) then
        e_star = m
        ok = .true.
        return
      end if

      if (nm == target_nodes .and. fa * fm <= 0.0_dp) then
        b = m
        fb = fm
      else if (nm == target_nodes) then
        a = m
        fa = fm
      else
        ! Si cambia nodos, tomar subintervalo que mantenga nodos objetivo.
        if (nm > target_nodes) then
          b = m
          call endpoint_nodes_numerov(b, l, fb, nb)
        else
          a = m
          call endpoint_nodes_numerov(a, l, fa, na)
        end if
        if (na /= target_nodes .or. nb /= target_nodes) return
      end if
    end do

    e_star = 0.5_dp * (a + b)
    ok = .true.
  end subroutine refine_bisection_numerov

  subroutine find_state_energy_numerov(l, target_nodes, e_lo, e_hi, tol, e_star, ok)
    integer, intent(in) :: l, target_nodes
    real(dp), intent(in) :: e_lo, e_hi, tol
    real(dp), intent(out) :: e_star
    logical, intent(out) :: ok

    integer, parameter :: n_scan = 20000
    integer :: i, n_cur
    real(dp) :: de, e_cur, f_cur, score, best_score
    real(dp) :: pen_nodes

    ok = .false.
    e_star = 0.5_dp * (e_lo + e_hi)
    best_score = huge(1.0_dp)

    de = (e_hi - e_lo) / real(n_scan - 1, dp)

    do i = 1, n_scan
      e_cur = e_lo + real(i - 1, dp) * de
      call endpoint_nodes_numerov(e_cur, l, f_cur, n_cur)
      if (.not. ieee_is_finite(f_cur)) cycle

      pen_nodes = 1.0_dp + 10.0_dp * abs(real(n_cur - target_nodes, dp))
      score = abs(f_cur) * pen_nodes

      if (score < best_score) then
        best_score = score
        e_star = e_cur
        ok = .true.
      end if
    end do

    if (.not. ok) return

    ! Refinamiento local simple alrededor del mejor valor encontrado.
    do i = 1, 20
      de = max(tol, 0.5_dp * de)
      call endpoint_nodes_numerov(max(e_lo, e_star - de), l, f_cur, n_cur)
      if (ieee_is_finite(f_cur)) then
        pen_nodes = 1.0_dp + 10.0_dp * abs(real(n_cur - target_nodes, dp))
        score = abs(f_cur) * pen_nodes
        if (score < best_score) then
          best_score = score
          e_star = max(e_lo, e_star - de)
        end if
      end if

      call endpoint_nodes_numerov(min(e_hi, e_star + de), l, f_cur, n_cur)
      if (ieee_is_finite(f_cur)) then
        pen_nodes = 1.0_dp + 10.0_dp * abs(real(n_cur - target_nodes, dp))
        score = abs(f_cur) * pen_nodes
        if (score < best_score) then
          best_score = score
          e_star = min(e_hi, e_star + de)
        end if
      end if
    end do

    ok = .true.
  end subroutine find_state_energy_numerov

  subroutine normalize_wave(r, y)
    real(dp), intent(in) :: r(n_grid)
    real(dp), intent(inout) :: y(n_grid)

    real(dp) :: norm
    integer :: j

    norm = 0.0_dp
    do j = 1, n_grid - 1
      norm = norm + 0.5_dp * (y(j) * y(j) + y(j + 1) * y(j + 1)) * (r(j + 1) - r(j))
    end do

    if (norm > 0.0_dp .and. ieee_is_finite(norm)) then
      y = y / sqrt(norm)
    end if
  end subroutine normalize_wave

  subroutine write_energies(fname, lbl, ll, nnodes, e_au, e_ev)
    character(len=*), intent(in) :: fname
    character(len=*), dimension(n_states), intent(in) :: lbl
    integer, dimension(n_states), intent(in) :: ll, nnodes
    real(dp), dimension(n_states), intent(in) :: e_au, e_ev

    integer :: unit, i

    unit = 21
    open(unit=unit, file=fname, status='replace', action='write')
    write(unit, '(a)') '# label l nodes E_au E_eV'
    do i = 1, n_states
      write(unit, '(a,1x,i2,1x,i2,1x,f20.12,1x,f20.12)') trim(lbl(i)), ll(i), nnodes(i), e_au(i), e_ev(i)
    end do
    close(unit)
  end subroutine write_energies

  subroutine write_waves(fname, lbl, r, w)
    character(len=*), intent(in) :: fname
    character(len=*), dimension(n_states), intent(in) :: lbl
    real(dp), dimension(n_grid), intent(in) :: r
    real(dp), dimension(n_grid, n_states), intent(in) :: w

    integer :: unit, i

    unit = 22
    open(unit=unit, file=fname, status='replace', action='write')
    write(unit, '(a)') '# r R_2s R_3s R_4s R_2p R_3p R_3d'
    do i = 1, n_grid
      write(unit, '(f16.8,1x,6(1x,es20.10))') r(i), w(i, 1), w(i, 2), w(i, 3), w(i, 4), w(i, 5), w(i, 6)
    end do
    close(unit)
  end subroutine write_waves

end program litio_numerov
