program litio_shooting
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: n_states = 6
  integer, parameter :: n_grid = 6000

  real(dp), parameter :: alpha = 2.535930_dp
  real(dp), parameter :: r_min = 1.0d-4
  real(dp), parameter :: r_max = 35.0_dp
  real(dp), parameter :: e_min = -10.0_dp
  real(dp), parameter :: e_max = -1.0d-4
  real(dp), parameter :: e_tol = 1.0d-10
  real(dp), parameter :: ev_per_au = 27.211386_dp

  character(len=8), dimension(n_states) :: labels
  integer, dimension(n_states) :: nvals, lvals, target_nodes

  real(dp), dimension(n_states) :: energies_au, energies_ev
  real(dp), dimension(n_grid) :: r_grid, r_tmp
  real(dp), dimension(n_grid, n_states) :: waves

  integer :: i, k, rank_l
  logical :: ok

  labels = [character(len=8) :: '2s', '3s', '4s', '2p', '3p', '3d']
  nvals = [2, 3, 4, 2, 3, 3]
  lvals = [0, 0, 0, 1, 1, 2]

  do i = 1, n_states
    target_nodes(i) = nvals(i) - lvals(i) - 1
  end do

  call build_grid(r_grid)

  write(*, '(a)') 'Resolviendo espectro del Litio con metodo de shooting...'
  do i = 1, n_states
    rank_l = 1
    do k = 1, i - 1
      if (lvals(k) == lvals(i)) rank_l = rank_l + 1
    end do

    call find_state_energy_rank(lvals(i), rank_l, e_min, e_max, e_tol, energies_au(i), ok)
    if (.not. ok) then
      write(*, '(a,a)') 'No se encontro estado para: ', trim(labels(i))
      stop 1
    end if
    energies_ev(i) = ev_per_au * energies_au(i)

    call integrate_full_state(energies_au(i), lvals(i), r_tmp)
    call normalize_wave(r_grid, r_tmp)
    waves(:, i) = r_tmp

    write(*, '(a,1x,a,2x,a,f12.6,2x,a,f12.6)') 'Estado', trim(labels(i)), 'E(a.u.)=', energies_au(i), 'E(eV)=', energies_ev(i)
  end do

  call write_energies('energias_shooting_fortran.dat', labels, lvals, target_nodes, energies_au, energies_ev)
  call write_waves('funciones_shooting_fortran.dat', labels, r_grid, waves)

  write(*, '(a)') 'Archivos generados:'
  write(*, '(a)') '  energias_shooting_fortran.dat'
  write(*, '(a)') '  funciones_shooting_fortran.dat'

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

  subroutine rhs_system(r, y1, y2, e, l, dy1, dy2)
    real(dp), intent(in) :: r, y1, y2, e
    integer, intent(in) :: l
    real(dp), intent(out) :: dy1, dy2

    dy1 = y2
    dy2 = 2.0_dp * (v_eff(r, l) - e) * y1
  end subroutine rhs_system

  subroutine rk4_step(r, h, y1, y2, e, l)
    real(dp), intent(in) :: r, h, e
    integer, intent(in) :: l
    real(dp), intent(inout) :: y1, y2

    real(dp) :: k1y1, k1y2, k2y1, k2y2, k3y1, k3y2, k4y1, k4y2
    real(dp) :: yt1, yt2

    call rhs_system(r, y1, y2, e, l, k1y1, k1y2)

    yt1 = y1 + 0.5_dp * h * k1y1
    yt2 = y2 + 0.5_dp * h * k1y2
    call rhs_system(r + 0.5_dp * h, yt1, yt2, e, l, k2y1, k2y2)

    yt1 = y1 + 0.5_dp * h * k2y1
    yt2 = y2 + 0.5_dp * h * k2y2
    call rhs_system(r + 0.5_dp * h, yt1, yt2, e, l, k3y1, k3y2)

    yt1 = y1 + h * k3y1
    yt2 = y2 + h * k3y2
    call rhs_system(r + h, yt1, yt2, e, l, k4y1, k4y2)

    y1 = y1 + (h / 6.0_dp) * (k1y1 + 2.0_dp * k2y1 + 2.0_dp * k3y1 + k4y1)
    y2 = y2 + (h / 6.0_dp) * (k1y2 + 2.0_dp * k2y2 + 2.0_dp * k3y2 + k4y2)
  end subroutine rk4_step

  subroutine endpoint_and_nodes(e, l, endpoint, nodes)
    real(dp), intent(in) :: e
    integer, intent(in) :: l
    real(dp), intent(out) :: endpoint
    integer, intent(out) :: nodes

    real(dp) :: y1, y2, h, r, y1_prev
    integer :: j

    h = (r_max - r_min) / real(n_grid - 1, dp)

    y1 = r_min ** real(l + 1, dp)
    y2 = real(l + 1, dp) * r_min ** real(l, dp)
    nodes = 0

    do j = 1, n_grid - 1
      r = r_min + real(j - 1, dp) * h
      y1_prev = y1

      call rk4_step(r, h, y1, y2, e, l)

      if (y1_prev * y1 < 0.0_dp) nodes = nodes + 1

      if (.not. ieee_is_finite(y1) .or. .not. ieee_is_finite(y2)) then
        endpoint = huge(1.0_dp)
        return
      end if

      if (abs(y1) > 1.0d60 .or. abs(y2) > 1.0d60) then
        y1 = y1 * 1.0d-60
        y2 = y2 * 1.0d-60
      end if
    end do

    endpoint = y1
  end subroutine endpoint_and_nodes

  subroutine find_state_energy_rank(l, rank_target, e_lo, e_hi, tol, e_star, ok)
    integer, intent(in) :: l, rank_target
    real(dp), intent(in) :: e_lo, e_hi, tol
    real(dp), intent(out) :: e_star
    logical, intent(out) :: ok

    integer, parameter :: n_scan = 12000
    integer :: i, nodes_mid
    real(dp) :: de, ea, eb, em, fa, fb, fm, last_e, last_f, cur_e, cur_f
    real(dp) :: last_root
    integer :: iter, max_iter
    integer :: roots_found

    ok = .false.
    e_star = 0.0_dp
    max_iter = 180
    last_root = huge(1.0_dp)
    roots_found = 0

    de = (e_hi - e_lo) / real(n_scan - 1, dp)
    last_e = e_lo
    call endpoint_and_nodes(last_e, l, last_f, nodes_mid)

    do i = 2, n_scan
      cur_e = e_lo + real(i - 1, dp) * de
      call endpoint_and_nodes(cur_e, l, cur_f, nodes_mid)

      if (ieee_is_finite(last_f) .and. ieee_is_finite(cur_f) .and. last_f * cur_f <= 0.0_dp) then
        ea = last_e
        eb = cur_e
        fa = last_f
        fb = cur_f

        do iter = 1, max_iter
          em = 0.5_dp * (ea + eb)
          call endpoint_and_nodes(em, l, fm, nodes_mid)

          if (abs(fm) < 1.0d-12 .or. abs(eb - ea) < tol) exit

          if (fa * fm <= 0.0_dp) then
            eb = em
            fb = fm
          else
            ea = em
            fa = fm
          end if
        end do

        e_star = 0.5_dp * (ea + eb)
        call endpoint_and_nodes(e_star, l, fm, nodes_mid)

        if (abs(e_star - last_root) > 1.0d-7) then
          roots_found = roots_found + 1
          last_root = e_star
        end if

        if (roots_found == rank_target) then
          ok = .true.
          return
        end if
      end if

      last_e = cur_e
      last_f = cur_f
    end do
  end subroutine find_state_energy_rank

  subroutine integrate_full_state(e, l, r_sol)
    real(dp), intent(in) :: e
    integer, intent(in) :: l
    real(dp), intent(out) :: r_sol(n_grid)

    real(dp) :: y1, y2, h, r
    integer :: j

    h = (r_max - r_min) / real(n_grid - 1, dp)

    y1 = r_min ** real(l + 1, dp)
    y2 = real(l + 1, dp) * r_min ** real(l, dp)
    r_sol(1) = y1

    do j = 2, n_grid
      r = r_min + real(j - 2, dp) * h
      call rk4_step(r, h, y1, y2, e, l)
      r_sol(j) = y1
    end do
  end subroutine integrate_full_state

  subroutine normalize_wave(r, y)
    real(dp), intent(in) :: r(n_grid)
    real(dp), intent(inout) :: y(n_grid)

    real(dp) :: norm
    integer :: j

    norm = 0.0_dp
    do j = 1, n_grid - 1
      norm = norm + 0.5_dp * (y(j) * y(j) + y(j + 1) * y(j + 1)) * (r(j + 1) - r(j))
    end do

    if (norm > 0.0_dp) then
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

end program litio_shooting
