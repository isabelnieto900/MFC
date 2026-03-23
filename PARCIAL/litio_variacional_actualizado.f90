module li_variational_model
  implicit none
  private
  public :: energy_li, ionization_energy_li, optimize_alpha_beta

  integer, parameter :: dp = kind(1.0d0)

contains

  pure function energy_li(alpha, beta, z) result(e)
    real(dp), intent(in) :: alpha, beta, z
    real(dp) :: e, y

    y = 2.0_dp * alpha * beta / (2.0_dp * alpha + beta)**5
    e = alpha**2 - 2.0_dp * z * alpha + 5.0_dp * alpha / 8.0_dp + beta**2 / 8.0_dp - z * beta / 4.0_dp + &
        y * (8.0_dp * alpha**4 + 20.0_dp * alpha**3 * beta + 12.0_dp * alpha**2 * beta**2 + &
             10.0_dp * alpha * beta**3 + beta**4)
  end function energy_li

  pure function ionization_energy_li(alpha, beta, z) result(ei)
    real(dp), intent(in) :: alpha, beta, z
    real(dp) :: ei, y

    y = 2.0_dp * alpha * beta / (2.0_dp * alpha + beta)**5
    ei = beta**2 / 8.0_dp - z * beta / 4.0_dp + &
         y * (8.0_dp * alpha**4 + 20.0_dp * alpha**3 * beta + 12.0_dp * alpha**2 * beta**2 + &
              10.0_dp * alpha * beta**3 + beta**4)
  end function ionization_energy_li

  subroutine optimize_alpha_beta(z, alpha0, beta0, alpha_opt, beta_opt, e_opt, iterations)
    real(dp), intent(in) :: z, alpha0, beta0
    real(dp), intent(out) :: alpha_opt, beta_opt, e_opt
    integer, intent(out) :: iterations

    real(dp) :: p(2), step(2), trial(2)
    real(dp) :: fbest, ftrial, tol
    integer :: i
    logical :: improved

    p = [max(alpha0, 1.0e-6_dp), max(beta0, 1.0e-6_dp)]
    step = [0.4_dp, 0.4_dp]
    tol = 1.0e-10_dp
    fbest = energy_li(p(1), p(2), z)

    do iterations = 1, 5000
      improved = .false.

      do i = 1, 2
        trial = p
        trial(i) = max(p(i) + step(i), 1.0e-8_dp)
        ftrial = energy_li(trial(1), trial(2), z)
        if (ftrial < fbest) then
          p = trial
          fbest = ftrial
          improved = .true.
        else
          trial = p
          trial(i) = max(p(i) - step(i), 1.0e-8_dp)
          ftrial = energy_li(trial(1), trial(2), z)
          if (ftrial < fbest) then
            p = trial
            fbest = ftrial
            improved = .true.
          end if
        end if
      end do

      if (.not. improved) step = 0.5_dp * step
      if (maxval(step) < tol) exit
    end do

    alpha_opt = p(1)
    beta_opt = p(2)
    e_opt = fbest
  end subroutine optimize_alpha_beta

end module li_variational_model

program litio_variacional_actualizado
  use li_variational_model
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: hartree_to_ev = 27.211386245988_dp

  real(dp) :: z
  real(dp) :: alpha0, beta0, alpha_opt, beta_opt, e0, ei
  integer :: it

  print *, 'Valor de Z (por ejemplo 3 para Li):'
  read(*, *) z

  alpha0 = 1.5_dp
  beta0 = 1.5_dp

  call optimize_alpha_beta(z, alpha0, beta0, alpha_opt, beta_opt, e0, it)
  ei = ionization_energy_li(alpha_opt, beta_opt, z)

  print *, 'Iteraciones            = ', it
  print *, 'Alpha optimo           = ', alpha_opt
  print *, 'Beta optimo            = ', beta_opt
  print *, 'Energia estado base [eV]= ', e0 * hartree_to_ev
  print *, 'Energia ionizacion [eV] = ', ei * hartree_to_ev

end program litio_variacional_actualizado
