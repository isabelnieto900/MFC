module nr_compat_litio
  implicit none
  private
  public :: zbrak, zbrent, tridag

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: itmax = 100
  real(dp), parameter :: eps = 3.0d-8

contains

  subroutine zbrak(fx, x1, x2, n, xb1, xb2, nb)
    implicit none

    integer, intent(in) :: n
    real(dp), intent(in) :: x1, x2
    real(dp), intent(out) :: xb1(:), xb2(:)
    integer, intent(inout) :: nb

    integer :: nbb, i, nbmax
    real(dp) :: x, fp, fc, dx
    real(dp), external :: fx

    nbb = 0
    nbmax = nb
    dx = (x2 - x1) / real(n, dp)
    x = x1
    fp = fx(x)

    do i = 1, n
      x = x + dx
      fc = fx(x)

      if (fc * fp <= 0.0d0) then
        nbb = nbb + 1
        if (nbb > size(xb1) .or. nbb > size(xb2)) exit

        xb1(nbb) = x - dx
        xb2(nbb) = x

        if (nbb == nbmax) then
          nb = nbb
          return
        end if
      end if

      fp = fc
    end do

    nb = nbb
  end subroutine zbrak


  real(dp) function zbrent(func, x1, x2, tol)
    implicit none

    real(dp), intent(in) :: x1, x2, tol
    real(dp), external :: func

    integer :: iter
    real(dp) :: a, b, c, d, e, min1, min2
    real(dp) :: fa, fb, fc, p, q, r, s, tol1, xm

    a = x1
    b = x2
    c = x2
    d = 0.0d0
    e = 0.0d0

    fa = func(a)
    fb = func(b)

    if ((fa > 0.0d0 .and. fb > 0.0d0) .or. (fa < 0.0d0 .and. fb < 0.0d0)) then
      stop 'Root must be bracketed in zbrent'
    end if

    fc = fb

    do iter = 1, itmax
      if ((fb > 0.0d0 .and. fc > 0.0d0) .or. (fb < 0.0d0 .and. fc < 0.0d0)) then
        c = a
        fc = fa
        d = b - a
        e = d
      end if

      if (abs(fc) < abs(fb)) then
        a = b
        b = c
        c = a
        fa = fb
        fb = fc
        fc = fa
      end if

      tol1 = 2.0d0 * eps * abs(b) + 0.5d0 * tol
      xm = 0.5d0 * (c - b)

      if (abs(xm) <= tol1 .or. fb == 0.0d0) then
        zbrent = b
        return
      end if

      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s = fb / fa

        if (a == c) then
          p = 2.0d0 * xm * s
          q = 1.0d0 - s
        else
          q = fa / fc
          r = fb / fc
          p = s * (2.0d0 * xm * q * (q - r) - (b - a) * (r - 1.0d0))
          q = (q - 1.0d0) * (r - 1.0d0) * (s - 1.0d0)
        end if

        if (p > 0.0d0) q = -q
        p = abs(p)

        min1 = 3.0d0 * xm * q - abs(tol1 * q)
        min2 = abs(e * q)

        if (2.0d0 * p < min(min1, min2)) then
          e = d
          d = p / q
        else
          d = xm
          e = d
        end if
      else
        d = xm
        e = d
      end if

      a = b
      fa = fb

      if (abs(d) > tol1) then
        b = b + d
      else
        b = b + sign(tol1, xm)
      end if

      fb = func(b)
    end do

    stop 'Maximum number of iterations exceeded in zbrent'
  end function zbrent


  subroutine tridag(a, b, c, r, u, n)
    implicit none

    integer, intent(in) :: n
    real(dp), intent(in) :: a(n), b(n), c(n), r(n)
    real(dp), intent(out) :: u(n)

    integer :: j
    real(dp) :: bet
    real(dp), allocatable :: gam(:)

    allocate(gam(n))

    if (b(1) == 0.0d0) then
      deallocate(gam)
      stop 'Error 1 in tridag'
    end if

    bet = b(1)
    u(1) = r(1) / bet

    do j = 2, n
      gam(j) = c(j - 1) / bet
      bet = b(j) - a(j) * gam(j)

      if (bet == 0.0d0) then
        deallocate(gam)
        stop 'Error 2 in tridag'
      end if

      u(j) = (r(j) - a(j) * u(j - 1)) / bet
    end do

    do j = n - 1, 1, -1
      u(j) = u(j) - gam(j + 1) * u(j + 1)
    end do

    deallocate(gam)
  end subroutine tridag

end module nr_compat_litio
