module nr_powell_compat
  implicit none
  private
  public :: powell

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: itmax = 200
  real(dp), parameter :: tiny = 1.0d-25

contains

  subroutine powell(p, xi, n, ftol, iter, fret, func)
    implicit none

    integer, intent(in) :: n
    real(dp), intent(inout) :: p(n), xi(n, n)
    real(dp), intent(in) :: ftol
    integer, intent(out) :: iter
    real(dp), intent(out) :: fret

    interface
      function func(x, n) result(f)
        import dp
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        real(dp) :: f
      end function func
    end interface

    integer :: i, j, ibig
    real(dp) :: del, fp, fptt, t
    real(dp), allocatable :: pt(:), ptt(:), xit(:)

    allocate(pt(n), ptt(n), xit(n))

    fret = func(p, n)
    pt = p

    do iter = 1, itmax
      fp = fret
      ibig = 1
      del = 0.0d0

      do i = 1, n
        xit = xi(:, i)
        fptt = fret
        call linmin(p, xit, n, fret, func)
        if (fptt - fret > del) then
          del = fptt - fret
          ibig = i
        end if
      end do

      if (2.0d0 * (fp - fret) <= ftol * (abs(fp) + abs(fret)) + tiny) then
        deallocate(pt, ptt, xit)
        return
      end if

      ptt = 2.0d0 * p - pt
      xit = p - pt
      pt = p

      fptt = func(ptt, n)

      if (fptt < fp) then
        t = 2.0d0 * (fp - 2.0d0 * fret + fptt) * (fp - fret - del)**2 - del * (fp - fptt)**2
        if (t < 0.0d0) then
          call linmin(p, xit, n, fret, func)
          xi(:, ibig) = xi(:, n)
          xi(:, n) = xit
        end if
      end if
    end do

    deallocate(pt, ptt, xit)
    stop 'powell exceeding maximum iterations.'
  end subroutine powell


  subroutine linmin(p, xi, n, fret, func)
    implicit none

    integer, intent(in) :: n
    real(dp), intent(inout) :: p(n), xi(n)
    real(dp), intent(out) :: fret

    interface
      function func(x, n) result(f)
        import dp
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        real(dp) :: f
      end function func
    end interface

    integer :: i
    real(dp), parameter :: gold = 1.618034d0
    real(dp) :: ax, bx, cx, fa, fb, fc
    real(dp) :: xmin, fmin
    real(dp) :: p0(n), dir(n)

    p0 = p
    dir = xi

    ax = 0.0d0
    bx = 1.0d0
    fa = f1d(ax)
    fb = f1d(bx)

    if (fb > fa) then
      call swap(ax, bx)
      call swap(fa, fb)
    end if

    cx = bx + gold * (bx - ax)
    fc = f1d(cx)

    do i = 1, 100
      if (fb <= fc) exit
      ax = bx
      fa = fb
      bx = cx
      fb = fc
      cx = bx + gold * (bx - ax)
      fc = f1d(cx)
    end do

    call brent_min(ax, bx, cx, 1.0d-10, xmin, fmin)

    p = p0 + xmin * dir
    xi = xmin * dir
    fret = fmin

  contains

    function f1d(alpha) result(f)
      real(dp), intent(in) :: alpha
      real(dp) :: f
      real(dp) :: xt(n)

      xt = p0 + alpha * dir
      f = func(xt, n)
    end function f1d

    subroutine brent_min(ax, bx, cx, tol, xmin, fmin)
      real(dp), intent(in) :: ax, bx, cx, tol
      real(dp), intent(out) :: xmin, fmin

      integer, parameter :: itmax_b = 200
      real(dp), parameter :: cgold = 0.3819660d0
      real(dp), parameter :: zeps = 1.0d-12
      integer :: iter
      real(dp) :: a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm, e

      a = min(ax, cx)
      b = max(ax, cx)
      x = bx
      w = bx
      v = bx
      fw = f1d(w)
      fv = fw
      fx = fw
      d = 0.0d0
      e = 0.0d0

      do iter = 1, itmax_b
        xm = 0.5d0 * (a + b)
        tol1 = tol * abs(x) + zeps
        tol2 = 2.0d0 * tol1

        if (abs(x - xm) <= (tol2 - 0.5d0 * (b - a))) then
          xmin = x
          fmin = fx
          return
        end if

        if (abs(e) > tol1) then
          r = (x - w) * (fx - fv)
          q = (x - v) * (fx - fw)
          p = (x - v) * q - (x - w) * r
          q = 2.0d0 * (q - r)
          if (q > 0.0d0) p = -p
          q = abs(q)
          etemp = e
          e = d

          if (abs(p) >= abs(0.5d0 * q * etemp) .or. p <= q * (a - x) .or. p >= q * (b - x)) then
            if (x >= xm) then
              e = a - x
            else
              e = b - x
            end if
            d = cgold * e
          else
            d = p / q
            u = x + d
            if (u - a < tol2 .or. b - u < tol2) d = sign(tol1, xm - x)
          end if
        else
          if (x >= xm) then
            e = a - x
          else
            e = b - x
          end if
          d = cgold * e
        end if

        if (abs(d) >= tol1) then
          u = x + d
        else
          u = x + sign(tol1, d)
        end if

        fu = f1d(u)

        if (fu <= fx) then
          if (u >= x) then
            a = x
          else
            b = x
          end if
          v = w
          fv = fw
          w = x
          fw = fx
          x = u
          fx = fu
        else
          if (u < x) then
            a = u
          else
            b = u
          end if
          if (fu <= fw .or. w == x) then
            v = w
            fv = fw
            w = u
            fw = fu
          else if (fu <= fv .or. v == x .or. v == w) then
            v = u
            fv = fu
          end if
        end if
      end do

      xmin = x
      fmin = fx
    end subroutine brent_min

    subroutine swap(a, b)
      real(dp), intent(inout) :: a, b
      real(dp) :: tmp
      tmp = a
      a = b
      b = tmp
    end subroutine swap

  end subroutine linmin

end module nr_powell_compat
