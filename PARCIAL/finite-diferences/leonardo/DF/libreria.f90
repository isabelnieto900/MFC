module herramientas_fisica
    implicit none
    contains

    subroutine zbrak(fx, x1, x2, n, xb1, xb2, nbb)
        real(8), external :: fx
        real(8), intent(in) :: x1, x2
        integer, intent(in) :: n
        integer, intent(inout) :: nbb
        real(8), dimension(:), intent(out) :: xb1, xb2
        integer :: i, nbins
        real(8) :: x, fp, fc, dx

        nbins = nbb
        nbb = 0
        dx = (x2 - x1) / n
        x = x1
        fp = fx(x)
        do i = 1, n
            x = x + dx
            fc = fx(x)
            if (fc * fp < 0.0d0) then
                nbb = nbb + 1
                xb1(nbb) = x - dx
                xb2(nbb) = x
                if (nbb == nbins) return
            end if
            fp = fc
        end do
    end subroutine zbrak

    function zbrent(fx, x1, x2, tol)
        real(8), external :: fx
        real(8), intent(in) :: x1, x2, tol
        real(8) :: zbrent
        integer, parameter :: itmax = 100
        real(8), parameter :: eps = 3.0d-15
        integer :: iter
        real(8) :: a, b, c, d, e, min1, min2, fa, fb, fc, p, q, r, s, tol1, xm

        a = x1; b = x2; fa = fx(a); fb = fx(b)
        if ((fa > 0.0d0 .and. fb > 0.0d0) .or. (fa < 0.0d0 .and. fb < 0.0d0)) return
        c = a; fc = fa; d = b - a; e = d
        do iter = 1, itmax
            if ((fb > 0.0d0 .and. fc > 0.0d0) .or. (fb < 0.0d0 .and. fc < 0.0d0)) then
                c = a; fc = fa; d = b - a; e = d
            end if
            if (abs(fc) < abs(fb)) then
                a = b; b = c; c = a; fa = fb; fb = fc; fc = fa
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
                    q = fa / fc; r = fb / fc
                    p = s * (2.0d0 * xm * q * (q - r) - (b - a) * (r - 1.0d0))
                    q = (q - 1.0d0) * (r - 1.0d0) * (s - 1.0d0)
                end if
                if (p > 0.0d0) q = -q
                p = abs(p)
                min1 = 3.0d0 * xm * q - abs(tol1 * q)
                min2 = abs(e * q)
                if (2.0d0 * p < min(min1, min2)) then
                    e = d; d = p / q
                else
                    d = xm; e = d
                end if
            else
                d = xm; e = d
            end if
            a = b; fa = fb
            if (abs(d) > tol1) then
                b = b + d
            else
                b = b + sign(tol1, xm)
            end if
            fb = fx(b)
        end do
        zbrent = b
    end function zbrent

    subroutine tridag(a, b, c, r, u, n)
        integer, intent(in) :: n
        real(8), dimension(n), intent(in) :: a, b, c, r
        real(8), dimension(n), intent(out) :: u
        real(8), dimension(n) :: gam
        real(8) :: bet
        integer :: j
        bet = b(1)
        u(1) = r(1) / bet
        do j = 2, n
            gam(j) = c(j-1) / bet
            bet = b(j) - a(j) * gam(j)
            u(j) = (r(j) - a(j) * u(j-1)) / bet
        end do
        do j = n-1, 1, -1
            u(j) = u(j) - gam(j+1) * u(j+1)
        end do
    end subroutine tridag

    subroutine norma(n, v_in, v_out)
        integer, intent(in) :: n
        real(8), dimension(n), intent(in) :: v_in
        real(8), dimension(n), intent(out) :: v_out
        real(8) :: s
        integer :: i
        s = 0.0d0
        do i = 1, n
            s = s + v_in(i)**2
        end do
        s = sqrt(s)
        do i = 1, n
            v_out(i) = v_in(i) / s
        end do
    end subroutine norma
end module herramientas_fisica