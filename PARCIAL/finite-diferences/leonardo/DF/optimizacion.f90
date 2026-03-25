module optimizacion_mod
    implicit none
    contains

    subroutine powell(p, xi, n, ftol, iter, fret, func)
        integer, intent(in) :: n
        real(8), dimension(n), intent(inout) :: p
        real(8), dimension(n, n), intent(inout) :: xi
        real(8), intent(in) :: ftol
        integer, intent(out) :: iter
        real(8), intent(out) :: fret
        interface
            function func(x)
                real(8), dimension(:), intent(in) :: x
                real(8) :: func
            end function func
        end interface

        integer, parameter :: itmax = 200
        real(8), dimension(n) :: pt, ptt, xit
        real(8) :: fptt, del, fp
        integer :: i, ibig

        fret = func(p)
        pt = p
        iter = 0

        do
            iter = iter + 1
            fp = fret
            ibig = 0
            del = 0.0d0
            do i = 1, n
                xit = xi(:,i)
                fptt = fret
                call linmin(p, xit, n, fret, func)
                if (abs(fptt - fret) > del) then
                    del = abs(fptt - fret)
                    ibig = i
                end if
            end do
            if (2.0d0 * abs(fp - fret) <= ftol * (abs(fp) + abs(fret)) + 1.0d-20) return
            if (iter == itmax) return
            ptt = 2.0d0 * p - pt
            xit = p - pt
            fptt = func(ptt)
            if (fptt < fp) then
                if (2.0d0*(fp-2.0d0*fret+fptt)*(fp-fret-del)**2 < del*(fp-fptt)**2) then
                    call linmin(p, xit, n, fret, func)
                    xi(:,ibig) = xi(:,n)
                    xi(:,n) = xit
                end if
            end if
            pt = p
        end do
    end subroutine powell

    subroutine linmin(p, xi, n, fret, func)
        integer, intent(in) :: n
        real(8), dimension(n), intent(inout) :: p
        real(8), dimension(n), intent(inout) :: xi
        real(8), intent(out) :: fret
        interface
            function func(x)
                real(8), dimension(:), intent(in) :: x
                real(8) :: func
            end function func
        end interface

        real(8) :: ax, xx, bx, fa, fx, fb, xmin
        real(8), dimension(n) :: pcom, xicom

        pcom = p
        xicom = xi
        ax = 0.0d0
        xx = 1.0d0
        call mnbrack(ax, xx, bx, fa, fx, fb, pcom, xicom, n, func)
        fret = brent_opt(ax, xx, bx, pcom, xicom, n, xmin, func)
        xi = xi * xmin
        p = p + xi
    end subroutine linmin

    subroutine mnbrack(ax, bx, cx, fa, fb, fc, pcom, xicom, n, func)
        real(8), intent(inout) :: ax, bx
        real(8), intent(out) :: cx, fa, fb, fc
        real(8), dimension(n), intent(in) :: pcom, xicom
        integer, intent(in) :: n
        interface
            function func(x)
                real(8), dimension(:), intent(in) :: x
                real(8) :: func
            end function func
        end interface
        real(8) :: gold, glimit, tiny, r, q, u, fu

        gold = 1.618034d0; glimit = 100.0d0; tiny = 1.0d-20
        fa = func(pcom + ax * xicom)
        fb = func(pcom + bx * xicom)
        if (fb > fa) then
            u = ax; ax = bx; bx = u
            u = fa; fa = fb; fb = u
        end if
        cx = bx + gold * (bx - ax)
        fc = func(pcom + cx * xicom)
        do while (fb > fc)
            r = (bx - ax) * (fb - fc)
            q = (bx - cx) * (fb - fa)
            u = bx - ((bx-cx)*q - (bx-ax)*r) / (2.0d0*sign(max(abs(q-r), tiny), q-r))
            fu = func(pcom + u * xicom)
            if ((bx-u)*(u-cx) > 0.0d0) then
                if (fu < fc) then
                    ax = bx; bx = u; fa = fb; fb = fu; return
                else if (fu > fb) then
                    cx = u; fc = fu; return
                end if
            end if
            ax = bx; bx = cx; cx = cx + gold * (cx - bx)
            fa = fb; fb = fc; fc = func(pcom + cx * xicom)
        end do
    end subroutine mnbrack

    function brent_opt(ax, bx, cx, pcom, xicom, n, xmin, func)
        real(8), intent(in) :: ax, bx, cx
        real(8), dimension(n), intent(in) :: pcom, xicom
        integer, intent(in) :: n
        real(8), intent(out) :: xmin
        real(8) :: brent_opt
        interface
            function func(x)
                real(8), dimension(:), intent(in) :: x
                real(8) :: func
            end function func
        end interface
        integer, parameter :: itmax = 100
        real(8), parameter :: cgld = 0.3819660d0, eps = 3.0d-8
        real(8) :: a, b, d, e, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm
        integer :: iter

        a = min(ax, cx); b = max(ax, cx); x = bx; w = bx; v = bx
        fw = func(pcom + x * xicom); fv = fw; fx = fw
        e = 0.0d0
        do iter = 1, itmax
            xm = 0.5d0 * (a + b)
            tol1 = eps * abs(x) + 1.0d-10; tol2 = 2.0d0 * tol1
            if (abs(x - xm) <= (tol2 - 0.5d0 * (b - a))) exit
            if (abs(e) > tol1) then
                r = (x - w) * (fx - fv); q = (x - v) * (fx - fw)
                p = (x - v) * q - (x - w) * r; q = 2.0d0 * (q - r)
                if (q > 0.0d0) p = -p
                q = abs(q)
                etemp = e; e = d
                if (abs(p) >= abs(0.5d0*q*etemp) .or. p <= q*(a-x) .or. p >= q*(b-x)) then
                    e = merge(a-x, b-x, x >= xm); d = cgld * e
                else
                    d = p / q; u = x + d
                    if (u - a < tol2 .or. b - u < tol2) d = sign(tol1, xm - x)
                end if
            else
                e = merge(a-x, b-x, x >= xm); d = cgld * e
            end if
            u = merge(x + d, x + sign(tol1, d), abs(d) >= tol1)
            fu = func(pcom + u * xicom)
            if (fu <= fx) then
                if (u >= x) then; a = x; else; b = x; endif
                v = w; fv = fw; w = x; fw = fx; x = u; fx = fu
            else
                if (u < x) then; a = u; else; b = u; endif
                if (fu <= fw .or. w == x) then
                    v = w; fv = fw; w = u; fw = fu
                else if (fu <= fv .or. v == x .or. v == w) then
                    v = u; fv = fu
                end if
            end if
        end do
        xmin = x; brent_opt = fx
    end function brent_opt
end module optimizacion_mod