program benchmark
    implicit none
    real(8), parameter :: b=0.02d0, d=0.015d0, r=0.1d0, p0=0.01d0, t_max=50.0d0
    integer, parameter :: n_steps = 10000000
    integer, parameter :: n_vis = 100
    real(8), parameter :: h = t_max / n_steps
    real(8), parameter :: h_vis = t_max / n_vis
    real(8) :: t_euler, t_taylor, t_trapecio
    real(8), allocatable :: p_euler(:), p_taylor(:), p_trapecio(:)
    integer :: unit, i
    real(8) :: t_val
    
    print *, 'Benchmark Fortran - 10^7 pasos'
    
    t_euler = euler(p0, h, n_steps)
    print '(A,ES12.4,A)', 'Euler:    ', t_euler, ' s'
    
    t_taylor = taylor2(p0, h, n_steps)
    print '(A,ES12.4,A)', 'Taylor2:  ', t_taylor, ' s'
    
    t_trapecio = trapecio(p0, h, n_steps)
    print '(A,ES12.4,A)', 'Trapecio: ', t_trapecio, ' s'
    
    if (t_euler < 0.0d0 .or. t_taylor < 0.0d0 .or. t_trapecio < 0.0d0) stop
    
    open(newunit=unit, file='benchmark_fortran.dat', status='replace')
    write(unit, '(A,F20.10)') 'Euler ', t_euler
    write(unit, '(A,F20.10)') 'Taylor2 ', t_taylor
    write(unit, '(A,F20.10)') 'Trapecio ', t_trapecio
    close(unit)
    
    ! Calcular soluciones para graficar
    allocate(p_euler(0:n_vis), p_taylor(0:n_vis), p_trapecio(0:n_vis))
    
    call euler_solucion(p0, h_vis, n_vis, p_euler)
    call taylor2_solucion(p0, h_vis, n_vis, p_taylor)
    call trapecio_solucion(p0, h_vis, n_vis, p_trapecio)
    
    open(newunit=unit, file='soluciones_fortran.dat', status='replace')
    write(unit, '(A)') '# t Exacta Euler Taylor2 Trapecio'
    do i = 0, n_vis
        t_val = i * h_vis
        write(unit, '(5F20.10)') t_val, exacta(t_val), p_euler(i), &
                                  p_taylor(i), p_trapecio(i)
    end do
    close(unit)
    
    deallocate(p_euler, p_taylor, p_trapecio)

contains

    pure function f(p) result(fp)
        real(8), intent(in) :: p
        real(8) :: fp
        fp = r * b * (1.0d0 - p)
    end function

    pure function f_prime(p) result(fpp)
        real(8), intent(in) :: p
        real(8) :: fpp
        fpp = -r * b * f(p)
    end function

    function euler(p_init, h, n) result(tiempo)
        real(8), intent(in) :: p_init, h
        integer, intent(in) :: n
        real(8) :: tiempo, p, t1, t2
        integer :: i
        p = p_init
        call cpu_time(t1)
        do i = 1, n
            p = p + h * f(p)
        end do
        call cpu_time(t2)
        tiempo = t2 - t1
        if (p < 0.0d0 .or. p > 2.0d0) print *, 'guard:', p
    end function

    function taylor2(p_init, h, n) result(tiempo)
        real(8), intent(in) :: p_init, h
        integer, intent(in) :: n
        real(8) :: tiempo, p, t1, t2
        integer :: i
        p = p_init
        call cpu_time(t1)
        do i = 1, n
            p = p + h * f(p) + (h * h / 2.0d0) * f_prime(p)
        end do
        call cpu_time(t2)
        tiempo = t2 - t1
        if (p < 0.0d0 .or. p > 2.0d0) print *, 'guard:', p
    end function

    function trapecio(p_init, h, n) result(tiempo)
        real(8), intent(in) :: p_init, h
        integer, intent(in) :: n
        real(8) :: tiempo, p, factor, t1, t2
        integer :: i
        factor = h * r * b / 2.0d0
        p = p_init
        call cpu_time(t1)
        do i = 1, n
            p = (p + (h / 2.0d0) * f(p) + factor) / (1.0d0 + factor)
        end do
        call cpu_time(t2)
        tiempo = t2 - t1
        if (p < 0.0d0 .or. p > 2.0d0) print *, 'guard:', p
    end function

    pure function exacta(t) result(p)
        real(8), intent(in) :: t
        real(8) :: p
        p = 1.0d0 - (1.0d0 - p0) * exp(-r * b * t)
    end function

    subroutine euler_solucion(p_init, h, n, sol)
        real(8), intent(in) :: p_init, h
        integer, intent(in) :: n
        real(8), intent(out) :: sol(0:n)
        real(8) :: p
        integer :: i
        p = p_init
        sol(0) = p
        do i = 1, n
            p = p + h * f(p)
            sol(i) = p
        end do
    end subroutine

    subroutine taylor2_solucion(p_init, h, n, sol)
        real(8), intent(in) :: p_init, h
        integer, intent(in) :: n
        real(8), intent(out) :: sol(0:n)
        real(8) :: p
        integer :: i
        p = p_init
        sol(0) = p
        do i = 1, n
            p = p + h * f(p) + (h * h / 2.0d0) * f_prime(p)
            sol(i) = p
        end do
    end subroutine

    subroutine trapecio_solucion(p_init, h, n, sol)
        real(8), intent(in) :: p_init, h
        integer, intent(in) :: n
        real(8), intent(out) :: sol(0:n)
        real(8) :: p, factor
        integer :: i
        p = p_init
        factor = h * r * b / 2.0d0
        sol(0) = p
        do i = 1, n
            p = (p + (h / 2.0d0) * f(p) + factor) / (1.0d0 + factor)
            sol(i) = p
        end do
    end subroutine

end program
