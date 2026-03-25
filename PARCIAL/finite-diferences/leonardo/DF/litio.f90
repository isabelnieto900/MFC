module datos_fisicos
    implicit none
    real(8) :: z_global
end module datos_fisicos

program litio_variacional
    use optimizacion_mod
    use datos_fisicos
    implicit none
    
    integer, parameter :: ndim = 2
    real(8), parameter :: ftol = 1.0e-8
    integer :: i, iter
    real(8) :: fret
    real(8), dimension(ndim) :: p
    real(8), dimension(ndim, ndim) :: xi
    integer(8) :: t_start, t_stop, t_rate
    real(8) :: tiempo_total
    
    ! Inicialización de la matriz de direcciones
    xi = 0.0d0
    do i = 1, ndim
        xi(i,i) = 1.0d0
    end do
    
    ! Valores iniciales para Alfa y Beta
    p = (/ 1.5d0, 1.5d0 /)
    
    write(*,*) 'CUAL ES EL VALOR DE Z'
    read(*,*) z_global
    call system_clock(t_start, t_rate)
    ! LLAMADA CORREGIDA: Se pasa la función 'func' al final
    call powell(p, xi, ndim, ftol, iter, fret, func)
    call system_clock(t_stop)
    tiempo_total = dble(t_stop - t_start) / dble(t_rate)

    write(*,*) '--------------------------------------'
    write(*,*) 'ITERACIONES:        =', iter
    write(*,*) 'ALFA, BETA:         =', p(1), p(2)
    write(*,*) 'ENER. ESTADO BASE (eV) =', 27.211d0 * fret
    write(*,*) 'ENER. IONIZACION  (eV) =', 27.211d0 * ip_func(p(1), p(2))
    write(*,*) '--------------------------------------'
    write(*, '(A, F10.6, A)') "Tiempo de ejecución: ", tiempo_total, " segundos"

contains

    function func(x)
        real(8), dimension(:), intent(in) :: x
        real(8) :: func
        real(8) :: a, b, y
        
        a = x(1)
        b = x(2)
        y = 2.0d0 * a * b / ((2.0d0 * a + b)**5)
        
        func = a**2 - 2.0d0*z_global*a + (5.0d0/8.0d0)*a + (1.0d0/8.0d0)*b**2 - (z_global*b/4.0d0) + &
               y * (8.0d0*a**4 + 20.0d0*a**3*b + 12.0d0*a**2*b**2 + 10.0d0*a*b**3 + b**4)
    end function func

    function ip_func(alfa, beta)
        real(8), intent(in) :: alfa, beta
        real(8) :: ip_func
        real(8) :: a, b, a2, b2, y
        
        a = alfa
        b = beta
        a2 = a * a
        b2 = b * b
        y = 2.0d0 * a * b / ((2.0d0 * a + b)**5)
        
        ip_func = (1.0d0/8.0d0)*b2 - (z_global*b/4.0d0) + &
                  y * (8.0d0*a2*a2 + 20.0d0*a2*a*b + 12.0d0*a2*b2 + 10.0d0*a*b*b2 + b2*b2)
    end function ip_func
end program litio_variacional