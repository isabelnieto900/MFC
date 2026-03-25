module herramientas_sturm
    implicit none
    contains

    ! ------------------------------------------------------------------
    ! Cuenta el número de autovalores estrictamente menores que e_test
    ! ------------------------------------------------------------------
    function contar_sturm(n, diag, subdiag, e_test) result(count)
        integer, intent(in) :: n
        real(8), dimension(n), intent(in) :: diag, subdiag
        real(8), intent(in) :: e_test
        integer :: count, i
        real(8) :: q

        count = 0
        ! Primer elemento de la secuencia
        q = diag(1) - e_test
        if (q < 0.0d0) count = count + 1

        ! Recurrencia de Sturm (evita el overflow naturalmente)
        do i = 2, n
            if (abs(q) < 1.0d-30) q = 1.0d-30 ! Prevenir división por cero
            q = diag(i) - e_test - (subdiag(i-1)**2) / q
            if (q < 0.0d0) count = count + 1
        end do
    end function contar_sturm

    ! ------------------------------------------------------------------
    ! Encuentra la energía del estado "n_estado" usando bisección pura
    ! ------------------------------------------------------------------
    function biseccion_estado(n, diag, subdiag, n_estado, e_min_in, e_max_in, tol) result(raiz)
        integer, intent(in) :: n, n_estado
        real(8), dimension(n), intent(in) :: diag, subdiag
        real(8), intent(in) :: e_min_in, e_max_in, tol
        real(8) :: raiz, e_min, e_max, e_mid
        integer :: num_menores

        e_min = e_min_in
        e_max = e_max_in

        do while (abs(e_max - e_min) > tol)
            e_mid = 0.5d0 * (e_min + e_max)
            num_menores = contar_sturm(n, diag, subdiag, e_mid)

            if (num_menores >= n_estado) then
                ! El autovalor buscado está a la izquierda (o es e_mid)
                e_max = e_mid
            else
                ! El autovalor buscado está a la derecha
                e_min = e_mid
            end if
        end do
        raiz = 0.5d0 * (e_min + e_max)
    end function biseccion_estado

    ! ------------------------------------------------------------------
    ! Resolutor Tridiagonal (Thomas Algorithm) - De tu código original
    ! ------------------------------------------------------------------
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

    ! ------------------------------------------------------------------
    ! Normalización de función de onda - De tu código original
    ! ------------------------------------------------------------------
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

end module herramientas_sturm