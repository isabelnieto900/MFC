program sturm_litio
    use herramientas_sturm
    implicit none
    
    integer, parameter :: npt = 6000
    integer, parameter :: max_estados = 4  ! Cuántos niveles quieres encontrar
    
    real(8) :: diag(npt), subdiag(npt), diag_shift(npt)
    real(8) :: superdiag(npt)
    real(8) :: vk(npt), vkm1(npt)
    real(8) :: h_step, ele, d_max, r, tol, energia
    integer :: i, estado, iter
    
    d_max = 25.0d0
    h_step = d_max / dble(npt)
    tol = 1.0d-12
    
    write(*,*) 'MOMENTUM ANGULAR (l) ='
    read(*,*) ele
    
    ! ---------------------------------------------------------
    ! 1. Construir la Matriz Hamiltoniana H
    ! H(i,i) = 1/h^2 + V(r)
    ! H(i,i+-1) = -1/(2h^2)
    ! ---------------------------------------------------------
    do i = 1, npt
        r = h_step * dble(i)
        diag(i) = (1.0d0 / (h_step**2)) + pot_local(r, ele)
        subdiag(i) = -1.0d0 / (2.0d0 * h_step**2)
        superdiag(i) = subdiag(i)
    end do
    
    open(10, file='LITIO_STURM.DAT')
    
    write(*,*) '---------------------------------------------------'
    write(*,*) ' Buscando energias usando la Secuencia de Sturm...'
    write(*,*) '---------------------------------------------------'
    
    ! ---------------------------------------------------------
    ! 2. Búsqueda y reconstrucción de cada estado
    ! ---------------------------------------------------------
    do estado = 1, max_estados
        
        ! Buscamos el autovalor "estado" en un intervalo seguro (ej: -50 a 0 au)
        energia = biseccion_estado(npt, diag, subdiag, estado, -50.0d0, 0.0d0, tol)
        
        write(*,*) 'N=', estado, ' E(au)=', energia, ' E(eV)=', 27.211d0 * energia
        
        ! --- Método de la potencia inversa para la función de onda ---
        ! (H - E*I) * v_new = v_old
        do i = 1, npt
            diag_shift(i) = diag(i) - energia
        end do
        
        vk = 1.0d0 ! Guess inicial plano
        
        do iter = 1, 15
            call tridag(subdiag, diag_shift, superdiag, vk, vkm1, npt)
            call norma(npt, vkm1, vk)
        end do
        
        do i = 2, npt
            if (abs(vk(i)) > 1.0d-5) then
                if (vk(i) < 0.0d0) then
                    vk = -vk  ! Si es negativo, lo volteamos
                end if
                exit ! Salimos del buscador de fase
            end if
        end do
        
        vk = vk / sqrt(h_step)

        ! Guardar en archivo
        do i = 1, npt, 10
            write(10,*) dble(i) * h_step, vk(i)
        end do
        write(10,*) ' ' 
        
    end do
    
    close(10)
    write(*,*) '---------------------------------------------------'
    write(*,*) 'Resultados guardados en LITIO_STURM.DAT'

contains

    ! Función potencial exacta a la que tenías [cite: 166]
    function pot_local(x, l)
        real(8), intent(in) :: x, l
        real(8) :: pot_local, alfa
        alfa = 2.535930d0
        pot_local = -(1.0d0/x) - (2.0d0/x) * (1.0d0 + alfa*x) * exp(-2.0d0*alfa*x)
        pot_local = pot_local + l * (l + 1.0d0) / (2.0d0 * x**2)
    end function pot_local

end program sturm_litio