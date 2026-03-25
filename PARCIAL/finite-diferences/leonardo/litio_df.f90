program litio_espectro
    use herramientas_fisica
    implicit none
    
    integer, parameter :: npt = 6000, nbmax = 1000
    real(8) :: xb1(nbmax), xb2(nbmax), c(npt)
    real(8) :: energia(nbmax), diag(npt), superd(npt), subd(npt)
    real(8) :: vk(npt), vkm1(npt)
    real(8) :: h, ele, d, x1, x2, tol, raiz, vp
    integer :: nn, nb, jj, i, j, k, n_divs
    
    real(8), external :: determ
    common /params/ h, ele, nn
    
    nn = npt
    d = 25.0d0
    
    write(*,*) 'MOMENTUM ANGULAR (l) ='
    read(*,*) ele
    
    h = d / dble(npt)
    x1 = -10.0d0
    x2 = 0.0d0
    
    ! Inicialización de búsqueda
    nb = nbmax
    n_divs = 20000  
    
    write(*,*) 'Buscando intervalos...'
    call zbrak(determ, x1, x2, n_divs, xb1, xb2, nb)
    
    if (nb == 0) then
        write(*,*) 'Error: No se encontraron intervalos.'
        stop
    end if
    
    jj = 0
    write(*,*) 'Buscando autovalores...'
    do i = 1, nb
        tol = 1.0d-12
        raiz = zbrent(determ, xb1(i), xb2(i), tol)
        
        if (raiz < -1.0d-8) then
            jj = jj + 1
            energia(jj) = raiz
            write(*,*) 'N=', jj, ' E(au)=', energia(jj), ' E(eV)=', 27.211d0 * energia(jj)
        end if
    end do
    
    ! Guardar el potencial previamente para optimizar
    do i = 1, npt
        c(i) = pot_local(h * dble(i), ele)
    end do
    
    open(10, file='LITIO.DAT')
    
    do j = 1, jj
        vp = energia(j)
        
        ! 1. Construcción de la matriz H - E (Tridiagonal)
        do k = 1, npt
            ! Esta expresión ya incluye la resta de la energía (-vp)
            diag(k) = 2.0d0 + 2.0d0 * h**2 * (c(k) - vp)
        end do
        superd = -1.0d0
        subd = -1.0d0
        
        ! 2. Iteración Inversa para hallar el autovector
        vk = 1.0d0  ! Semilla inicial
        do k = 1, 15
            ! Resolvemos el sistema (H-E) * v_k = v_{k-1}
            call tridag(subd, diag, superd, vk, vkm1, npt)
            call norma(npt, vkm1, vk)
        end do

        ! 3. CORRECCIÓN DE FASE (Evita gráficas invertidas)
        ! Si el primer punto es negativo, reflejamos toda la función
        if (vk(1) < 0.0d0) then
            vk = -vk
        end if
        vk = vk / sqrt(h)
        ! 4. Escritura en el archivo .DAT
        do i = 1, npt, 10
            write(10,*) dble(i) * h, vk(i)
        end do
        write(10,*) ' '  ! Separador para Gnuplot o Python
    end do
    
    close(10)
    write(*,*) '--------------------------------------------'
    write(*,*) 'Proceso finalizado. Datos en LITIO.DAT'

contains

    function pot_local(x, l)
        real(8), intent(in) :: x, l
        real(8) :: pot_local, alfa
        alfa = 2.535930d0  ! <-- VUELVE AL ALFA DE C++
        pot_local = -(1.0d0/x) - (2.0d0/x) * (1.0d0 + alfa*x) * exp(-2.0d0*alfa*x)
        pot_local = pot_local + l * (l + 1.0d0) / (2.0d0 * x**2)
    end function pot_local

end program litio_espectro

! --- Función Determinante ---
real(8) function determ(x)
    implicit none
    real(8), intent(in) :: x
    real(8) :: h, ele
    integer :: nn
    common /params/ h, ele, nn
    
    real(8) :: h2, p0, p1, p2, pot_val, r, alfa
    integer :: i

    alfa = 2.535930d0 ! <-- VUELVE AL ALFA DE C++
    h2 = h * h
    p0 = 0.0d0
    p1 = 1.0d0 

    do i = 1, nn ! <-- HASTA nn (INCLUSIVE COMO EN C++)
        r = h * dble(i)
        pot_val = -(1.0d0/r) - (2.0d0/r)*(1.0d0 + alfa*r)*exp(-2.0d0*alfa*r) + &
                  ele*(ele + 1.0d0)/(2.0d0*r**2)
        
        ! <-- VUELVE EL IF QUE TIENE C++
        if (i == 1) then
            p1 = (2.0d0 + 2.0d0 * h2 * pot_val - 2.0d0 * x * h2)
        else
            p2 = (2.0d0 + 2.0d0 * h2 * (pot_val - x)) * p1 - p0
            p0 = p1
            p1 = p2
        end if
        
        if (abs(p1) > 1.0d20) then
            p1 = p1 * 1.0d-20
            p0 = p0 * 1.0d-20
        end if
    end do
    determ = p1
end function determ