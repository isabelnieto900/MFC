! poisson.f90 — Solución de la ecuación de Poisson con LAPACK dgbsv
!
!   ∇²V = (x²+y²)·exp(x·y)  en Ω = [0,2]×[0,1]
!   Solución analítica: V(x,y) = exp(x·y)
!
! Compilación:
!   gfortran -O2 poisson.f90 -o poisson_f -llapack -lblas
!
! Uso:
!   ./poisson_f N

program poisson_fd
    implicit none

    integer          :: N, argc
    character(len=32):: arg

    argc = command_argument_count()
    if (argc < 1) then
        write(*,*) "Uso: ./poisson_f N"
        stop
    end if
    call get_command_argument(1, arg)
    read(arg,*) N

    call solve_poisson(N)

end program poisson_fd

! ─────────────────────────────────────────────────────────────────────────────
subroutine solve_poisson(N)
    implicit none
    integer, intent(in) :: N

    ! Declaraciones explícitas de funciones externas
    double precision    :: f_rhs, V_exact, bc_left, bc_right, bc_bot, bc_top

    integer             :: M, sz, kl, ku, ldab
    integer             :: i, j, k, ki, kj, info
    double precision    :: hx, hy, ax, ay, ac
    double precision    :: x, y, diff, errL2, errMax, mem_MB
    double precision    :: t0, t1, elapsed_ms

    ! Variables para LAPACK dgbsv
    integer, allocatable          :: ipiv(:)
    double precision, allocatable :: AB(:,:), b_vec(:)

    ! Para métricas
    character(len=64)  :: fname_csv, fname_met

    M    = N
    sz   = N * M
    hx   = 2.0d0 / dble(N + 1)
    hy   = 1.0d0 / dble(N + 1)
    ax   = 1.0d0 / (hx*hx)
    ay   = 1.0d0 / (hy*hy)
    ac   = -2.0d0*ax - 2.0d0*ay

    ! Ancho de banda: kl = ku = N  (vecinos en y a distancia N)
    kl   = N
    ku   = N
    ldab = 2*kl + ku + 1

    allocate(AB(ldab, sz))
    allocate(b_vec(sz))
    allocate(ipiv(sz))
    AB     = 0.0d0
    b_vec  = 0.0d0

    ! Tiempo de inicio (cpu_time)
    call cpu_time(t0)

    ! ── Construcción de la matriz en formato banda LAPACK ──────────────────
    ! AB(kl+ku+1 + i - j, j) = A(i,j)   → offset = kl+ku+1
    do j = 0, M-1
        y = dble(j+1)*hy
        do i = 0, N-1
            x = dble(i+1)*hx
            k = j*N + i + 1          ! índice Fortran (1-base)

            ! diagonal
            AB(kl+ku+1, k) = ac

            ! vecino izquierdo (k-1)
            if (i > 0) then
                ki = k - 1
                AB(kl+ku+1 + k - ki, ki) = ax    ! AB(kl+ku+2, k-1)
            else
                b_vec(k) = b_vec(k) - ax * bc_left(y)
            end if

            ! vecino derecho (k+1)
            if (i < N-1) then
                ki = k + 1
                AB(kl+ku+1 + k - ki, ki) = ax    ! AB(kl+ku, k+1)
            else
                b_vec(k) = b_vec(k) - ax * bc_right(y)
            end if

            ! vecino inferior (k-N)
            if (j > 0) then
                kj = k - N
                AB(kl+ku+1 + k - kj, kj) = ay
            else
                b_vec(k) = b_vec(k) - ay * bc_bot(x)
            end if

            ! vecino superior (k+N)
            if (j < M-1) then
                kj = k + N
                AB(kl+ku+1 + k - kj, kj) = ay
            else
                b_vec(k) = b_vec(k) - ay * bc_top(x)
            end if

            ! fuente
            b_vec(k) = b_vec(k) + f_rhs(x, y)
        end do
    end do

    ! ── Resolver con LAPACK dgbsv ──────────────────────────────────────────
    call dgbsv(sz, kl, ku, 1, AB, ldab, ipiv, b_vec, sz, info)

    call cpu_time(t1)
    elapsed_ms = (t1 - t0) * 1.0d3

    if (info /= 0) then
        write(*,*) "[Fortran] Error en dgbsv: info =", info
        stop
    end if

    ! ── Errores ──────────────────────────────────────────────────────────────
    errL2  = 0.0d0
    errMax = 0.0d0
    do j = 0, M-1
        y = dble(j+1)*hy
        do i = 0, N-1
            x    = dble(i+1)*hx
            k    = j*N + i + 1
            diff = dabs(b_vec(k) - V_exact(x, y))
            errL2  = errL2 + diff*diff
            if (diff > errMax) errMax = diff
        end do
    end do
    errL2 = dsqrt(errL2 / dble(N*M))

    ! Memoria estimada: entradas de la banda con valores ~ 5*sz × 20 bytes
    mem_MB = dble(ldab * sz) * 8.0d0 / (1024.0d0*1024.0d0)

    write(*,'(A,I0,A,I0,A,I0,A)') "Fortran/LAPACK  malla ", N, "x", N, &
                                    "  (", N*N, " incognitas)"
    write(*,'(A,F12.4,A)') "  Tiempo    : ", elapsed_ms, " ms"
    write(*,'(A,F10.4,A)') "  Mem. est. : ", mem_MB, " MB"
    write(*,'(A,ES12.4)')  "  Error L2  : ", errL2
    write(*,'(A,ES12.4)')  "  Error Max : ", errMax

    ! ── CSV de la solución ────────────────────────────────────────────────────
    write(fname_csv,'(A,I0,A)') "sol_", N, "_fortran.csv"
    open(unit=10, file=trim(fname_csv), status="replace")
    write(10,'(A)') "x,y,V_num,V_exact"
    do j = 0, M-1
        y = dble(j+1)*hy
        do i = 0, N-1
            x = dble(i+1)*hx
            k = j*N + i + 1
            write(10,'(ES16.8,A,ES16.8,A,ES16.8,A,ES16.8)') &
                x,",",y,",",b_vec(k),",",V_exact(x,y)
        end do
    end do
    close(10)
    write(*,'(A,A)') "  CSV: ", trim(fname_csv)

    ! ── Métricas acumulativas ────────────────────────────────────────────────
    fname_met = "metrics_fortran.csv"
    ! Simplemente append; el notebook limpiará antes si es necesario
    open(unit=11, file=trim(fname_met), status="unknown", position="append")
    write(11,'(I0,A,I0,A,F12.4,A,F10.6,A,ES14.6,A,ES14.6)') &
        N,",",N*N,",",elapsed_ms,",",mem_MB,",",errL2,",",errMax
    close(11)

    deallocate(AB, b_vec, ipiv)

end subroutine solve_poisson

! ─────────────────────────────────────────────────────────────────────────────
double precision function f_rhs(x, y)
    double precision, intent(in) :: x, y
    f_rhs = (x*x + y*y) * dexp(x*y)
end function

double precision function V_exact(x, y)
    double precision, intent(in) :: x, y
    V_exact = dexp(x*y)
end function

double precision function bc_left(y)
    double precision, intent(in) :: y
    bc_left = 1.0d0
end function

double precision function bc_right(y)
    double precision, intent(in) :: y
    bc_right = dexp(2.0d0*y)
end function

double precision function bc_bot(x)
    double precision, intent(in) :: x
    bc_bot = 1.0d0
end function

double precision function bc_top(x)
    double precision, intent(in) :: x
    bc_top = dexp(x)
end function
