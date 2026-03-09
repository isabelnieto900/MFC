! poisson_mumps.f90 — Solución de la ecuación de Poisson con MUMPS (disperso)
!
!   ∇²V = (x²+y²)·exp(x·y)  en Ω = [0,2]×[0,1]
!   Solución analítica: V(x,y) = exp(x·y)
!
! Compilación:
!   mpif90 -O2 poisson_mumps.f90 -o poisson_mumps -ldmumps -lmumps_common -lpord -llapack -lblas
!
! Uso:
!   ./poisson_mumps N

program poisson_mumps_main
    use mpi
    implicit none

    integer          :: N, argc, ierr
    character(len=32):: arg

    call MPI_Init(ierr)

    argc = command_argument_count()
    if (argc < 1) then
        write(*,*) "Uso: ./poisson_mumps N"
        call MPI_Finalize(ierr)
        stop
    end if
    call get_command_argument(1, arg)
    read(arg,*) N

    call solve_poisson(N)

    call MPI_Finalize(ierr)

end program poisson_mumps_main

! ─────────────────────────────────────────────────────────────────────────────
subroutine solve_poisson(N)
    use mpi
    implicit none
    include 'dmumps_struc.h'

    integer, intent(in) :: N

    double precision    :: f_rhs, V_exact, bc_left, bc_right, bc_bot, bc_top

    type(DMUMPS_STRUC)  :: id
    integer             :: Nx, Ny, sz, nnz_max, nnz
    integer             :: i, j, k, ierr
    double precision    :: hx, hy, ax, ay, ac
    double precision    :: x, y, diff, errL2, errMax, mem_MB
    double precision    :: t0, t1, elapsed_ms

    integer,          allocatable, target :: irn(:), jcn(:)
    double precision, allocatable, target :: aval(:), rhs(:), sol_exact(:)

    character(len=64) :: fname_csv, fname_met
    logical           :: file_exists

    Nx = N
    Ny = N
    sz = Nx * Ny

    hx = 2.0d0 / dble(Nx + 1)
    hy = 1.0d0 / dble(Ny + 1)
    ax = 1.0d0 / (hx*hx)
    ay = 1.0d0 / (hy*hy)
    ac = -2.0d0*ax - 2.0d0*ay

    ! Reservar COO con 5 entradas por fila como máximo
    nnz_max = 5 * sz
    allocate(irn(nnz_max), jcn(nnz_max), aval(nnz_max))
    allocate(rhs(sz), sol_exact(sz))
    rhs   = 0.0d0
    nnz   = 0

    call cpu_time(t0)

    ! ── Construcción en formato COO (índices 1-base para MUMPS) ─────────────
    do j = 0, Ny-1
        y = dble(j+1)*hy
        do i = 0, Nx-1
            x = dble(i+1)*hx
            k = j*Nx + i + 1          ! índice 1-base

            sol_exact(k) = V_exact(x, y)

            ! Diagonal
            nnz = nnz + 1
            irn(nnz) = k;  jcn(nnz) = k;  aval(nnz) = ac

            ! Vecino izquierdo
            if (i > 0) then
                nnz = nnz + 1
                irn(nnz) = k;  jcn(nnz) = k-1;  aval(nnz) = ax
            else
                rhs(k) = rhs(k) - ax * bc_left(y)
            end if

            ! Vecino derecho
            if (i < Nx-1) then
                nnz = nnz + 1
                irn(nnz) = k;  jcn(nnz) = k+1;  aval(nnz) = ax
            else
                rhs(k) = rhs(k) - ax * bc_right(y)
            end if

            ! Vecino inferior
            if (j > 0) then
                nnz = nnz + 1
                irn(nnz) = k;  jcn(nnz) = k-Nx;  aval(nnz) = ay
            else
                rhs(k) = rhs(k) - ay * bc_bot(x)
            end if

            ! Vecino superior
            if (j < Ny-1) then
                nnz = nnz + 1
                irn(nnz) = k;  jcn(nnz) = k+Nx;  aval(nnz) = ay
            else
                rhs(k) = rhs(k) - ay * bc_top(x)
            end if

            ! Fuente
            rhs(k) = rhs(k) + f_rhs(x, y)
        end do
    end do

    ! ── Configurar MUMPS (modo secuencial, PAR=1) ────────────────────────────
    id%JOB  = -1  ! inicializar
    id%SYM  = 0   ! no simétrico (genérico)
    id%PAR  = 1   ! proceso anfitrión también trabaja

    call DMUMPS(id)

    ! Suprimir salida de diagnóstico
    id%ICNTL(1) = -1
    id%ICNTL(2) = -1
    id%ICNTL(3) = -1
    id%ICNTL(4) =  0

    id%N   = sz
    id%NNZ = nnz
    id%IRN => irn(1:nnz)
    id%JCN => jcn(1:nnz)
    id%A   => aval(1:nnz)
    id%RHS => rhs

    ! Analizar + factorizar + resolver
    id%JOB = 6
    call DMUMPS(id)

    call cpu_time(t1)
    elapsed_ms = (t1 - t0) * 1.0d3

    if (id%INFO(1) /= 0) then
        write(*,*) "[MUMPS] Error INFO(1)=", id%INFO(1), " INFO(2)=", id%INFO(2)
        id%JOB = -2;  call DMUMPS(id)
        stop
    end if

    ! Memoria real usada por MUMPS (en MB)
    mem_MB = dble(id%INFOG(22))   ! [MB] según documentación MUMPS

    ! Finalizar MUMPS
    id%JOB = -2
    call DMUMPS(id)

    ! ── Errores ─────────────────────────────────────────────────────────────
    errL2  = 0.0d0
    errMax = 0.0d0
    do k = 1, sz
        diff = dabs(rhs(k) - sol_exact(k))   ! rhs contiene la solución tras DMUMPS
        errL2  = errL2 + diff*diff
        if (diff > errMax) errMax = diff
    end do
    errL2 = dsqrt(errL2 / dble(sz))

    write(*,'(A,I0,A,I0,A,I0,A)') "Fortran/MUMPS  malla ", N, "x", N, &
                                    "  (", sz, " incognitas)"
    write(*,'(A,F12.4,A)') "  Tiempo    : ", elapsed_ms, " ms"
    write(*,'(A,F10.4,A)') "  Mem. est. : ", mem_MB, " MB"
    write(*,'(A,ES12.4)')  "  Error L2  : ", errL2
    write(*,'(A,ES12.4)')  "  Error Max : ", errMax

    ! ── CSV de la solución ───────────────────────────────────────────────────
    write(fname_csv,'(A,I0,A)') "sol_", N, "_fortran.csv"
    open(unit=10, file=trim(fname_csv), status="replace")
    write(10,'(A)') "x,y,V_num,V_exact"
    do j = 0, Ny-1
        y = dble(j+1)*hy
        do i = 0, Nx-1
            x = dble(i+1)*hx
            k = j*Nx + i + 1
            write(10,'(ES16.8,A,ES16.8,A,ES16.8,A,ES16.8)') &
                x,",",y,",",rhs(k),",",sol_exact(k)
        end do
    end do
    close(10)
    write(*,'(A,A)') "  CSV solucion: ", trim(fname_csv)

    ! ── Métricas ─────────────────────────────────────────────────────────────
    fname_met = "metrics_fortran.csv"
    ! Escribir cabecera solo si el archivo no existe
    inquire(file=trim(fname_met), exist=file_exists)
    open(unit=12, file=trim(fname_met), status="unknown", position="append")
    if (.not. file_exists) then
        write(12,'(A)') "N,unknowns,time_ms,mem_MB,errL2,errMax"
    end if
    write(12,'(I0,A,I0,A,F12.4,A,F10.6,A,ES14.6,A,ES14.6)') &
        N,",",sz,",",elapsed_ms,",",mem_MB,",",errL2,",",errMax
    close(12)

    deallocate(irn, jcn, aval, rhs, sol_exact)

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
