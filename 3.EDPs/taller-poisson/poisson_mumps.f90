! poisson_mumps.f90 — Solución de la ecuación de Poisson con MUMPS (disperso PARALELO)
!
!   ∇²V = (x²+y²)·exp(x·y)  en Ω = [0,2]×[0,1]
!   Solución analítica: V(x,y) = exp(x·y)
!
! Compilación:
!   mpif90 -O2 -I/usr/include poisson_mumps.f90 -o poisson_mumps \
!          -ldmumps -lmumps_common -lpord -llapack -lblas
!
! Uso (P procesos MPI):
!   mpirun -np P ./poisson_mumps N

program poisson_mumps_main
    use mpi
    implicit none

    integer          :: N, argc, ierr, myrank
    character(len=32):: arg

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

    argc = command_argument_count()
    if (argc < 1) then
        if (myrank == 0) write(*,*) "Uso: mpirun -np P ./poisson_mumps N"
        call MPI_Finalize(ierr)
        stop
    end if
    call get_command_argument(1, arg)
    read(arg,*) N

    call solve_poisson(N, myrank)

    call MPI_Finalize(ierr)

end program poisson_mumps_main

! ─────────────────────────────────────────────────────────────────────────────
subroutine solve_poisson(N, myrank)
    use mpi
    implicit none
    include 'dmumps_struc.h'

    integer, intent(in) :: N, myrank

    double precision    :: f_rhs, V_exact, bc_left, bc_right, bc_bot, bc_top

    type(DMUMPS_STRUC)  :: id
    integer             :: Nx, Ny, sz, nnz_max, nnz
    integer             :: i, j, k, ierr
    double precision    :: hx, hy, ax, ay, ac
    double precision    :: x, y, diff, errL2, errMax, mem_MB
    double precision    :: t0, t1, elapsed_ms
    integer             :: mem_before_kb

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

    ! ── Solo el proceso 0 construye la matriz ────────────────────────────────
    if (myrank == 0) then

        mem_before_kb = 0
        block
            integer :: up, ios
            character(len=80) :: ln
            integer :: kb2
            open(newunit=up, file='/proc/self/status', status='old', action='read', iostat=ios)
            if (ios == 0) then
                do
                    read(up,'(A)',iostat=ios) ln
                    if (ios /= 0) exit
                    if (ln(1:6) == 'VmRSS:') then
                        read(ln(7:),*) kb2
                        mem_before_kb = kb2
                        exit
                    end if
                end do
                close(up)
            end if
        end block

        nnz_max = 5 * sz
        allocate(irn(nnz_max), jcn(nnz_max), aval(nnz_max))
        allocate(rhs(sz), sol_exact(sz))
        rhs = 0.0d0
        nnz = 0

        do j = 0, Ny-1
            y = dble(j+1)*hy
            do i = 0, Nx-1
                x = dble(i+1)*hx
                k = j*Nx + i + 1

                sol_exact(k) = V_exact(x, y)

                nnz = nnz + 1
                irn(nnz) = k;  jcn(nnz) = k;  aval(nnz) = ac

                if (i > 0) then
                    nnz = nnz + 1
                    irn(nnz) = k;  jcn(nnz) = k-1;  aval(nnz) = ax
                else
                    rhs(k) = rhs(k) - ax * bc_left(y)
                end if

                if (i < Nx-1) then
                    nnz = nnz + 1
                    irn(nnz) = k;  jcn(nnz) = k+1;  aval(nnz) = ax
                else
                    rhs(k) = rhs(k) - ax * bc_right(y)
                end if

                if (j > 0) then
                    nnz = nnz + 1
                    irn(nnz) = k;  jcn(nnz) = k-Nx;  aval(nnz) = ay
                else
                    rhs(k) = rhs(k) - ay * bc_bot(x)
                end if

                if (j < Ny-1) then
                    nnz = nnz + 1
                    irn(nnz) = k;  jcn(nnz) = k+Nx;  aval(nnz) = ay
                else
                    rhs(k) = rhs(k) - ay * bc_top(x)
                end if

                rhs(k) = rhs(k) + f_rhs(x, y)
            end do
        end do

    end if  ! myrank == 0

    ! ── Todos los procesos: inicializar MUMPS ────────────────────────────────
    id%JOB  = -1  ! inicializar
    id%SYM  = 0   ! no simétrico (genérico)
    id%PAR  = 1   ! host process también trabaja como worker

    call DMUMPS(id)

    ! Suprimir salida de diagnóstico en todos los rangos
    id%ICNTL(1) = -1
    id%ICNTL(2) = -1
    id%ICNTL(3) = -1
    id%ICNTL(4) =  0

    ! Solo rango 0 asigna la matriz y el RHS (entrada centralizada)
    if (myrank == 0) then
        id%N   = sz
        id%NNZ = nnz
        id%IRN => irn(1:nnz)
        id%JCN => jcn(1:nnz)
        id%A   => aval(1:nnz)
        id%RHS => rhs
    end if

    ! ── Todos los procesos: analizar + factorizar + resolver ─────────────────
    t0 = MPI_Wtime()
    id%JOB = 6
    call DMUMPS(id)
    t1 = MPI_Wtime()
    elapsed_ms = (t1 - t0) * 1.0d3

    ! ── Solo rango 0: postproceso, errores e I/O ─────────────────────────────
    if (myrank == 0) then

        if (id%INFO(1) /= 0) then
            write(*,*) "[MUMPS] Error INFO(1)=", id%INFO(1), " INFO(2)=", id%INFO(2)
            id%JOB = -2;  call DMUMPS(id)
            stop
        end if

        ! Memoria real del proceso (VmRSS de /proc/self/status)
        block
            integer :: unit_proc
            character(len=80) :: line
            integer :: kb
            mem_MB = 0.0d0
            open(newunit=unit_proc, file='/proc/self/status', status='old', action='read', iostat=kb)
            if (kb == 0) then
                do
                    read(unit_proc,'(A)',iostat=kb) line
                    if (kb /= 0) exit
                    if (line(1:6) == 'VmRSS:') then
                        read(line(7:),*) kb
                        mem_MB = dble(kb - mem_before_kb) / 1024.0d0
                        exit
                    end if
                end do
                close(unit_proc)
            end if
        end block

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
        write(*,'(A,F10.4,A)') "  Mem. RSS  : ", mem_MB, " MB"
        write(*,'(A,ES12.4)')  "  Error L2  : ", errL2
        write(*,'(A,ES12.4)')  "  Error Max : ", errMax

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

        fname_met = "metrics_fortran.csv"
        inquire(file=trim(fname_met), exist=file_exists)
        open(unit=12, file=trim(fname_met), status="unknown", position="append")
        if (.not. file_exists) then
            write(12,'(A)') "N,unknowns,time_ms,mem_MB,errL2,errMax"
        end if
        write(12,'(I0,A,I0,A,F12.4,A,F10.6,A,ES14.6,A,ES14.6)') &
            N,",",sz,",",elapsed_ms,",",mem_MB,",",errL2,",",errMax
        close(12)

        deallocate(irn, jcn, aval, rhs, sol_exact)

    end if  ! myrank == 0

    ! Todos los procesos finalizan MUMPS
    id%JOB = -2
    call DMUMPS(id)

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
