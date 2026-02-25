! Diferencias Finitas - Ec. de Schrodinger: V(x) = dsin(xi)
PROGRAM PotSin
  implicit none
  integer::i,j,N,Nf
  real(8),allocatable::A(:,:),d(:),e(:),tem(:)
  real(8)::h,xi,h_2,Rmin,Rmax
  Rmin=-5.0d0; Rmax=5.0d0
  write(6,*)'6 primeros autovalores (E = lambda/2):'
  N=50; Nf=1000
  do while(N.le.Nf)
     h=(Rmax-Rmin)/N; h_2=1.d0/h**2
     allocate(A(N,N),d(N),e(N),tem(N))
     A=0.d0; d=0.d0; e=0.d0
     do i=1,N
        xi=Rmin+i*h
        d(i)=2.d0*h_2 + dsin(xi)
        e(i)=-h_2
        A(i,i)=1.d0
     end do
     call diagotri(d,e,N,A,.true.)
     do i=1,N
        do j=i+1,N
           if(d(j).lt.d(i))then
              xi=d(i); d(i)=d(j); d(j)=xi
              tem(:)=A(:,i); A(:,i)=A(:,j); A(:,j)=tem(:)
           end if
        end do
     end do
     write(6,"(i5,6(2x,F10.6))")N,0.5d0*d(1:6)
     if(N.ge.Nf/2)then
        open(unit=1,file="data_sin")
        write(1,"('# N=',i5,2x,6(F12.8,1x))")N,0.5d0*d(1:6)
        tem=0.d0
        do i=1,N
           tem(:)=tem(:)+A(i,:)**2
        end do
        tem=h*tem
        do i=1,N
           xi=Rmin+i*h
           write(1,"(7(F12.7,1x))")xi,A(i,1:6)**2/tem(i)
        enddo
        close(1)
     end if
     deallocate(A,d,e,tem)
     N=2*N
  end do
END PROGRAM PotSin
