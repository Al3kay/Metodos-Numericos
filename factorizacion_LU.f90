program LU
integer, parameter :: n=3
real, dimension(n,n) :: A
real, dimension(n) :: b,bp,x
integer :: i

A(1,1)=80.0 ; A(1,2)=-20.0 ; A(1,3)=-20.0
A(2,1)=-20.0 ; A(2,2)=40.0 ; A(2,3)=-20.0
A(3,1)=-20.0 ; A(3,2)=-20.0 ; A(3,3)=130.0


b(1)=20.0 ; b(2)=20.0  ; b(3)=20.0

call LUfactor(n,A)
call solve(n,A,b,bp,x)

print*, 'Solución: '
do i=1,n
	write(*,*) i,x(i)
enddo

end program

subroutine lufactor(n,A)
implicit none

integer :: n,k,i,j
real :: a(n,n),em

!Este paso determina la matriz L y U, ellos se guardan en A, es decir
!la matriz A es reescrita y podemos encontrar implicitamente a L y U 
do k=1,n-1
	do i=k+1,n
		em=a(i,k)/a(k,k)
		a(i,k)=em
		do j=k+1,n
			a(i,j)=a(i,j)-em*a(k,j)
		enddo
	enddo
enddo

return
end subroutine

subroutine solve(n,A,b,bp,x)
implicit none
integer :: n,k,i,j
real :: a(n,n),bp(n),b(n),x(n),em,sumj

bp(1)=b(1)

!Forward step para calcular b'
do i=2,n
	bp(i)=b(i)
	do j=1,i-1
		bp(i)=bp(i)-a(i,j)*bp(j)
	enddo
enddo

!Backward para contrar x solucion
x(n)=bp(n)/a(n,n)
do i=n-1,1,-1
	do j=n,i+1,-1
		x(i)=x(i)-a(i,j)*x(j)
	enddo
	x(i)=x(i)+bp(i)
	x(i)=x(i)/a(i,i)
enddo
return
end subroutine 
