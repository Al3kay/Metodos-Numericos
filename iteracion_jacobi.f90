program LU
integer, parameter :: n=5
real, dimension(n,n) :: A
real, dimension(n) :: b,x2,x1
integer :: i

A(1,1)=4.0 ; A(1,2)=-1.0 ; A(1,3)=0.0 ; A(1,4)=1.0 ; A(1,5)=0.0
A(2,1)=-1.0 ; A(2,2)=4.0 ; A(2,3)=-1.0 ; A(2,4)=0.0 ; A(2,5)=1.0
A(3,1)=0.0 ; A(3,2)=-1.0 ; A(3,3)=4.0 ; A(3,4)=-1.0 ; A(3,5)=0.0
A(4,1)=1.0 ; A(4,2)=0.0 ; A(4,3)=-1.0 ; A(4,4)=4.0 ; A(4,5)=-1.0
A(5,1)=0.0 ; A(5,2)=1.0 ; A(5,3)=0.0 ; A(5,4)=-1.0 ; A(5,5)=4.0

b(1)=100.0 ; b(2)=100.0  ; b(3)=100.0 ; b(4)=100.0 ; b(5)=100.0

x1=0.0

call solve(n,A,b,x2,x1)

print*, 'Solución: '
do i=1,n
	write(*,*) i,x2(i)
enddo

end program

subroutine solve(n,A,b,x2,x1)
implicit none

integer :: n,k,i,j,kmax
real :: a(n,n),x2(n),x1(n),b(n)
real, dimension(:), allocatable :: R

allocate(R(n))

kmax=200
!Este paso determina la matriz L y U, ellos se guardan en A, es decir
!la matriz A es reescrita y podemos encontrar implicitamente a L y U 
do k=1,kmax
R=0.0
	do i=1,n
		R(i)=R(i)+b(i)
		do j=1,n
			R(i)=R(i)-a(i,j)*x1(j)
		enddo
		x2(i)=x1(i)+(R(i)/a(i,i))	
	enddo
x1=x2
enddo

deallocate(R)
return
end subroutine
