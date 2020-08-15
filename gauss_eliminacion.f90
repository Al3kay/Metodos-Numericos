program gausselim
integer, parameter :: n=3
real, dimension(n,n) :: A
real, dimension(n) :: b,x
integer :: i

A(1,1)=80.0 ; A(1,2)=-20.0 ; A(1,3)=-20.0
A(2,1)=-20.0 ; A(2,2)=40.0 ; A(2,3)=-20.0
A(3,1)=-20.0 ; A(3,2)=-20.0 ; A(3,3)=130.0


b(1)=20.0 ; b(2)=20.0  ; b(3)=20.0

call gauss(n,a,b,x)

print*, 'Solución: '
do i=1,n
	write(*,*) i,x(i)
enddo

end program

subroutine gauss(n,a,b,x)
implicit none
integer :: n,k,i,j
real :: a(n,n),b(n),x(n),em,sumj

do k=1,n-1
	do j=k+1,n 
		do i=k+1,n 	
	!--- Modifica todos los elementos arriba de la diagonal
			em=a(k,j)/a(k,k)
			a(i,j)=a(i,j)-a(i,k)*(em) 
		enddo
	enddo
enddo

do k=1,n-1
	do i=k+1,n
	!------ Modificia b, el último (n) es el deseado para x_n 
		b(i)=b(i)-(a(i,k)/a(k,k))*b(k)
	enddo
enddo

x(n)=b(n)/a(n,n)
do i=n-1,1,-1
	do j=n,i+1,-1
		x(i)=x(i)-a(i,j)*x(j)
	enddo
	x(i)=x(i)+b(i)
	x(i)=x(i)/a(i,i)
enddo


return
end subroutine 
