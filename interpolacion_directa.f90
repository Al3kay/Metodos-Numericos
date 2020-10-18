program directpol
implicit none

real :: x(3),f(3),a(3,3),b(3),c(3),fxp
integer :: i
integer, parameter :: ndim=3,ndeg=2,n=3
real, parameter :: xp=3.44

x(1)=3.4 ; x(2)=3.5 ; x(3)=3.6
f(1)=0.294118 ; f(2)=0.285714 ; f(3)=0.277778

call directt(ndim,ndeg,n,x,fxp,f,xp,a,b,c)

do i=1,n
write(*,*) i,x(i),f(i)
enddo
print*,''
print*,'xp=',xp,'f(xp)=',fxp
print*,''
print*, c(1),'+',c(2),'x +',c(3),'x2'
end program

subroutine directt(ndim,ndeg,n,x,fxp,f,xp,a,b,c)
implicit none

integer :: ndim,ndeg,n,i
real :: x(3),f(3),a(3,3),b(3),c(3),xp,fxp

do i=1,n
a(i,1)=1.0
a(i,2) =x(i)
a(i,3) =x(i)**2
b(i)=f(i)
end do

call gauss (n,a,b,c)

fxp=c(1)+c(2)*xp+c(3)*xp**2
return
end subroutine

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
