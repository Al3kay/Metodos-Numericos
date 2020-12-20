program diferencias
implicit none

real :: xk,yk,dx
integer :: reason,ndata,n,i,j,k,imin,imax
real, dimension(:), allocatable :: xx,yy,dydx
character*512 :: IFILE

ndata=0
IFILE='data.dat'

open(1,file=trim(IFILE))
do
	read(1,*,IOSTAT=reason) xk,yk	
	if(reason > 0) then		
		write(*,*) 'Something wrong'
		exit
	else if (reason < 0) then
		write(*,*) 'EoF reached'
		exit
	else
		ndata=ndata+1
	end if
enddo
close(1)

open(2,file=trim(IFILE))
print*, 'Number of Data:',ndata

allocate(xx(ndata),yy(ndata),dydx(ndata))

do n=1,ndata
	read(2,*) xk,yk
	xx(n)=xk
	yy(n)=yk
enddo
close(2)

dx=xx(2)-xx(1)
dydx=0.0

imin=2 ; imax=ndata-1

!INTERNOS
do i=imin,imax
dydx(i)=(yy(i+1)-yy(i-1))/(2.*dx)
enddo

!x=imin-1
dydx(imin-1)=(yy(imin)-yy(imin-1))/dx

!x=imax+1
dydx(imax+1)=(yy(imax+1)-yy(imax))/dx

open(2,file='dydx.dat')
do i=imin-1,imax+1
write(2,*) xx(i),yy(i),dydx(i)
enddo

do i=imin-1,imax+1 
print*,'Error absoluto from index',i, abs(dydx(i)-(-1.0/(xx(i)*xx(i))))
enddo

deallocate(xx,yy)
end program
