program integracion
implicit none

real :: xk,yk,h,Itrap,Isimp
integer :: reason,ndata,n,i,j,k
real, dimension(:), allocatable :: xx,yy
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

allocate(xx(ndata),yy(ndata))

do n=1,ndata
	read(2,*) xk,yk
	xx(n)=xk
	yy(n)=yk
enddo
close(2)

h=xx(2)-xx(1)

Isimp=0.0 ; Itrap=0.0

do i=1,ndata-1
Itrap=Itrap+0.5*h*(yy(i+1)+yy(i))
enddo

do i=2,ndata-1,2
Isimp=Isimp+(h/3.0)*(yy(i-1)+4.0*yy(i)+yy(i+1))
enddo

print*, 'Error from Trapezoid rule', abs(Itrap-log(abs(xx(ndata)/xx(1))))
print*, 'Error from Simpson 1/3 rule',abs(Isimp-log(abs(xx(ndata)/xx(1))))
deallocate(xx,yy)
end program
