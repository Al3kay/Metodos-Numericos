program lagrange
implicit none

real :: xk,yk,p,xp,yp
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


xp=3.44
yp=0.0
do i=1,ndata
p=1.0
	do j=1,ndata
	if(i.ne.j) then
		p=p*(xp-xx(j))/(xx(i)-xx(j))
	end if
	enddo
yp=yp+p*yy(i)
enddo

print*,''
print*,'xp=',xp,'yp=',yp
print*,''
print*,'Error absoluto', ABS(1.0/xp - yp)
deallocate(xx,yy)
end program
