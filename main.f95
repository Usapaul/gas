program gas

use init
use runmethods

implicit none

! pl -- целое число отвечает за параметр точности real
! N -- количество ячеек (всё поле -- квадрат NxN) -- в модуле init

integer :: i, j, k

character(15) :: filename, kstr

!------------------------------------------------

call init_const() ! Просто инициализация начальных значений, например, N

call init_values() ! Присвоение нач. знач. в момент времени = 0 (Pr0, Rho0,..)

filename = 'time0.dat'
! open(100,file=filename,status='replace')
! do i = 1, N
! 	do j = 1, N
! 		write(100,*) i, j, W(i,j,1)
! 	end do
! end do
! close(100,status='keep')

open(100,file=filename,status='replace')
do i = 1, N
	write(100,*) W(1:N,i,1)
end do
close(100,status='keep')


! Цикл, в котором каждая итерация получает значения для уровня k+1 по времени 
! из имеющихся значений для k-того уровня:
do k = 1, maxTimelvl
	write(*,*) 'NUMBER: ', k
	Pr = fill_pressure(W,U,V,Rho)
	! Один шаг разбивается на два полушага step1 и step2 (модуль runmethods):
	call step1(Pr,U,V,Rho,W)
	call step2(Pr,U,V,Rho,W)

	if (mod(k,maxTimelvl/10) == 0) then
		write(kstr,'(i5)') k
		filename = 'time'//trim(adjustl(kstr))//'.dat'
		open(100,file=filename,status='replace')
		do i = 1, N
			write(100,*) W(1:N,i,1)
		end do
		close(100,status='keep')
	end if
end do


end program gas
