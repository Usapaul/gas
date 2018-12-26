module runmethods

use init

implicit none

contains



function fill_pressure(W,U,V,Rho) result(Pr)
	! Процедура заполняет массив Pr (давление), используя значения
	! в массивах W, U, V, Rho (энергия, компоненты скорости и плотность).
	! Используется уравнение состояния для идеального одноатомного газа:
	! p = rho * eps * (gamma - 1), gamma = 5/3
	implicit none

	real(pl), dimension(0:,0:,1:), intent(in) :: W, U, V, Rho
	real(pl), dimension(0:size(W,1)-1,0:size(W,2)-1,1:size(W,3)) :: Pr

	integer :: N

	!--------------------------------------------
	N = size(W,1) - 1

	! ур-е сост-я ideal_1atom_gas: p = rho * eps * (gamma - 1), gamma = 5/3
	! eps (удельная энергия) = w - (u**2 + v**2)/2
	Pr(1:N,1:N,1) = (W(1:N,1:N,1) - (U(1:N,1:N,1)**2 + V(1:N,1:N,1)**2)/2) * &
					& Rho(1:N,1:N,1) * 2._pl/3 ! (gamma-1 = 5/3-1 = 2/3)

	call fill_layers(Pr)

end function fill_pressure



subroutine fill_layers(A)
	! Процедура из данных в массиве A с целочисленными индексами (1-ый слой)
	! заполняет остальные три слоя (значения на гранях ячеек и в углах),
	! используя линейную аппроксимацию -- среднее значение смежных ячеек
	implicit none

	real(pl), dimension(0:,0:,1:), intent(inout) :: A

	integer :: N

	!--------------------------------------------
	N = size(A,1) - 1

	! Значение в полушаге сетки (грань ячейки) -- среднее между
	! значениями в ячейках. Для угловых -- среднее по диагонали
	A(1:N-1,1:N,2) = (A(2:N,1:N,1) + A(1:N-1,1:N,1)) / 2
	A(1:N,1:N-1,3) = (A(1:N,2:N,1) + A(1:N,1:N-1,1)) / 2
	A(1:N-1,1:N-1,4) = (A(2:N,2:N,1) + A(1:N-1,1:N-1,1)) / 2
	
	!--------------------------------------------
	! А теперь пытаюсь разобраться, чем заполнять значения на краях...

	! В середине верхних и нижних границ краевых ячеек:
	A(0,1:N,2) = A(1,1:N,1)
	A(N,1:N,2) = A(N,1:N,1)
	
	! В середине боковых границ краевых ячеек:
	A(1:N,0,3) = A(1:N,1,1)
	A(1:N,N,3) = A(1:N,N,1)

	! В углах краевых ячеек:
	A(0,0,4) = A(1,1,1)
	A(N,0,4) = A(N,1,1)
	A(N,N,4) = A(N,N,1)
	A(0,N,4) = A(1,N,1)

	A(0,1:N-1,4) = (A(1,2:N,1) + A(1,1:N-1,1)) / 2
	A(N,1:N-1,4) = (A(N,2:N,1) + A(N,1:N-1,1)) / 2
	A(1:N-1,0,4) = (A(2:N,1,1) + A(1:N-1,1,1)) / 2
	A(1:N-1,N,4) = (A(2:N,N,1) + A(1:N-1,N,1)) / 2	

end subroutine fill_layers



subroutine step1(Pr,U,V,Rho,W)
	! Процедура выполняет первый "полушаг" -- по методу PIC
	implicit none

	real(pl), dimension(0:,0:,1:), intent(inout) :: Pr, U, V, Rho, W

	integer :: i, j
	integer :: N

	!--------------------------------------------
	N = size(Pr,1) - 1

	! Rho остается без изменений: 1/tau * (rho(i,j) - rho(i,j)^n) = 0 

	call fill_layers(U)
	call fill_layers(V)
	call fill_layers(Rho)

	! По формулам из L07.pdf вычисляются U, V и W:
	forall (i = 1 : N, j = 1 : N)
		U(i,j,1) = U(i,j,1) - tau/h * (Pr(i,j,2) - Pr(i-1,j,2)) / Rho(i,j,1)
		V(i,j,1) = V(i,j,1) - tau/h * (Pr(i,j,3) - Pr(i,j-1,3)) / Rho(i,j,1)
		W(i,j,1) = W(i,j,1) - tau/h / Rho(i,j,1) * &
			& ( Pr(i,j,2) * U(i,j,2) - Pr(i-1,j,2) * U(i-1,j,2) + &
			&   Pr(i,j,3) * V(i,j,3) - Pr(i,j-1,3) * V(i,j-1,3)     ) 
	end forall

end subroutine step1



subroutine step2(Pr,U,V,Rho,W)
	! Процедура выполняет второй "полушаг" -- по методу PIC
	implicit none

	real(pl), dimension(0:,0:,1:), intent(inout) :: Pr, U, V, Rho, W

	! В логическом массиве sgn будут записаны знаки скоростей (см. L07.pdf).
	! .TRUE. -- если знак "+", .FALSE. -- если знак "-"
	logical, dimension(0:size(Pr,1)-1,0:size(Pr,2)-1,1:size(Pr,3)) :: sgn

	real(pl), dimension(0:size(Pr,1),0:size(Pr,1),2:3) :: PP ! Потоки "П"

	! Понадобится сохранить плотность на n-ном шаге по времени, это есть
	! плотность с тильдой. Сохраняю ее в массив Rho0
	real(pl), dimension(0:size(Rho,1),0:size(Rho,2),1:size(Rho,3)) :: Rho0
	real(pl), dimension(0:size(Pr,1),0:size(Pr,1),1:size(Pr,3)) :: F ! см. ниже


	integer :: i, j
	integer :: N

	!--------------------------------------------
	N = size(Pr,1) - 1
	sgn = .FALSE.
	Rho0 = Rho

	forall (i = 1 : N-1, j = 1 : N)
		sgn(i,j,2) = U(i,j,1) + U(i+1,j,1) >= 0
	end forall

	forall (i = 1 : N, j = 1 : N-1)
		sgn(i,j,3) = V(i,j,1) + V(i,j+1,1) >= 0
	end forall

	!--------------------------------------------
	forall (i = 2 : N-2, j = 2 : N-2)
		U(i,j,2) = velocity(U(i-1:i+2,j,1),sgn(i,j,2))
		V(i,j,3) = velocity(V(i,j-1:j+2,1),sgn(i,j,3))
	end forall

	forall (i = 0 : 1) U(i,1:N,2) = U(2,1:N,2)
	forall (i = N-1 : N) U(i,1:N,2) = U(N-2,1:N,2)
	forall (j = 0 : 1) V(1:N,j,3) = V(1:n,2,3)
	forall (j = N-1 : N) V(1:N,j,3) = V(1:N,N-2,3)
	
	!--------------------------------------------
	! Потокам присваиваются значения (см. формулу ниже) только
	! если скорости (рассчитанные в предыдущем блоке данной процедуры)
	! имеют такой же знак, какой был сохранен в массив sgn.
	! Если не такой же, то считается, что поток равен нулю
	where (sgn(:,:,2) .neqv. (U(:,:,2) >= 0))
		PP(:,:,2) = h * tau * U(:,:,2) * Rho(:,:,2)
	else where
		PP(:,:,2) = 0
	end where	

	where (sgn(:,:,3) .neqv. (V(:,:,3) >= 0))
		PP(:,:,3) = h * tau * V(:,:,3) * Rho(:,:,3)
	else where
		PP(:,:,3) = 0
	end where

	!--------------------------------------------
	! Из аппроксимации уравнения неразрывности получаем новую
	! плотность (с индексом шага по времени равным n+1). Но старую
	! (с тильдой) сохранили уже в Rho0 -- пригодится далее
	forall (i = 1 : N, j = 1 : N)
		Rho(i,j,1) = Rho0(i,j,1) - 1 / h**2 * &
			& (PP(i,j,2) - PP(i-1,j,2) + PP(i,j,3) - PP(i,j-1,3))
	end forall

	!--------------------------------------------
	! Теперь аппроксимация оставшихся трех уравнений переноса, где
	! по условию за F у нас обозначено {u,v,w}. Из неё получаем
	! новые значения F
	call sub_for_F(U,PP,sgn,Rho,Rho0)
	call sub_for_F(V,PP,sgn,Rho,Rho0)
	call sub_for_F(W,PP,sgn,Rho,Rho0)

end subroutine step2



pure real(pl) function velocity(v,sgn)
	! Функция вычисляет скорость V(i+1/2) по имеющимся данным
	! V(i-1), V(i), V(i+1), V(i+2) и знаку sgn (см. L07.pdf)
	implicit none
	
	real(pl), dimension(-1:2), intent(in) :: v
	logical, intent(in) :: sgn

	!--------------------------------------------
	if (sgn) then
		! sign = 1
		velocity = v(0) + (v(1) - v(-1)) / 4
	else
		! sign = -1
		velocity = v(1) - (v(2) - v(0)) / 4
	end if

end function velocity



subroutine sub_for_F(F,PP,sgn,Rho,Rho0)
	! Процедура, согласно формулам в L07.pdf, получает F^n+1 из F-тильда.
	! Для этого нужно знать значения плотностей с тильдой и Rho^n+1,
	! а также массив PP со значениями для потоков и массив со знаками sgn
	implicit none

	real(pl), dimension(0:,0:,1:), intent(inout) :: F
	real(pl), dimension(0:,0:,1:), intent(in) :: PP, Rho, Rho0
	logical, dimension(0:,0:,1:), intent(in) :: sgn

	integer :: i, j, N
	!--------------------------------------------
	N = size(F,1) - 1

	where (sgn(1:N-1,1:N,2))
		F(1:N-1,1:N,2) = F(1:N-1,1:N,1)
	else where
		F(1:N-1,1:N,2) = F(2:N,1:N,1)
	end where
	F(N,1:N,2) = F(N-1,1:N,2)
	F(0,1:N,2) = F(1,1:N,2)

	where (sgn(1:N,1:N-1,3))
		F(1:N,1:N-1,3) = F(1:N,1:N-1,1)
	else where
		F(1:N,1:N-1,3) = F(1:N,2:N,1)
	end where
	F(1:N,N,3) = F(1:N,N-1,3)
	F(1:N,0,3) = F(1:N,1,3)

	!--------------------------------------------
	forall (i = 1 : N, j = 1 : N)
		F(i,j,1) = 1 / Rho(i,j,1) * &
			& ( Rho0(i,j,1) * F(i,j,1) - 1 / h**2 * &
			&  (PP(i,j,2) * F(i,j,2) - PP(i-1,j,2) * F(i-1,j,2) + &
			&   PP(i,j,3) * F(i,j,3) - PP(i,j-1,3) * F(i,j-1,3)  ) )
	end forall

end subroutine sub_for_F



end module runmethods