!!!!!!!!!! CALCULAR ENERGIA DISIPADA EN UN OSCILADOR AMORTIGUADO !!!!!!!!!!

program programa_energia
  use modulo_romberg, only : dp, integrar_romberg
  implicit none

  real(dp), parameter :: gamma = 0.5_dp                             ! Coeficiente de amortiguamiento (s^-1)
  real(dp), parameter :: w = 2.0_dp * 4.0_dp * atan(1.0_dp)         ! Frecuencia angular
  real(dp), parameter :: w0_sq = w**2 + gamma**2                    ! Frecuencia natural al cuadrado
  real(dp), parameter :: a = 0.0_dp                                 ! Tiempo inicial (s)
  real(dp), parameter :: b = 3.0_dp                                 ! Tiempo final (s)
  real(dp), parameter :: tol = 1.0e-6_dp                            ! Tolerancia corregida
  integer, parameter :: kmax = 12                                   ! nivel maximo
  real(dp), allocatable :: S(:,:)                                   ! Arreglo de resultados
  integer  :: k_final, n_eval, k, j                                 ! Variables para el bucle
  real(dp) :: e_dis, e_0, e_T, e_ref, diferencia                    ! Variables para el bucle

  allocate(S(0:kmax, 0:kmax))                                       ! Reservar memoria para el arreglo

!! CALCULO DE E_DIS Y COMPARACION CON E_0-E_T DE B.3 Y B.4!!

  call integrar_romberg(a, b, tol, kmax, 1, S, k_final, n_eval)     !llamamos a la funcion integrar_romberg para calcular la energia disipada
  e_dis = S(k_final, k_final)                                       !Energia disipada de la punta diagonal principal
  e_0 = 0.5_dp * ((-gamma)**2) + 0.5_dp * w0_sq * (1.0_dp**2)       !Energia del oscilador al tiempo inicial, t=0 : x(0)= 1.0 and v(0) = -gamma
  e_T = calcular_energia_exacta(b, gamma, w, w0_sq)                 !Energia del oscilador al tiempo final

  e_ref = e_0 - e_T                                                 ! Energia exacta
  diferencia = abs(e_dis - e_ref)                                   ! Diferencia de energia exacta

!! ELABORACION DE TABLA_ENERGIA.DAT DE B.5!!

  open(unit = 10, file = "tabla_energia.dat", status = "replace", action = "write")
  do k = 0, k_final
    write(10, "(I3)", advance = "NO") k
    do j = 0, k
      write(10, "(1X, ES25.15)", advance = "NO") S(k,j)             !advance="NO" permite escribir varios numeros en la misma fila (NUEVO)
    end do
    write(10, *)
  end do
  close(10)

!! ESCRIBIR ARCHIVO RESULTADO_ENERGIA.DAT !!

  open(unit=20, file="resultado_energia.dat", status="replace", action="write")
  write(20, "(I5, 1X, I10, 1X, ES25.15, 1X, ES25.15, 1X, ES25.15)") k_final, n_eval, e_dis, e_ref, diferencia
  close(20)

  deallocate(S)                                                          !Liberar memoria del arreglo

contains
  real(dp) function calcular_energia_exacta(t, gamma, w, w0_sq)          !funcion para la energia exacta en un tiempo = t
    real(dp), intent(in) :: t, gamma, w, w0_sq                           !Variables de entrada
    real(dp) :: xt, vt                                                   !variables para la energia exacta
    xt = exp(-gamma*t) * cos(w*t)                                        !calcular la posicion exacta con la ecuacion
    vt = -exp(-gamma*t) * (gamma*cos(w*t) + w*sin(w*t))                  !calcular la velocidad exacta con la ecuacion
    calcular_energia_exacta = 0.5_dp * vt**2 + 0.5_dp * w0_sq * xt**2    !calcular la energia exacta

  end function calcular_energia_exacta

end program programa_energia
