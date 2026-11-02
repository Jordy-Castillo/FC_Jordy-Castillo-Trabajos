!!!!!!!!!! CONSTRUIR FUNCION DE INTEGRACION DE ROMBERG Y LOS INTEGRANDOS !!!!!!!!!!

module modulo_romberg
  implicit none
  integer, parameter :: dp = kind(1.0d0)                          !precisión

!! INTEGRAMOS LA POTENCIA DISIPADA POR FRICCION !!

contains
  real(dp) function f_energia(t)                                  ! Funcion integrando 2gamma*v(t)^2
    real(dp), intent(in) :: t                                     ! Tiempo (s)
    real(dp), parameter  :: gamma = 0.5_dp                        ! Coeficiente de amortiguamiento (s^-1)
    real(dp), parameter  :: w = 2.0_dp *4.0_dp* atan(1.0_dp)      ! Frecuencia angular de oscilación (rad/s) con PI = 4*arctan(1)
    real(dp) :: v                                                 ! Velocidad (m/s)

    v = -exp(-gamma*t) * (gamma*cos(w*t) + w*sin(w*t))            !Calculo de la Velocidad en funcion del tiempo v(t)
    f_energia = 2.0_dp * gamma * v**2                             !Calculo de la energia
  end function f_energia

!! INTEGRAMOS LA FUNCION PARA EL CALCULO PI > 4/(1+X^2) !!

  real(dp) function f_pi(x)                                       ! Funcion integrando 4/(1+x^2)
    real(dp), intent(in) :: x                                     ! Posicion
    f_pi = 4.0_dp / (1.0_dp + x**2)                               !Calculo de la funcion
  end function f_pi

!! CREAMOS UNA FUNCION QUE SELECCIONE EL INTENGRANDO !!

  real(dp) function eval_f(x, tipo_funcion)                       !funcion integrando 4/(1+x^2) --> es como un selector
    real(dp), intent(in) :: x                                     ! Posicion
    integer, intent(in) :: tipo_funcion                          ! Tipo de funcion (energia o de pi)

    if (tipo_funcion == 1) then                                   !Si el tipo de funcion es 1 llama al integrando 2gamma*v(t)^2
      eval_f = f_energia(x)
    else if (tipo_funcion == 2) then                              !Si el tipo de funcion es 2 llama al integrando 4/(1+x^2)
      eval_f = f_pi(x)
    else
      eval_f = 0.0_dp                                             !Si el tipo de funcion no es valido devuelve 0
    end if
  end function eval_f

!! SUBRUTINA PARA LA INTEGRACION DE ROMBERG !!

  subroutine integrar_romberg(a, b, tol, kmax, tipo_funcion, S, k_final, n_eval)
    real(dp), intent(in)  :: a, b, tol                            !Rango de integracion (a y b), tolerancia
    integer, intent(in)   :: kmax, tipo_funcion                    !Numero de niveles y tipo de funcion
    real(dp), intent(out) :: s(0:, 0:)                            !Arreglo de resultados
    integer, intent(out)  :: k_final, n_eval                      !Nivel final y numero de evaluaciones
    integer  :: k, j, i, puntos_nuevos                             !Variables para el bucle
    real(dp) :: h, suma, x_nuevo                                  !Variables para el bucle

    !INICIAR EL NIVEL K = 0
    h = b - a
    s(0,0) = 0.5_dp * h * (eval_f(a, tipo_funcion) + eval_f(b, tipo_funcion))    !regla del trapecio cuando K=0
    n_eval = 2                                                                   !inicia un contador en 2 evaluaciones (f(a) y f(b))
    k_final = kmax                                                               !inicia el nivel maximo en el nivel maximo

    !BUCLE PARA EL NIVEL K = 1, 2, ..., KMAX
    do k = 1, kmax                                                               !realiza los refinamientos de la malla
      h = h/2.0_dp                                                               !se reduce el ancho del intervalo
      suma = 0.0_dp                                                              !se reinicia la suma
      puntos_nuevos = 2**(k-1)                                                   !se calcula el numero de puntos nuevos en el nivel k

      do i = 1, puntos_nuevos
        x_nuevo = a + real(2*i-1, dp) * h                                        !puntos medios impares del intervalo
        suma = suma + eval_f(x_nuevo, tipo_funcion)                              !se suma la funcion en los puntos medios impares
      end do

      n_eval = n_eval + puntos_nuevos                                            !incrementa el contador de evaluaciones
      S(k,0) = 0.5_dp * S(k-1,0) + h * suma                                      !aplicar formula de refinamiento del trapecio

      do j = 1, k
        S(k,j) = ((4.0_dp**j) * S(k, j-1) - S(k-1, j-1)) / ((4.0_dp**j) - 1.0_dp)  !extrapolacion de richardson para las columnas superiores
      end do

      if (abs(S(k,k) - S(k-1, k-1)) < tol ) then                                 !criterio de parada
        k_final = k                                                              !se actualiza el nivel final
        exit                                                                     !salir del bucle
      end if
    end do
  end subroutine integrar_romberg
end module modulo_romberg
