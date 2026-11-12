!!!!!!!!!! AÑADIR NUEVOS INTEGRANDOS AL MODULO !!!!!!!!!!

module modulo_romberg
  implicit none
  integer, parameter :: dp = kind(1.0d0)                          !definir precisión

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! INTEGRADORES NUEVOS PARA LA PARTE C (NUEVA PARTE) !!

  !!!! CON CAMBIO DE VARIABLE EN [0,1] !!!!
  real(dp) function f_sigma_t(t)                                  !x(t) = t / (1 - t), dx = dt / (1 - t)^2
    real(dp), intent(in) :: t
    real(dp) :: x, exp_val                                        !Varaibles

    if (t <= 0.0_dp .or. t >= 1.0_dp) then                        !controlar extremos para evitar indeterminaciones
       f_sigma_t = 0.0_dp                                         !no hay cambio de variable
       return
    end if

    x = t / (1.0_dp - t)                                          !definir la variable para su cambio de variable

    if (x > 700.0_dp) then                                        !proteger ante un desbordamiento si (x) es muy grande
       f_sigma_t = 0.0_dp                                         !no hay cambio de variable
       return
    end if

    exp_val = exp(x)                                              !re-definimos la exponencial
    f_sigma_t = (x**3 / (exp_val - 1.0_dp)) / ((1.0_dp - t)**2)   !definimos la funcion original (I)
  end function f_sigma_t

  !!!! SIN CAMBIO DE VARIABLE EN [0,100] !!!!
  real(dp) function f_sigma_x(x)                                  !f(x) = x^3 / (e^x - 1)
    real(dp), intent(in) :: x

    if (x <= 0.0_dp) then                                         !controlar extremos para evitar indeterminaciones
       f_sigma_x = 0.0_dp                                         !no hay cambio de variable
       return
    end if

    if (x > 700.0_dp) then                                        !proteger ante un desbordamiento si (x) es muy grande
       f_sigma_x = 0.0_dp                                         !no hay cambio de variable
       return
    end if

    f_sigma_x = (x**3) / (exp(x) - 1.0_dp)                        !definimos la funcion original
  end function f_sigma_x


  !!!! SELECTOR DE INTEGRANDOS PARA ESTE TALLER !!!!
  real(dp) function eval_f(x, tipo_funcion)                       !creamos una funcion evaluadora
    real(dp), intent(in) :: x                                     !variable
    integer, intent(in) :: tipo_funcion                           !tipo de funcion

    select case (tipo_funcion)
    case (1)
       eval_f = f_energia(x)                                      !funcion de energia con cambio de variable [0, 1]
    case (2)
       eval_f = f_pi(x)                                           !funcion de romberg sin cambio de variable [0, 100]
    case (3)
       eval_f = f_sigma_t(x)                                      !funcion con cambio de variable [0, 1]
    case (4)
       eval_f = f_sigma_x(x)                                      !funcion sin cambio de variable [0, 100]
    case default
       eval_f = 0.0_dp                                            !funcion hay indeterminacion o no es valida
    end select
  end function eval_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
