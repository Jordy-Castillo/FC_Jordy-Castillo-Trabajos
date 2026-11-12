!!!!!!!!! DEFINIR LA FUNCION Y SUS METODOS !!!!!!!!!

module modulo_ceros
  implicit none
  integer, parameter :: dp = kind(1.0d0)                                        !definir precisión

!! FUNCION ORIGINAL !!
contains
  real(dp) function F(x, n)
    real(dp), intent(in) :: x                                                   !raíz de la función
    real(dp), intent(in) :: n                                                   !número de la función
    F = x - n * (1.0_dp - exp(-x))                                              !definimos la funcion que busca el cero
  end function F

!! FUNCION DERIVADA !!

  real(dp) function DF_F(x,n)
    real(dp), intent(in) :: x                                                   !raíz de la función
    real(dp), intent(in) :: n                                                   !número de la función
    DF_F = 1.0_dp - n * exp(-x)                                                 !definimos la derivada de la funcion (calculada a mano)
  end function DF_F

!! METODO BISECCION !!

  subroutine biseccion(n, a_int, b_int, tol, max_inter, raiz, iteracciones, u_hist)
    real(dp), intent(in)  :: n, a_int, b_int, tol                  !definir los parametros de entrada
    integer, intent(in)   :: max_inter                             !definir el maximo de iteraciones
    real(dp), intent(out) :: raiz                                  !definir la raiz que usaremos
    integer, intent(out)  :: iteracciones                          !definir el numero de iteraciones
    integer, intent(in), optional :: u_hist                        !definimos el "historial" de iteracciones con "Optional" para que no sea obligatorio
    real(dp) :: a, b, c                                            !variables para proteger los datos originales
    integer :: i                                                   !contador arbitrario

    a = a_int                                                                   !#re-definimos los datos de entrada para poder modificarlos
    b = b_int                                                                   !#recordando que las variables de lectura no se pueden modificar
    iteracciones = 0

    if (F(a, n) * F(b, n) >= 0.0_dp) then                                       !verificar el cambio del signo incial
       stop "Bisección: f(a) y f(b) no tienen signos distintos"
    end if

    do i = 1, max_inter                                                         !punto medio del intervalo
      c = 0.5_dp * (a + b)                                                      !conteo de las iteracciones
      iteracciones = i

      if (present(u_hist)) then                                                 !##"present" indetifica si u_hist fue dado o no por el programa
         write(u_hist, '(i2, i4, es24.16)') 1, i, c                             !##lo hacemos para resutilizar el codigo mas adelante
      end if

      if ((b - c) <= tol) then                                                  !criterio de parada
        raiz = c                                                                !guardar la raiz
        return                                                                  !terminar la subrutina
      end if

      if (F(a,n) * F(c,n) <= 0.0_dp) then                                       !actualizar el intervalo de la cota de error
        b = c
      else
        a = c
      end if
    end do

    raiz = c                                                                    !si se llega aqui, no se cumplio el criterio de parada
    print *, "Bisección: No convergió en max_inter iteraciones"
  end subroutine biseccion


!! METODO NEWTON-RAPHSON !!

  subroutine newton(n, x_int, tol, max_inter, raiz, iteracciones, u_hist)
    real(dp), intent(in)  :: n, x_int, tol                          !definir los parametros de entrada
    integer, intent(in)   :: max_inter                              !definir el maximo de iteraciones
    real(dp), intent(out) :: raiz                                   !definir la raiz que usaremos
    integer, intent(out)  :: iteracciones                           !definir el numero de iteracciones
    integer, intent(in), optional :: u_hist                         !definimos el "historial" de iteracciones con "Optional" para que no sea obligatorio
    real(dp) :: x_cu, x_nex, paso, dfx                              !variables para proteger los datos originales
    integer :: i                                                    !contador arbitrario

    x_cu = x_int                                                    !#re-definimos los datos de entrada para poder modificarlos
    x_nex = x_cu                                                    !#recordando que las variables de lectura no se pueden modificar
    iteracciones = 0

    do i = 1, max_inter                                             !contador arbitrario comienzo con ciclo do
      dfx = DF_F(x_cu,n)                                            !definimos la funcion derivada
      if (abs(dfx) < tiny(1.0_dp)) then
         stop "Newton-Raphson: Derivada Nula o muy cercana a cero"
      end if

      paso = -f(x_cu,n) / dfx                                        !calculo del paso con correcion de Newton
      x_nex = x_cu + paso                                            !calculo de la nueva aproximacion
      iteracciones = i                                               !conteo de las iteracciones

      if (present(u_hist)) then                                                 !##"present" indetifica si u_hist fue dado o no por el programa
         write(u_hist, '(i2, i4, es24.16)') 2, i, x_nex
      end if

      if (abs(paso) < tol) then                                                 !criterio de parada por paso
         raiz = x_nex                                                           !actualizar la aproximacion
         return
      end if
      x_cu = x_nex
    end do

    raiz = x_nex
    print *, "Newton-Raphson; No convergió en max_inter iteraciones"
  end subroutine newton


!! METODO DE LA SECANTE !!

  subroutine secante(n, x_int, x_int_1, tol, max_inter, raiz, iteracciones, u_hist)
    real(dp), intent(in)  :: n, x_int, x_int_1, tol                !definir los parametros de entrada
    integer, intent(in)   :: max_inter                             !definir el maximo de  iteraciones
    real(dp), intent(out) :: raiz                                  !definir la raiz que usaremos
    integer, intent(out)  :: iteracciones                          !definir el numero de iteracciones
    integer, intent(in), optional :: u_hist
    real(dp) :: x_cu, x_nex, paso                                  !variables para proteger los datos originales
    real(dp) :: x_aux, f_aux, f_cu, denom                          !variables para proteger los datos originales
    integer :: i                                                   !contador arbitrario

    x_aux = x_int                                                  !#re-definimos los datos de entrada para poder modificarlos
    x_cu  = x_int_1                                                !#recordando que las variables de lectura no se pueden modificar
    x_nex = x_cu                                                   !inicializar x_nex para que nunca quede sin valor
    f_aux = F(x_aux, n)                                            !inicializar f_aux para que nunca quede sin valor
    iteracciones = 0                                               !inicar contador en 0

    do i = 1, max_inter                                            !contador arbitrario
      f_cu  = F(x_cu, n)                                           !calcular las funciones en los puntos de f_cu
      denom = f_cu - f_aux                                         !calculo del denominador

      if (abs(denom) < tiny(1.0_dp)) then
         stop "Secante: Dominador Nulo o muy cercano a cero"
      end if

      paso = -f_cu * (x_cu - x_aux) / denom
      x_nex = x_cu + paso
      iteracciones = i                                                   !En la secante, el primer valor nuevo corresponde a x_2

      if (present(u_hist)) then
         write(u_hist, '(i2, i4, es24.16)') 3, i, x_nex
      end if

      if (abs(paso) < tol) then                                                 !criterio de parada por paso
         raiz = x_nex
         return
      end if

      x_aux = x_cu                                                              !actualizar los puntos de la secante
      f_aux = f_cu
      x_cu  = x_nex
    end do

    raiz = x_nex
    print *, "Secante: No convergió en max_inter iteraciones"
  end subroutine secante

end module modulo_ceros
