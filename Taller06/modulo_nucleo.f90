!!!!!!!!!! Módulo de la fórmula semiempírica !!!!!!!!!! 
 

!!! DECLARAMOS LAS CONSTANTES !!!

module modulo_nucleo
  implicit none
  integer, parameter :: dp = kind(1.0d0)      !doble precisión
  real(dp), parameter :: a_V = 15.8_dp        !coeficiente de volumen
  real(dp), parameter :: a_S = 18.3_dp        !coeficiente de superficie 
  real(dp), parameter :: a_C = 0.71_dp        !coeficiente de Coulomb
  real(dp), parameter :: a_A = 92.7_dp        !Coeficiente de asimetria
  real(dp), parameter :: M_n = 939.57_dp      !masa del neutrón
  real(dp), parameter :: M_H = 938.78_dp      !masa del atomo de hidrogeno (H)


!!! DECLARAMOS LAS 4 FUNCIONES DE LOS TERMINOS INDIVIDUALES !!!

contains

  !!! Función del Volumen !!!
  real(dp) function termino_vol(a)
    implicit none
    integer, intent(in) :: a
    termino_vol = -a_V + real(0*a, dp)          !Evita la advertencia de variable no usada (fallo continuo que tuve)
  end function termino_vol

  !!! Función de la Superficie !!!
  real(dp) function termino_sup(a)
    implicit none
    integer, intent(in) :: a
    termino_sup = a_S * (real(a, dp)**(-1.0_dp/3.0_dp))
  end function termino_sup

  !!! Función del Coulomb !!!
  real(dp) function termino_coul(z, a)
    implicit none
    integer, intent(in) :: z, a
    termino_coul = a_C * (real(z, dp)**2.0_dp) * (real(a, dp)**(-4.0_dp/3.0_dp))
  end function termino_coul

  !!! Función de la Asimetria !!!
  real(dp) function termino_asim(z, a)
    implicit none
    integer, intent(in) :: z, a
    termino_asim = (a_A * ((real(z, dp) - real(a, dp)/2.0_dp)**2.0_dp) / real(a, dp)**2.0_dp)
  end function termino_asim


!!! DECLARAMOS LA FUNCION QUE SUMA LOS 4 TERMINOS !!!

  real(dp) function E_B(z, a)
    implicit none
    integer, intent(in) :: z, a
    E_B = termino_vol(a) + termino_sup(a) + termino_coul(z, a) + termino_asim(z, a)
  end function E_B

!!! DECLARAMOS LA FUNCIÓN QUE EVALUA LA MASA ATOMICA !!!

  real(dp) function M_masa(z, a)
    implicit none
    integer, intent(in) :: z, a
    M_masa = real(z, dp) * M_H + (real(a, dp) - real(z, dp))*M_n + real(a, dp)*E_B(z, a)
  end function M_masa

!! CREAMOS UNA SUBRUTINA PARA ENCONTRAR EL Z DE MENOR E_B !!!

  subroutine z_estable(a, z_optimo, EB_min)
    implicit none
    integer, intent(in) :: a                     !variable de entrada
    integer, intent(out) :: z_optimo             !variable de salida
    real(dp), intent(out) :: EB_min              !variable de salida

    integer :: z                                 !variable local
    real(dp) :: EB_actual                        !variable local
    z_optimo = 1                                 !iniciar con un valor arbitrario
    EB_min = E_B(1, a)                           !iniciar con un valor arbitrario
   
    do z = 2, a                                  !creamos un ciclo en Z = 2,3,4... A
      EB_actual = E_B(z, a)
      if (EB_actual < EB_min) then               !si EB_actual es menor que EB_min entonces;
        EB_min = EB_actual                       !actualizamos el EB_min
        z_optimo = z                             !actualizamos el Z_optimo
      endif
    end do
  end subroutine z_estable

  end module modulo_nucleo
