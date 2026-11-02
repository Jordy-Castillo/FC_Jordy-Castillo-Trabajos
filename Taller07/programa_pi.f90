
program programa_pi
  use modulo_romberg, only: dp, integrar_romberg
  implicit none

  real(dp), parameter :: a = 0.0_dp                !Rango de integracion inferior
  real(dp), parameter :: b = 1.0_dp                !Rango de integracion superior
  real(dp), parameter :: tol = 1.0e-5_dp           !Tolerancia
  integer, parameter :: kmax = 10                  !nivel maximo
  real(dp) :: pi_exacto, pi_numerico, error_pi     !Variables para el calculo de pi
  integer :: k_final, n_eval                       !Variables para el bucle
  real(dp), allocatable :: S(:,:)                  !Arreglo de resultados

  allocate(S(0:kmax, 0:kmax))                                       !Reservar memoria para el arreglo
  call integrar_romberg(a, b, tol, kmax, 2, S, k_final, n_eval)     !llamar a la funcion integrar_romberg para calcular

!! CALCULO DEL VALOR NUMERICO Y DEL ERROR !!
  pi_exacto = 4.0_dp*atan(1.0_dp)
  pi_numerico = S(k_final, k_final)
  error_pi = abs(pi_numerico - pi_exacto)

!! ESCRIBIR ARCHIVO RESULTADO_PI.DAT DE D.2 !!
  open(unit=10, file="resultado_pi.dat", status="replace", action="write")
  write(10, "(I5, 1X, I10, 1X, ES25.15, 1X, ES25.15)") k_final, n_eval, pi_numerico, error_pi
  close(10)

  deallocate(S)                                                     !Liberar memoria del arreglo

end program programa_pi
