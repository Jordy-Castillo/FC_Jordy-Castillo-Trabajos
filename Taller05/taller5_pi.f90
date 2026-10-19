! Taller 5 -- Convergencia del estimador de pi por acierto y error.
! Salidas: convergencia.dat (N, media, desviación) y pendiente.dat (pendiente m).

!-----------------------------------------------------------------------------------------------!

program taller5_pi
  implicit none
  integer, parameter :: dp = kind(1.0d0)                       !precisión de punto flotante
  integer, parameter :: M = 20                                 !corridas por cada N
  integer :: k, j, n_sem, u_dat, n_puntos                      !contadores
  integer(8) :: i, n_total, n_aciertos                         !contadores
  integer, allocatable :: semilla(:)                           !semillas de random_seed
  real(dp) :: u, v, x, y, pi_est, suma, suma2, media, desv     !variables generadas
  real(dp) :: X_k, Y_k, sx, sy, sxx, sxy, pendiente            !variables minimos cuadrado

!-----------------------------------------------------------------------------------------------!

  call random_seed(size=n_sem)                                              !configurar generador
  allocate(semilla(n_sem))                                                  !reservar una "memoria" para las semillas

  sx = 0.0_dp; sy = 0.0_dp; sxx = 0.0_dp; sxy = 0.0_dp; n_puntos = 0        !sumas del ajuste
  open(newunit=u_dat,file='convergencia.dat',status='replace')

!-----------------------------------------------------------------------------------------------!

  do k = 2, 7                                                  !N = 10^2, ..., 10^7 (crear ciclo)
     n_total = 10_8**k                                         !10_8 asegura que sea entero
     suma = 0.0_dp; suma2 = 0.0_dp                             !sumas sobre las M corridas (reiniciando)
     do j = 1, M
        semilla = 1000*k + j                                   !una semilla distinta por corrida
        call random_seed(put=semilla)

!-----------------------------------------------------------------------------------------------!
! TODO A1: obtenga pi_est con n_total puntos uniformes en el cuadrado (sección 2)

        n_aciertos = 0_8                                         !acumulador
        do i = 1_8, n_total
          call random_number(u)                                  !u en [0,1)
          call random_number(v)                                  !v en [0,1)
          x = 2.0_dp*u - 1.0_dp                                  !cambio de intervalo a [-1,1]
          y = 2.0_dp*v - 1.0_dp                                  !cambio de intervalo a [-1,1]
          if (x*x + y*y < 1.0_dp) then
            n_aciertos = n_aciertos + 1_8                        !este es el punto dentro del circulo
          end if
        end do

        pi_est = 4.0_dp* real(n_aciertos, dp) / real(n_total, dp)   !estimacion de pi
        !al convertir n_aciertos y n_total a real(dp), se evita la division entera

!-----------------------------------------------------------------------------------------------!
! TODO A2: acumule lo necesario para la media y la desviación de las M corridas.

        suma = suma + pi_est                        !acumular pi_est
        suma2 = suma2 + pi_est*pi_est               !acumular pi_est^2
     end do

!-----------------------------------------------------------------------------------------------!
! TODO A3: calcule la media y la desviación estándar muestral de las M estimaciones.

     media = suma/real(M,dp)                                                                  !calculamos la media
     desv = sqrt((real(M, dp) / (real(M, dp) - 1.0_dp)) * (suma2 / real(M, dp) - media**2))   !se usa la ecuacion del enunciado

     write(u_dat,'(i10,2es25.16e3)') n_total,media,desv          !guardar archivos
     print '(i10,2f12.6)', n_total,media,desv                    !imprimir en pantalla

!-----------------------------------------------------------------------------------------------!
! TODO A4: acumule lo necesario para ajustar log10(desv) contra log10(N).

     X_k = log10(real(n_total,dp))
     Y_k = log10(desv)
     sx = sx + X_k
     sy = sy + Y_k
     sxx = sxx + X_k*X_k
     sxy = sxy + X_k*Y_k
     n_puntos = n_puntos + 1

  end do
  close(u_dat)

!-----------------------------------------------------------------------------------------------!
 ! TODO A5: calcule la pendiente del ajuste por mínimos cuadrados.

  pendiente = (real(n_puntos, dp)*sxy - sx*sy)/(real(n_puntos, dp)*sxx - sx*sx)   !ajuste por minimos cuadrados

  open(newunit=u_dat,file='pendiente.dat',status='replace')                       !escribimos pendiente en archivo
  write(u_dat,'(es25.16e3)') pendiente
  close(u_dat)

  print '(a,f9.4)', 'Pendiente de log10(s) contra log10(N): ', pendiente          !imprimimos
end program taller5_pi
