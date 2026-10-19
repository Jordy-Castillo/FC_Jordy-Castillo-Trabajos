! Taller 5 -- Integral de sin(x) en [0,pi/2) por valor medio y acierto y error.
! Salida: seno.dat (media_vm, desv_vm, media_ae, desv_ae).

!-----------------------------------------------------------------------------------------------!

program taller5_seno
  implicit none
  integer, parameter :: dp = kind(1.0d0)                                           !precisión de punto flotante
  integer, parameter :: N = 1000, R = 1000                                         !muestras por corrida y repeticiones
  real(dp), parameter :: a = 0.0_dp, b = 2.0_dp*atan(1.0_dp), H = 1.0_dp           !rangos de integración
  integer :: rep, i, n_bajo, n_sem, u_dat                                          !contadores a usar
  integer, allocatable :: semilla(:)                                               !semillas random
  real(dp) :: u, v, x, y, suma_f, I_vm, I_ae                                       !variables
  real(dp) :: s_vm, s2_vm, s_ae, s2_ae, media_vm, desv_vm, media_ae, desv_ae       !variables

!-----------------------------------------------------------------------------------------------!

!El tamaño del arreglo de semilla depende de la implementación.
  call random_seed(size=n_sem)                                         !configurar generador
  allocate(semilla(n_sem))                                             !conservar una "memoria" para las semillas
  semilla = 314                                                        !se elige una semilla fija
  call random_seed(put=semilla)                                        !se fija una vez, antes de las R corridas

  !acumuladores para las R corridas
  s_vm = 0.0_dp
  s2_vm = 0.0_dp
  s_ae = 0.0_dp
  s2_ae = 0.0_dp


!-----------------------------------------------------------------------------------------------!
! TODO B1: obtenga I_vm por valor medio usando N puntos (sección 4).

  do rep = 1, R
     suma_f = 0.0_dp                                                    !creamos un contado de la suma
     do i = 1, N                                                        !creamos ciclo para N
       call random_number(u)                                            !u en [0,1)
       x = a + (b-a)*u                                                  !realizamos la combinacion lineal para [0, pi/2)
       suma_f = suma_f + sin(x)                                         !acumulamos la suma
     end do
     I_vm = (b-a)*(suma_f/real(N,dp))                                   !calculamos integral con la ecuacion

!-----------------------------------------------------------------------------------------------!
! TODO B2: obtenga I_ae por acierto y error en [a,b) x [0,H) (sección 4).

     n_bajo= 0                                                          !contador de puntos debajo de la curva
     do i = 1,N
       call random_number(u)
       call random_number(v)
       x = a + (b-a)*u                                                  !combinacion lineal para x para [0, pi/2)
       y = v*H                                                          !combinacion lineal para y para [0,1)
       if (y < sin(x)) then                                             !los puntos deben cumplir esta condicion
         n_bajo = n_bajo + 1                                            !punto debajo de la curva de sin(x)
       end if
     end do

     I_ae = (b-a)*H*((real(n_bajo,dp))/real(N,dp))                      !calculamos integral con la ecuacion

!-----------------------------------------------------------------------------------------------!
! TODO B3: acumule las dos integrales y sus cuadrados para las R corridas.

     s_vm = s_vm + I_vm                                                 !sumatoria para ajuste
     s2_vm = s2_vm + I_vm*I_vm
     s_ae = s_ae + I_ae
     s2_ae = s2_ae + I_ae*I_ae

  end do

!-----------------------------------------------------------------------------------------------!
! TODO B4: calcule media y desviación estándar muestral para cada método.

  media_vm = s_vm/real(R,dp)
  desv_vm = sqrt((real(R, dp) / (real(R, dp) - 1.0_dp)) * (s2_vm / real(R, dp) - media_vm**2))

  media_ae = s_ae/real(R,dp)
  desv_ae = sqrt((real(R, dp) / (real(R, dp) - 1.0_dp)) * (s2_ae / real(R, dp) - media_ae**2))

  !Cuatro números sin encabezado para la lectura automática (guardar archivo)
  open(newunit=u_dat,file='seno.dat',status='replace')
  write(u_dat,'(4es25.16e3)') media_vm,desv_vm,media_ae,desv_ae
  close(u_dat)

  !imprimos los resultados
  print '(a,2f10.5)', 'Valor medio (media, desviacion): ',media_vm,desv_vm
  print '(a,2f10.5)', 'Acierto y error (media, desviacion): ',media_ae,desv_ae

end program taller5_seno
