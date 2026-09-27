
program ajuste_pendulo
use, intrinsic :: iso_fortran_env, only: real64, iostat_end
    implicit none

!Definimos las constantes que usaremos en el código (PARTE 1)                         
    real(real64) , parameter :: PI = 3.14159265358979323846                 !Valor de PI para usar
    integer :: u_in, u_out, ierr, id, n_osc                                 !unidades del archivo e indicadores aceptados
    integer :: n                                                            !Contador de filas/datos válidos
    real(real64) :: l_cm, ang, t_s                                          !Variables leidas de la tabla: Longitud (cm) y Período (s)
    real(real64) :: periodo, x, y, y_pred                                   !Variables transformadas para el ajuste: x = L (m), y = T^2 (s^2)
    real(real64) :: sx, sy, sxx, sxy, denom                                 !definimos las variables (acumuladores) del minimo cuadrados
    real(real64) :: a, b, g, r2                                                !definimos las variables de la pendiente, intercepto, gravedad, r^2
    real(real64) :: y_mean, ss_res, ss_tot                  !definimos variables para el calculo de R^2

! -------------------------------------------------------
! PASO 1: Abrir pendulo_limpio.dat y comprobar que el archivo se abrió correctamente.

open(newunit=u_in, file='pendulo_limpio.dat', status='old', action='read', iostat=ierr)
if (ierr /= 0) error stop "No se pudo abrir el archivo pendulo_limpio.dat"


! -------------------------------------------------------
! PASO 2: Leer un número desconocido de filas mediante un ciclo y iostat.


!iniciamos contadores de filas y acumuladores en cero.
n = 0
sx = 0.0_real64
sy = 0.0_real64
sxx = 0.0_real64
sxy = 0.0_real64


   do
        read(u_in, *, iostat=ierr) id, l_cm, ang, n_osc, t_s              !Leemos los 5 campos del archivo
        
        if (ierr == iostat_end) exit                                      !Salir si se llega al final del archivo leido
        
        if (ierr /= 0) error stop "Fallo al leer los datos del archivo"   !Detenemos si hay un error de formato leido
     

        !PASO 3: Convertir cada longitud de centímetros a metros.
        x = l_cm / 100.0_real64                                                  
        
        !PASO 4: Convertir T y construir X=L and y=T^2
        periodo = t_s / real(n_osc, real64)
        y = periodo**2

        !PASO 5: Acumular N, sx, sy, sxx, sxy
        n = n + 1                                                        
        sx = sx + x
        sy = sy + y
        sxx = sxx + (x * x)
        sxy = sxy + (x * y)

    end do


! -------------------------------------------------------
! PASO 11: Detectar al menos estos errores: archivo inexistente, menos de dos datos y denominador nulo.
! este paso se hace antes ya que se necesita para continuar primero con el calculo de la pendeinte e intercepto y R^2

if (n<2) error stop "Error: Se requieren al menos 2 datos para el ajuste"
denom= real(n, real64)* sxx - sx*sx
if (abs(denom) <= tiny(denom)) error stop 'No se puede calcular la pendiente'


! -------------------------------------------------------
! PASO 6: Calcular la pendiente (a) y el intercepto (b)

a = (real(n, real64) * sxy - sx * sy) / denom              !definimos la pendiente segun las fmormulas
b = (sy - a * sx) / real(n, real64)                        !definimos el intercepto segun las formulas

! -------------------------------------------------------
!PASO 7: Estimar la gravedad
g = (4.0_real64 * (PI**2)) / a                             !definimos la gravedad segun el despeje de la pendiente


! -------------------------------------------------------
! PASO 8: Realizar una segunda lectura del archivo para calcular R^2

rewind(u_in)                                                   !Volvemos al inicio del archivo

    y_mean = sy / real(n, real64)
    ss_res = 0.0_real64
    ss_tot = 0.0_real64

    do
        read(u_in, *, iostat=ierr) id, l_cm, ang, n_osc, t_s                 !revisar todos las variables para evitar errores en la lectura
        if (ierr == iostat_end) exit
        if (ierr /= 0) error stop "Error de lectura en la segunda pasada"

        x = l_cm / 100.0_real64                                              !definimos que x=L segun nuestra transformacion
        periodo = t_s / real(n_osc, real64)                                  !definimos nuevamente que es el periodo
        y = periodo**2                                                       !definimos y=T^ segun nuestra transformacion
        y_pred = a * x + b                                                   !prediccion de la recta

        ss_res = ss_res + (y - y_pred)**2                                    !suma de residuos al cuadrado 
        ss_tot = ss_tot + (y - y_mean)**2                                    !variacion alrededor de la media

    end do

    close(u_in)

    if (ss_tot > 0.0_real64) then
        r2 = 1.0_real64 - (ss_res / ss_tot)
    else
        r2 = 0.0_real64
    end if


! -------------------------------------------------------
! PASO 9: Mostar claramente N, a, b, R^2 y g, incluyendo unidades

print "(A,I0)",           "Numero de datos N  = ", n
print "(A,F12.6,A)",      "Pendiente a        = ", a, " s^2/m"
print "(A,F12.6,A)",      "Intercepto b       = ", b, " s^2"
print "(A,F12.6)",        "Coeficiente R2     = ", r2
print "(A,F12.6,A)",      "Aceleracion g      = ", g, " m/s^2"


! -------------------------------------------------------
! PASO 10: Escribir esos cinco valores, en ese orden y sin encabezado, en resultados_ajuste.dat.

open(newunit=u_out, file='resultados_ajuste.dat', status='replace', action='write')
write(u_out, *) n, a, b, r2, g
close(u_out)


end program ajuste_pendulo
